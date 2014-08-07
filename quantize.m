function quantize(sID, eID)

%% Load list shot
run('/net/per610a/export/das11f/ledduy/plsang/nvtiep/libs/vlfeat-0.9.18/toolbox/vl_setup.m');
addpath(genpath('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/code/funcs'));
addpath(genpath('/net/per610a/export/das11f/ledduy/plsang/nvtiep/funcs'));

%% Set perameters
DB = 'INS2013';
disp(DB)
work_dir = fullfile('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS', DB);
build_quant = true;
build_bow = true;
grid_mode = false;

database.comp_sim= struct('query_obj','fg+bg_0.1','feat_detr','-vgg -hesaff', 'feat_desc', '-rootsift -noangle',...
  'clustering','akmeans','K',1000000,'num_samps',100000000,'iter',50,...
  'build_params',struct('algorithm', 'kdtree','trees', 8, 'checks', 800, 'cores', 20),...
  'video_sampling',1,'frame_sampling',1,'knn',3,'delta_sqr',6250,'db_agg','avg_pooling',...
  'vocab','full','trim','notrim','freq','clip','weight','idf','norm','l1',...
  'query_knn',3,'query_delta_sqr',6250,'query_num',-1,'query_agg','avg_pooling','dist','l2asym_ivf');

database.comp_sim= struct('query_obj','fg+bg_0.1','feat_detr','-hesaff', 'feat_desc', '-rootsift -noangle',...
  'clustering','akmeans','K',1000000,'num_samps',100000000,'iter',50,...
  'build_params',struct('algorithm', 'kdtree','trees', 8, 'checks', 800, 'cores', 10),...
  'video_sampling',1,'frame_sampling',1,'knn',1,'delta_sqr',6250,'db_agg','avg_pooling',...
  'vocab','full','trim','notrim','freq','clip','weight','idf','norm','nonorm',...
  'query_knn',3,'query_delta_sqr',6250,'query_num',-1,'query_agg','avg_pooling','dist','l2_ivf');

quant_struct = struct('quantize','kdtree','build_params',database.comp_sim.build_params,'knn',database.comp_sim.knn,'delta_sqr',database.comp_sim.delta_sqr);


% feature name
feature_detr = database.comp_sim.feat_detr;
feature_desc = database.comp_sim.feat_desc;
feature_config = sprintf('%s %s', feature_detr, feature_desc);
feature_name = strrep(feature_config, '-','');
feature_name = strrep(feature_name, ' ','_');
database.feat_mat_dir = sprintf('%s/%s_mat',work_dir,feature_name);

if quant_struct.knn>1 && quant_struct.delta_sqr~=-1
    if ~isempty(strfind(feature_name, 'root'))
        quant_struct.delta_sqr=quant_struct.delta_sqr/5e5;
    elseif ~isempty(strfind(feature_name, 'color'))
        quant_struct.delta_sqr=quant_struct.delta_sqr*2;
    elseif ~isempty(strfind(feature_name, 'mom'))
        quant_struct.delta_sqr=quant_struct.delta_sqr/1e3;
    end
end

% Clustering name
if ~isempty(strfind(database.comp_sim.clustering,'akmeans'))
	clustering_name = sprintf('%s_%d_%d_%d',database.comp_sim.clustering,...
		database.comp_sim.K,database.comp_sim.num_samps,database.comp_sim.iter); 
end
database.cluster_dir = fullfile(work_dir,[feature_name,'_cluster'],clustering_name);

% Quantize name
quantize_name = sprintf('v%d_f%d_%d', database.comp_sim.video_sampling, database.comp_sim.frame_sampling, quant_struct.knn);
bow_name = quantize_name;
if ~strcmp(quant_struct.quantize,'kdtree') 
    build_name = sprintf('%s',quant_struct.quantize); 
else
    build_name = sprintf('%s_%d_%d',quant_struct.build_params.algorithm,quant_struct.build_params.trees,quant_struct.build_params.checks); 
    if quant_struct.knn>1
        bow_name = sprintf('%s_%g', quantize_name,quant_struct.delta_sqr);
    end
end
database.build_dir = fullfile(database.cluster_dir,build_name);
database.bow_dir = fullfile(database.build_dir,bow_name);
database.subbow_dir = fullfile(database.bow_dir,'sub_bow');
if ~exist(database.subbow_dir,'dir')
    mkdir(database.subbow_dir);
end

% Quantize database
load(fullfile(work_dir,'meta','lst_shots.mat'));
database.quant_dir = fullfile(database.build_dir, [quantize_name,'_sub_quant']);
if ~exist(database.quant_dir,'dir')
	mkdir(database.quant_dir);
end

%filter out those sample video
test_ids = cellfun(@(x) isempty(strfind(x, 'shot0_')), lst_shots, 'UniformOutput', false);
lst_shots = lst_shots(cell2mat(test_ids));
num_clip = length(lst_shots);
frame_sampling = database.comp_sim.frame_sampling;

% Load codebook
if strcmp(quant_struct.quantize, 'kdtree') 
	tic;
	disp('Loading codebook file ...');
	cluster_filename = dir(fullfile(database.cluster_dir,'Cluster*.hdf5'));
	assert(length(cluster_filename) == 1);
	database.cluster_filename = cluster_filename(1).name;
	avg_big_bow_file = fullfile(database.bow_dir,'avg_pooling_raw_bow.mat');
	%max_big_bow_file = fullfile(database.bow_dir,'max_pooling_raw_bow.mat');
	raw_bow_file = fullfile(database.bow_dir,'raw_bow.mat');
	big_bow_info_file = fullfile(database.bow_dir,'raw_bow_info.mat');
	if exist(big_bow_info_file,'file')~=0
		disp('big bow files exist');
		%return;
	end
	centers = hdf5read(fullfile(database.cluster_dir, database.cluster_filename),'/clusters');
	dataset = single(centers);
	[feat_len,hist_len] = size(dataset);
	fprintf('Deduced cluster center info %d %d...\n', feat_len, hist_len);
end

if build_quant
	%% Build flann
	if ~isempty(strfind(quant_struct.quantize, 'kdtree'))
		kdtree_filename = fullfile(database.build_dir,'flann_kdtree.bin');
		kdsearch_filename = fullfile(database.build_dir,'flann_kdtree_search.mat');
		
		if exist(kdtree_filename,'file')
			fprintf('Loading kdtree ...');
			tic;
			kdtree = flann_load_index(kdtree_filename,dataset);
			load(kdsearch_filename);
			search_params.cores = quant_struct.build_params.cores;
			kdtree_time.load = toc;
		else
			tic;
			fprintf('Building kdtree ...');
			[kdtree,search_params,speedup] = flann_build_index(dataset,quant_struct.build_params); 
			kdtree_time.speedup = speedup;
			kdtree_time.build = toc;
			fprintf('%.0f \n',kdtree_time.build);
			fprintf('save kdtree ...');
			tic;
			flann_save_index(kdtree,kdtree_filename);
			save(kdsearch_filename,'search_params');
			kdtree_time.save = toc;
		end
	end
end

clear test_ids centers

% Quantize first
if build_quant
	disp('Building quant files separately');
	% Log
	logfile=fopen(sprintf('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/log/quantize_%d.txt', grid_mode),'a');
	fprintf(logfile, '%d - %d: Bat dau quantizing\n', sID, eID);
	fclose(logfile);
	
	for clip_id = sID:eID
		fprintf('\r%d(%d~%d) - %s', clip_id,sID,eID, lst_shots{clip_id});
		quant_file = fullfile(database.quant_dir, [lst_shots{clip_id},'.mat']);
		
		% Check consistency
		if exist(quant_file,'file')
			try
				%load(quant_file);
				%fprintf('\r%d/%d', clip_id,num_clip);
				continue;
			catch err
				%unix(['rm', quant_file]);
				logfile=fopen(sprintf('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/log/quantize_%d.txt', grid_mode),'a');
				fprintf(logfile, '%d - %d: Loi o line 153 -> Shot:%s - %s\n', sID, eID, lst_shots{clip_id}, err.identifier);
				fclose(logfile);
			end
		end
		clip_feat_file = fullfile(database.feat_mat_dir,[lst_shots{clip_id},'.mat']);
		if ~exist(clip_feat_file,'file')
			disp([clip_feat_file ' does not exist!']);
			continue;
		end

		%load feature
		try
			load(clip_feat_file, 'clip_desc');

			%quantize feature
			num_frame = length(clip_desc);
			selected_frame_id = 1:frame_sampling:num_frame;
			selected_frame_num = length(selected_frame_id);
			bins = cell(1,num_frame);
			sqrdists = cell(1,num_frame);
			for id = 1:selected_frame_num
				frame_id = selected_frame_id(id);
				if isempty(clip_desc{frame_id})
					continue;
				end
				frame_desc = clip_desc{frame_id}(1:feat_len,:);
				[bins{frame_id},sqrdists{frame_id}] = flann_search(kdtree,single(frame_desc),quant_struct.knn, search_params);
			end
			% Save quant file
			save(quant_file, 'bins','sqrdists');
		catch err
			logfile=fopen(sprintf('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/log/quantize_%d.txt', grid_mode),'a');
			fprintf(logfile, '%d - %d: Loi o line 186 -> Shot:%s NFrame=%d - %s\n', sID, eID, lst_shots{clip_id}, selected_frame_num, err.identifier);
			fclose(logfile);
		end
		clear bins sqrdists	frame_desc clip_desc
	end
	
	% Log
	logfile=fopen(sprintf('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/log/quantize_%d.txt', grid_mode),'a');
	fprintf(logfile, '%d - %d: Quantizing done!\n', sID, eID);
	fclose(logfile);
end

if build_bow
	disp('Building bow files separately');
	% Clear remained data
	if exist('kdtree', 'var')
		flann_free_index(kdtree);
	end
	clear bins sqrdists	frame_desc clip_desc dataset search_params frame_desc
	
	% Log
	logfile=fopen(sprintf('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/log/quantize_%d.txt', grid_mode),'a');
	fprintf(logfile, '%d - %d: Bat dau build bow\n', sID, eID);
	fclose(logfile);
	% Build Bag of Word
	for clip_id = sID:eID
		fprintf('\r%d(%d~%d) -%s', clip_id,sID,eID, lst_shots{clip_id});
		quant_file = fullfile(database.quant_dir, [lst_shots{clip_id},'.mat']);
		subbow_file = fullfile(database.subbow_dir, [lst_shots{clip_id},'.mat']);
		
		% Check error bow file
		if exist(subbow_file,'file')
			try
				%load(subbow_file);
				%fprintf('\r%d/%d', clip_id,num_clip);
				continue;
			catch err
				%unix(['rm ', subbow_file]);
				logfile=fopen(sprintf('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/log/quantize_%d.txt', grid_mode),'a');
				fprintf(logfile, '%d - %d: Loi o line 224 -> Shot:%s - %s\n', sID, eID, lst_shots{clip_id}, err.identifier);
				fclose(logfile);
			end
		end
		% in grid_mode, the file size should be small
		if grid_mode
			qfile = dir(quant_file);
			if length(qfile)~=1 || qfile(1).bytes > 800000
				continue;
			end
		end
		%% Build bow
		try
			load(quant_file);
			num_frame = length(bins);
			selected_frame_id = 1:frame_sampling:num_frame;
			selected_frame_num = length(selected_frame_id);
			% frame_bow to store bow of each frame
			frame_bow = zeros(hist_len,selected_frame_num);
			if quant_struct.knn>1 && quant_struct.delta_sqr ~= -1 
				frame_freq = zeros(hist_len,selected_frame_num);
			end
			for id = 1:selected_frame_num
				frame_id = selected_frame_id(id);
				if isempty(bins{frame_id})
					continue;
				end
				bin = reshape(bins{frame_id}(1:quant_struct.knn,:),1,[]);
				if quant_struct.knn>1 && quant_struct.delta_sqr ~= -1 
					sqrdist = sqrdists{frame_id}(1:quant_struct.knn,:);
					weis = exp(-sqrdist./(2*quant_struct.delta_sqr));
					weis = weis./repmat(sum(weis,1),size(weis,1),1);  % philbin, Lost in Quantization
					weis = reshape(weis,1,[]);
					frame_freq(:,id) = vl_binsum(frame_freq(:,id),double(ones(size(bin))),double(bin)); 
				else
					weis = ones(size(bin));
				end
				frame_bow(:,id) = vl_binsum(frame_bow(:,id),double(weis),double(bin)); 
			end 
			frame_bow = sparse(frame_bow);
			
			% Save quant_file and subbow_file
			if ~exist('frame_freq','var')
				save(subbow_file, 'frame_bow');
			else
				frame_freq = sparse(frame_freq);
				save(subbow_file, 'frame_bow','frame_freq');
			end
		catch err
			logfile=fopen(sprintf('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/log/quantize_%d.txt', grid_mode),'a');
			fprintf(logfile, '%d - %d: Loi o line 275 -> Shot:%s NFrame=%d - %s\n', sID, eID, lst_shots{clip_id}, selected_frame_num, err.identifier);
			fclose(logfile);
		end

		clear bins sqrdists frame_bow frame_freq weis bin bins sqrdists
	end
	
	% Log
	logfile=fopen(sprintf('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/log/quantize_%d.txt', grid_mode),'a');
	fprintf(logfile, '%d - %d: Build bow done!\n', sID, eID);
	fclose(logfile);
end
fprintf('\nFinished\n');
quit;