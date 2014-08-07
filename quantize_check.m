function quantize_check(sID, eID)

%% Set perameters
DB = 'INS2013';
disp(DB)
work_dir = fullfile('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS', DB);
grid_mode = true;

database.comp_sim= struct('query_obj','fg+bg_0.1','feat_detr','-perdoch -hesaff', 'feat_desc', '-rootsift',...
  'clustering','akmeans','K',1000000,'num_samps',100000000,'iter',50,...
  'build_params',struct('algorithm', 'kdtree','trees', 8, 'checks', 800, 'cores', 20),...
  'video_sampling',1,'frame_sampling',1,'knn',3,'delta_sqr',6250,'db_agg','avg_pooling',...
  'vocab','full','trim','notrim','freq','clip','weight','idf','norm','l1',...
  'query_knn',3,'query_delta_sqr',6250,'query_num',-1,'query_agg','avg_pooling','dist','l2asym_ivf');

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
if ~exist(database.subbow_dir,'dir'),
    mkdir(database.subbow_dir);
end

% Quantize database
load(fullfile(work_dir,'meta','lst_shots.mat'));
database.quant_dir = fullfile(database.build_dir, [quantize_name,'_sub_quant']);
if ~exist(database.quant_dir,'dir'),
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
end

disp('Building bow files separately');
% Clear remained data
if exist('kdtree', 'var')
	flann_free_index(kdtree);
end
clear test_ids bins sqrdists frame_desc clip_desc dataset search_params frame_desc

% Log
logfile=fopen(sprintf('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/log/quantize_%d.txt', grid_mode),'a');
fprintf(logfile, '%d - %d: Bat dau checking \n', sID, eID);
fclose(logfile);
% Checking bag of word
for clip_id = sID:eID
	fprintf('\r%d(%d~%d) -%s', clip_id,sID,eID, lst_shots{clip_id});
	quant_file = fullfile(database.quant_dir, [lst_shots{clip_id},'.mat']);
	subbow_file = fullfile(database.subbow_dir, [lst_shots{clip_id},'.mat']);
	caizhi_subbow_file = fullfile('/net/per610a/export/das11g/caizhizhu/ins/ins2013/perdoch_hesaff_rootsift_cluster/akmeans_1000000_100000000_50/kdtree_8_800/v1_f1_3_0.0125/sub_bow', [lst_shots{clip_id},'.mat']);
	
	qfile = dir(quant_file);
	sfile = dir(subbow_file);
	% Check error bow file
	if length(qfile)==1 && length(sfile)==1	
		try
			load(subbow_file);
			load(quant_file);
			%fprintf('\r%d/%d', clip_id,num_clip);
			% Compare with Caizhi data
			my_frame_bow = frame_bow;
			my_frame_freq= frame_freq;
			load(caizhi_subbow_file);
			% Check
			if sum(size(my_frame_bow)-size(frame_bow)) ~= 0 || sum(size(my_frame_freq)-size(frame_freq)) ~= 0 || abs(sum(sum(my_frame_bow-frame_bow)))>0.1 || abs(sum(sum(my_frame_bow-frame_bow)))>=0.1
				logfile=fopen(sprintf('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/log/quantize_%d.txt', grid_mode),'a');
				fprintf(logfile, '%d - %d: Loi o line 125 -> Shot:%s - Khac anh Caizhi\n', sID, eID, lst_shots{clip_id});
				fclose(logfile);
			end
		catch err
			%unix(['rm ', subbow_file]);
			logfile=fopen(sprintf('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/log/quantize_%d.txt', grid_mode),'a');
			fprintf(logfile, '%d - %d: Loi o line 131 -> Shot:%s - %s - Khong load duoc\n', sID, eID, lst_shots{clip_id}, err.identifier);
			fclose(logfile);
		end
	else
		logfile=fopen(sprintf('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/log/quantize_%d.txt', grid_mode),'a');
		fprintf(logfile, '%d - %d: Loi o line 136 -> Shot:%s - Khong ton tai image\n', sID, eID, lst_shots{clip_id});
		fclose(logfile);
	end
	%% Build bow
	clear bins sqrdists frame_bow frame_freq my_frame_bow my_frame_freq
end

% Log
logfile=fopen(sprintf('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/log/quantize_%d.txt', grid_mode),'a');
fprintf(logfile, '%d - %d: Build bow done!\n', sID, eID);
fclose(logfile);

fprintf('\nFinished\n');
quit;