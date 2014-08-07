
%% Load list shot
run('/net/per610a/export/das11f/ledduy/plsang/nvtiep/libs/vlfeat-0.9.18/toolbox/vl_setup.m');
addpath(genpath('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/code/funcs'));
addpath(genpath('/net/per610a/export/das11f/ledduy/plsang/nvtiep/funcs'));

%% Set perameters
DB = 'INS2013';
disp(DB)
clobber = false
save_raw_bow = true
save_max_avg_bow = true
idf_l1_norm = true
work_dir = fullfile('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS', DB);

%database.comp_sim= struct('query_obj','fg+bg_0.1','feat_detr','-perdoch -hesaff', 'feat_desc', '-rootsift',...
%  'clustering','akmeans','K',1000000,'num_samps',100000000,'iter',50,...
%  'build_params',struct('algorithm', 'kdtree','trees', 8, 'checks', 800, 'cores', 20),...
%  'video_sampling',1,'frame_sampling',1,'knn',3,'delta_sqr',6250,'db_agg','avg_pooling',...
%  'vocab','full','trim','notrim','freq','clip','weight','idf','norm','l1',...
%  'query_knn',3,'query_delta_sqr',6250,'query_num',-1,'query_agg','dist_avg','dist','l2asym_ivf');

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
if ~exist(database.subbow_dir,'dir'),
    mkdir(database.subbow_dir);
end

%% Quantize database
load(fullfile(work_dir,'meta','lst_shots.mat'));
database.quant_dir = fullfile(database.build_dir,quantize_name);
if ~exist(database.quant_dir,'dir'),
	mkdir(database.quant_dir);
end

%filter out those sample video
test_ids = cellfun(@(x) isempty(strfind(x, 'shot0_')), lst_shots, 'UniformOutput', false);
lst_shots = lst_shots(cell2mat(test_ids));
disp('Building quant files separately');
num_clip = length(lst_shots);

%% Quantize and build bag-of-word
frame_sampling = database.comp_sim.frame_sampling;
avg_big_bow_file = fullfile(database.bow_dir,'avg_pooling_raw_bow.mat');
%max_big_bow_file = fullfile(database.bow_dir,'max_pooling_raw_bow.mat');
raw_bow_file = fullfile(database.bow_dir,'raw_bow.mat');
big_bow_info_file = fullfile(database.bow_dir,'raw_bow_info.mat');

if exist(big_bow_info_file)
	disp('big bow files exist!');
	%return
end

% Load dictionary to get number of code book
cluster_filename = dir(fullfile(database.cluster_dir,'Cluster*.hdf5'));
assert(length(cluster_filename) == 1);
database.cluster_filename = cluster_filename(1).name;
centers = hdf5read(fullfile(database.cluster_dir, database.cluster_filename),'/clusters');
[feat_len,hist_len] = size(centers);

%% Load bow files
disp('Load bow files into a big one');
list_term_freq=struct('occu',zeros(hist_len+1,1),'frame',zeros(hist_len+1,1),'clip',zeros(hist_len+1,1));
list_clip_frame_num = zeros(num_clip,1);

list_id2clip_lut = lst_shots; %cellfun(@(x) x(1:end-4), lst_shots, 'UniformOutput',false);

% Clear unused data
clear test_ids centers

%% begin merge data
sID = 1;
eID = num_clip;

if save_max_avg_bow
	%clear list_max_pooling_bow
	%list_max_pooling_bow = sparse(hist_len,eID-sID+1);
	clear list_avg_pooling_bow
	list_avg_pooling_bow = sparse(hist_len,eID-sID+1);
end

if save_raw_bow
	clear list_frame_bow
	list_frame_bow = cell(1,eID-sID+1);
end

%% Load all part of raw bow
file_name = 'raw_bow.mat';
file_path = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/perdoch_hesaff_rootsift_cluster/akmeans_1000000_100000000_50/kdtree_8_800/v1_f1_3_0.0125/';
files = dir([file_path,file_name,'*']);

% Sort files
nfile = length(files);

for i=1:nfile-1
	[remat, retok] = regexp(files(i).name, 'from_(\d*).', 'match', 'tokens');
	idi = str2double(retok{end}{1});
	for j=i+1:nfile
		[remat, retok] = regexp(files(j).name, 'from_(\d*).', 'match', 'tokens');
		idj = str2double(retok{end}{1});
		if idi>idj
			tmp = files(i).name;
			files(i).name=files(j).name;
			files(j).name=tmp;
			idi=idj;
		end
	end
end

% Load

clip_id = 0;
for i=1:nfile
	clear list_frame_bow list_clip_freq
	load(fullfile(file_path,files(i).name));
	num_new_shot = size(list_frame_bow,2);
	
	for j = 1:num_new_shot
		clip_id = clip_id+1;
		fprintf('\r%d(%d~%d)', clip_id, i, nfile);
		num_frame = size(list_frame_bow{j},2);
		if num_frame == 0
			num_frame
		end
		clip_freq = sum(list_clip_freq{j},2);
		
		list_avg_pooling_bow(:,clip_id) = mean(list_frame_bow{j},2);
		
		list_clip_frame_num(clip_id) = num_frame;
		list_term_freq.occu = list_term_freq.occu+[clip_freq;sum(clip_freq)];
		list_term_freq.clip = list_term_freq.clip+[double(clip_freq>0);1];
		list_term_freq.frame = list_term_freq.frame+[sum(list_frame_bow{j}>0,2);size(list_frame_bow{j},2)];
	end
end

% make sure your memory can handle
fprintf('\n save big bow info');
save(big_bow_info_file, 'list_term_freq','list_clip_frame_num','list_id2clip_lut','-v6');   

fprintf('\n save big avg bow ..');
save(avg_big_bow_file, 'list_avg_pooling_bow','-v7.3');   

if idf_l1_norm 
	load(big_bow_info_file);
	term_freq = list_term_freq.clip;
	weight = get_wei(term_freq,'idf');
	
	%% Build for avg first
	db_bow = list_avg_pooling_bow;
	db_lut=list_id2clip_lut;
	clip_frame_num = list_clip_frame_num;
	clear term_freq list_clip_frame_num list_avg_pooling_bow list_term_freq
	
	wei_nonorm_bow_file = fullfile(database.bow_dir,'bow_clip_full_notrim_clip_idf_nonorm_avg_pooling.mat');
	wei_norm_bow_file = fullfile(database.bow_dir,'bow_clip_full_notrim_clip_idf_l1_avg_pooling.mat');
	fprintf('\n weighting ...');
	tic;
	for i=1:size(db_bow,2)
		db_bow(:,i) = db_bow(:,i).*weight;
	end
	fprintf('%.0fs\n', toc);
	save(wei_nonorm_bow_file,'db_bow','db_lut','weight','clip_frame_num','-v7.3');
	fprintf('\n Normalizing ...');
	tic;
	for i=1:size(db_bow,2)
		bow_norm = sum(db_bow(:,i))+eps;
		db_bow(:,i) = db_bow(:,i)./bow_norm;
	end 
	fprintf('%.0fs\n', toc);
	save(wei_norm_bow_file,'db_bow','db_lut','weight','clip_frame_num','-v7.3');
end
disp('finished');

