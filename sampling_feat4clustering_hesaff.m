clear all;
clc;

%% parameter setting
% directory setup
DB = 'MED2013';
database.num_sampled_features = 1e8;
switch DB
case 'MED2013'
    work_dir = fullfile('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/',DB);
end

feature_config = '-hesaff -rootsift -noangle';
feature_name = strrep(feature_config, '-', '');
feature_name = strrep(feature_name, ' ', '_');
database.mat_dir = fullfile(work_dir,[feature_name '_mat']);
database.cluster_dir = fullfile(work_dir,[feature_name '_cluster']);
if ~exist(database.cluster_dir,'dir'),
    mkdir(database.cluster_dir);
end;

clobber = false;
list_mat_file = ['/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/' DB '/meta/', feature_name, '_mat_dir.mat'];
if exist(list_mat_file, 'file')
	load(list_mat_file);
else
	mat_dir = dir(fullfile(database.mat_dir,'*.mat')); 
	mat_dir = {mat_dir(:).name}; 
	save(list_mat_file, 'mat_dir');
end

database.num_mat = length(mat_dir);
kp_length = 5;
desc_length = 128;
data_type = 'single';
database.totalfeatures = 0;
features_per_mat = zeros(1,database.num_mat);

% Get number of all features
fprintf('get total number of image feature for clustering...\n');
range = [1:20000:database.num_mat database.num_mat+1];
rid = 17 % ######### CHANGE HERE FROM 1 --> 24 to multi thread

nfeat_files = dir(['/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/' DB '/meta/',feature_name,'_*']);
if length(nfeat_files) ~= length(range) % check da co du thong tin so luong feat chua
	tic;
	for i = range(rid):range(rid+1)-1
		fprintf('\r%d/%d ',i,database.num_mat);
		mat_feat_file = fullfile(database.mat_dir,mat_dir{i});
		load(mat_feat_file);
		num_mat_frames = size(clip_kp,2);
		for j = 1:num_mat_frames,
			database.totalfeatures = database.totalfeatures + size(clip_kp{j},2);
		end
	end
	fid = fopen(['/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/' DB '/meta/',feature_name,'_',num2str(range(rid)),'_',num2str(database.totalfeatures), '_feats.txt'], 'w');
	fclose(fid);
	database.counting_time = toc;
	fprintf('total feature: %d, %.0f s',database.totalfeatures,database.counting_time);
	return;
end

re = 'vgg_hesaff_rootsift_noangle_\d*_(\d*)_feats.txt';
database.totalfeatures = 0;
for i=1:length(nfeat_files)
	[rematch retok] = regexp(nfeat_files(i).name, re, 'match', 'tokens');
    if ~isempty(retok)
        database.totalfeatures = database.totalfeatures + str2double(retok{1}{1});
    end
end
% sampling ratio 
database.sampling_ratio = database.num_sampled_features/database.totalfeatures;
assert(database.sampling_ratio<=1);

fprintf('Sampling image feature for clustering...');

nsampledfeatures = 0;
sampling_files = dir([database.cluster_dir,'/',feature_name,'*']);
if length(sampling_files) ~= length(range)-1 % check da co du thong tin so luong feat chua
	database.sampleddesc=zeros([desc_length,database.num_sampled_features],data_type);
	% Vet not so feature cho du 100M :D
    for i = range(rid):range(rid+1)-1
        if nsampledfeatures >= database.num_sampled_features
            break;
        end
        mat_name = mat_dir{i};
        mat_feat_file = fullfile(database.mat_dir, mat_name);
        load(mat_feat_file);
        if size(clip_kp,2) < 1
            continue;
        end

        eflag = cell2mat(cellfun(@(x) ~isempty(x), clip_kp, 'UniformOutput', false));

        switch data_type
        case 'uint8'
            sampled_desc = uint8(cell2mat(clip_desc(eflag)));
        case 'single'
            sampled_desc = single(cell2mat(clip_desc(eflag)));
        end
        nsamples = round(size(sampled_desc,2)*database.sampling_ratio);
        nsamples = min(database.num_sampled_features - nsampledfeatures,nsamples);
        if nsamples == 0
            continue;
        end
        rand_list = randperm(size(sampled_desc,2));
        rand_list = rand_list(1:nsamples);
        sampled_desc = sampled_desc(:,rand_list);
        database.sampleddesc(:,nsampledfeatures+1:nsampledfeatures+nsamples)=sampled_desc;
        nsampledfeatures = nsampledfeatures + nsamples;
        fprintf('\r%d/%d: sampled %d    ',i,database.num_mat,nsampledfeatures);
    end
	nsamples = nsampledfeatures;
	tmp_data=database.sampleddesc(:,1:nsamples);
	
	% Save keypoints and descriptor to mat file
	sampled_feat_matfile = fullfile(database.cluster_dir,sprintf('%s%d_%d.mat',feature_name,range(rid),nsamples));
    clear database;
	save(sampled_feat_matfile, 'tmp_data', '-v7.3');
else
	data=zeros([desc_length,database.num_sampled_features],data_type);
	for i=1:length(sampling_files)
		load(fullfile(database.cluster_dir,sampling_files(i).name));
		ntmp_sampling = size(tmp_data,2);
		if ntmp_sampling > database.num_sampled_features-nsampledfeatures
			ntmp_sampling = database.num_sampled_features-nsampledfeatures;
			data(:,nsampledfeatures+1:nsampledfeatures+ntmp_sampling) = tmp_data(:,1:ntmp_sampling);
		else
			data(:,nsampledfeatures+1:nsampledfeatures+ntmp_sampling) = tmp_data;
		end
		nsampledfeatures = nsampledfeatures+ntmp_sampling;
		clear tmp_data
	end
	% neu van con thieu feature
    for i = 1:database.num_mat
        if nsampledfeatures >= database.num_sampled_features
            break;
        end
        mat_name = mat_dir{i};
        mat_feat_file = fullfile(database.mat_dir, mat_name);
        load(mat_feat_file);
        if size(clip_kp,2) < 1
            continue;
        end

        eflag = cell2mat(cellfun(@(x) ~isempty(x), clip_kp, 'UniformOutput', false));

        switch data_type
        case 'uint8'
            sampled_desc = uint8(cell2mat(clip_desc(eflag)));
        case 'single'
            sampled_desc = single(cell2mat(clip_desc(eflag)));
        end
        nsamples = round(size(sampled_desc,2)*database.sampling_ratio);
        nsamples = min(database.num_sampled_features - nsampledfeatures,nsamples);
        if nsamples == 0
            continue;
        end
        rand_list = randperm(size(sampled_desc,2));
        rand_list = rand_list(1:nsamples);
        sampled_desc = sampled_desc(:,rand_list);
        data(:,nsampledfeatures+1:nsampledfeatures+nsamples)=sampled_desc;
        nsampledfeatures = nsampledfeatures + nsamples;
    end
	% Save all keypoints and descriptor
	sampled_feat_hdf5file = fullfile(database.cluster_dir,sprintf('%s%d_128D.hdf5',feature_name,nsampledfeatures));
	hdf5write(sampled_feat_hdf5file,'/data',data);
end
