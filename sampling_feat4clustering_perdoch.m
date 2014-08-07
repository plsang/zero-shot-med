%clear all;
clc;

%% parameter setting
% directory setup
DB = 'INS2013';
database.num_sampled_features = 1e8;
switch DB
case 'INS2013'
    work_dir = fullfile('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/',DB);
end

feature_config = '-perdoch -hesaff -rootsift';
feature_name = strrep(feature_config, '-', '');
feature_name = strrep(feature_name, ' ', '_');
database.mat_dir = fullfile(work_dir,[feature_name '_mat']);
database.cluster_dir = fullfile(work_dir,[feature_name '_cluster']);
if ~exist(database.cluster_dir,'dir'),
    mkdir(database.cluster_dir);
end;

clobber = false;
mat_dir = dir(fullfile(database.mat_dir,'*.mat'));
mat_dir = {mat_dir(:).name};
database.num_mat = length(mat_dir);
kp_length = 5;
desc_length = 128;
data_type = 'single';
database.totalfeatures = 0;
features_per_mat = zeros(1,database.num_mat);

% Get number of all features
tic;
fprintf('get total number of image feature for clustering...\n');
for i = 1:database.num_mat,
	fprintf('\r%d/%d ',i,database.num_mat);
	mat_feat_file = fullfile(database.mat_dir,mat_dir{i});
	load(mat_feat_file);
	num_mat_frames = size(clip_kp,2);
	for j = 1:num_mat_frames,
		database.totalfeatures = database.totalfeatures + size(clip_kp{j},2);
	end
end
fid = fopen(['/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/meta/nfeat_is_', num2str(database.totalfeatures), '.txt'], 'w');
fclose(fid);
database.counting_time = toc;
fprintf('total feature: %d, %.0f s',database.totalfeatures,database.counting_time);

% sampling ratio 
database.sampling_ratio = database.num_sampled_features/database.totalfeatures;
assert(database.sampling_ratio<=1);

fprintf('Sampling image feature for clustering...');
database.sampledkp=zeros([kp_length,database.num_sampled_features],'double');
database.sampleddesc=zeros([desc_length,database.num_sampled_features],data_type);
nsampledfeatures = 0;
mat_id = 0;
for j = 1:2
	% Vet not so feature cho du 100M :D
    for i = database.num_mat:-1:1
        if nsampledfeatures >= database.num_sampled_features
            break;
        end
        mat_name = mat_dir{i};
        mat_id = mat_id + 1;
        mat_feat_file = fullfile(database.mat_dir, mat_name);
        load(mat_feat_file);
        if size(clip_kp,2) < 1
            continue;
        end

        eflag = cell2mat(cellfun(@(x) ~isempty(x), clip_kp, 'UniformOutput', false));

        sampled_kp = cell2mat(clip_kp(eflag));
        switch data_type
        case 'uint8'
            sampled_desc = uint8(cell2mat(clip_desc(eflag)));
        case 'single'
            sampled_desc = single(cell2mat(clip_desc(eflag)));
        end
        nsamples = round(size(sampled_kp,2)*database.sampling_ratio);
        nsamples = min(database.num_sampled_features - nsampledfeatures,nsamples);
        if nsamples == 0
            continue;
        end
        rand_list = randperm(size(sampled_kp,2));
        rand_list = rand_list(1:nsamples);
        sampled_kp = sampled_kp(:,rand_list);
        sampled_desc = sampled_desc(:,rand_list);
        database.sampledkp(:,nsampledfeatures+1:nsampledfeatures+nsamples)=sampled_kp;
        database.sampleddesc(:,nsampledfeatures+1:nsampledfeatures+nsamples)=sampled_desc;
        nsampledfeatures = nsampledfeatures + nsamples;
        fprintf('\r%d/%d: sampled %d    ',i,database.num_mat,nsampledfeatures);
    end;
end
database.sampling_time = toc;
fprintf('%.0f',database.sampling_time);

nsamples = nsampledfeatures;
data=database.sampleddesc;
kp = single(database.sampledkp);
clear database;

% Save keypoints and descriptor
sampled_feat_matfile = fullfile(database.cluster_dir,sprintf('%s%d.mat',feature_name,nsamples));
save(sampled_feat_matfile, 'kp','-v7.3');
sampled_feat_hdf5file = strrep(sampled_feat_matfile,'.mat','_128D.hdf5');
hdf5write(sampled_feat_hdf5file,'/data',data);

