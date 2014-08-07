function extract_hesaffine_rootsift_noangle(startShotInd, endShotInd)
% example run
% extract_hesaff_rootsift_noangle('INS2013', 1, 1)
DB = 'MED2013';
switch DB
case 'MED2013'
	lst_shots_file = '/net/per610a/export/das11f/plsang/trecvidmed13/metadata/medmd_2014.mat';
	db_frame_dir = '/net/per610a/export/das11f/plsang/trecvidmed13/keyframes';
	db_feat_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/MED2013/hesaff_rootsift_noangle_mat';
	if ~exist(db_feat_dir, 'dir')
		mkdir(db_feat_dir);
	end
	% open list shot file
	load(lst_shots_file, 'MEDMD');
	lst_shots = MEDMD.EventBG.default.clips;
	nshot = length(lst_shots);
end

exe = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/tool/compute_descriptors_64bit.ln';
run('/net/per610a/export/das11f/ledduy/plsang/nvtiep/libs/vlfeat-0.9.18/toolbox/vl_setup.m');

feature_config = '-hesaff -rootsift -noangle';
feature_dir = strrep(feature_config, '-', '');
feature_dir = strrep(feature_dir, ' ', '_');
tmp_dir = '/tmp';
sift_dir = fullfile(tmp_dir, feature_dir);
if ~exist(sift_dir, 'dir')
	mkdir(sift_dir);
end


%matlabpool('open', 8);
for i=startShotInd:endShotInd
	fprintf('\r %d - %d - %d', i, startShotInd, endShotInd);
	shot_id = lst_shots{i};
	shot_file = MEDMD.lookup.(shot_id);
	shot_frame_dir = sprintf('%s/%s', db_frame_dir, shot_file(1:end-4));
	shot_feature_file = fullfile(db_feat_dir, [shot_id,'.mat']);
	if exist(shot_feature_file, 'file')
		continue;
	end
	% Load all frames of a shot
	frame_folders = dir([shot_frame_dir '/*.jpg']);
	frame_folders = {frame_folders.name};
    % Number of frames
    num_frame = length(frame_folders);
	clip_kp = cell(1,num_frame);
	clip_desc = cell(1,num_frame);
	clip_frame = cell(num_frame,1);
	
	% Extract feature using compute_descriptors_64bit.ln hesaff rootsift noangle
	for k=1:num_frame
		frame_name = frame_folders{k};
		frame_path = fullfile(shot_frame_dir, frame_name);
		sift_filename = fullfile(sift_dir, [shot_id,'_',frame_name]);
		sift_filename = strrep(sift_filename, 'jpg', 'txt');
		clip_frame{k}=frame_name(1:end-4);
		
		cmd = sprintf('%s %s -i %s -o1 %s', exe, strrep(feature_config,'root',''), frame_path, sift_filename);
		unix(cmd);
		if exist(sift_filename, 'file')
			[clip_kp{k},clip_desc{k}] = vl_ubcread(sift_filename, 'format', 'oxford');
			feat_len = size(clip_desc{k},1);

			if ~isempty(strfind(feature_config,'rootsift'))
				sift = double(clip_desc{k});
				clip_desc{k} = single(sqrt(sift./repmat(sum(sift), feat_len, 1)));
			end
		end
	end
	% save to .mat file
	save(shot_feature_file, 'clip_kp', 'clip_desc', 'clip_frame', '-v7.3');
	fileattrib(shot_feature_file, '+w', 'a');
	% remove temporary file
	unix(['rm ', sift_filename]);
end
%matlabpool('close');
unix(['rm -r ', sift_dir]);
quit;
end
