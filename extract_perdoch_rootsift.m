function extract_perdoch_rootsift(DB, startShotInd, endShotInd)
% example run
% extract_perdoch_rootsift('INS2013', 1, 1)

renew = false;
addpath('/net/per610a/export/das11f/ledduy/plsang/nvtiep/funcs/perdoch_hesaff');
feature_config = '-perdoch -hesaff -rootsift';

switch DB
case 'INS2013'
	lst_shots_file = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/meta/lst_shots.mat';
	db_frame_dir = '/net/per610a/export/das11g/caizhizhu/ins/ins2013/frames_png';
	db_feat_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/perdoch_hesaff_rootsift_mat';
end

% open list shot file
load(lst_shots_file);
nshot = length(lst_shots);

%matlabpool('open', 8);	
for i=startShotInd:endShotInd
	fprintf('\r %d - %d - %d', i, startShotInd, endShotInd);
	shot_name = lst_shots{i};
	shot_frame_dir = fullfile(db_frame_dir, shot_name);
	shot_feature_file = fullfile(db_feat_dir, [shot_name,'.mat']);
	if exist(shot_feature_file, 'file') && ~renew
		continue;
	end
	% Load all frames of a shot
	fid = fopen(fullfile(shot_frame_dir,'frames.txt'));
	frame_folders = textscan(fid, '%s');
	fclose(fid);
	frame_folders = frame_folders{1};
	
    % Number of frames
    num_frame = length(frame_folders);
	clip_kp = cell(1,num_frame);
	clip_desc = cell(1,num_frame);
	clip_frame = cell(num_frame,1);
	
	% Extract feature using perdoch hesaff rootsift
	for k=1:num_frame
		frame_name = frame_folders{k};
		frame_path = fullfile(shot_frame_dir, frame_folders{k});
		clip_frame{k}=frame_name(1:end-4);
		[clip_kp{k},clip_desc{k}] = mxhesaff(frame_path,~isempty(strfind(feature_config,'root')),false);
	end
	% save to .mat file
	save(shot_feature_file, 'clip_kp', 'clip_desc', 'clip_frame', '-v7.3'); 
end
%matlabpool('close');
quit;
end
