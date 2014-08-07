function extract_vgg_hesaff_rootsift_noangle(startShotInd, endShotInd)
% example run
% extract_vgg_hesaff_rootsift_noangle('INS2013', 1, 1)

DB = 'INS2013';
switch DB
case 'INS2013'
	lst_shots_file = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/meta/lst_shots.mat';
	db_frame_dir = '/net/per610a/export/das11g/caizhizhu/ins/ins2013/frames_png';
	db_feat_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/vgg_hesaff_rootsift_noangle_mat';
	if ~exist(db_feat_dir, 'dir')
		mkdir(db_feat_dir);
	end
end
renew = false;
detector = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/tool/h_affine.ln';
descriptor = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/tool/compute_descriptors.ln';
run('/net/per610a/export/das11f/ledduy/plsang/nvtiep/libs/vlfeat-0.9.18/toolbox/vl_setup.m');

feature_config = '-vgg -hesaff -rootsift -noangle';
feature_dir = strrep(feature_config, '-', '');
feature_dir = strrep(feature_dir, ' ', '_');
tmp_dir = '/tmp';
sift_dir = fullfile(tmp_dir, feature_dir);
sift_dir = [sift_dir, '_nvtiep'];
if ~exist(sift_dir, 'dir')
	mkdir(sift_dir);
end

% open list shot file
load(lst_shots_file);
test_ids = cellfun(@(x) isempty(strfind(x, 'shot0_')), lst_shots, 'UniformOutput', false);
lst_shots = lst_shots(cell2mat(test_ids));
clear test_ids
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
	
	% Extract feature using h_affine detector and compute_descriptors.ln with sift desc and noangle (VGG)
	for k=1:num_frame
		frame_name = frame_folders{k};
		frame_path = fullfile(shot_frame_dir, frame_name);
		kp_filename = fullfile(sift_dir, [shot_name,'_',frame_name]);
		kp_filename = strrep(kp_filename, 'png', 'hesaff');
		sift_filename = fullfile(sift_dir, [shot_name,'_',frame_name]);
		sift_filename = strrep(sift_filename, 'png', 'txt');
		clip_frame{k}=frame_name(1:end-4);
		% detect keypoint
		cmd = sprintf('%s -hesaff -i %s -o %s', detector, frame_path, kp_filename);
		unix(cmd);
		
		% check number of keypoint
		kp_fid = fopen(kp_filename,'r');
		fgetl(kp_fid);
		nkp = str2double(strtrim(fgetl(kp_fid)));
		fclose(kp_fid);
		if nkp == 0
			%remove temporary file
			unix(['rm ', kp_filename]);
			continue;
		end
		
		% describe keypoint
		while ~exist(sift_filename, 'file')
			cmd = sprintf('%s -sift -noangle -i %s -p1 %s -o1 %s', descriptor, frame_path, kp_filename, sift_filename);
			unix(cmd);
		end
		
		% read to memory
		[clip_kp{k},clip_desc{k}] = vl_ubcread(sift_filename, 'format', 'oxford');
		feat_len = size(clip_desc{k},1);

		if ~isempty(strfind(feature_config,'rootsift'))
			sift = double(clip_desc{k});
			clip_desc{k} = single(sqrt(sift./repmat(sum(sift), feat_len, 1)));
		end
		% remove temporary file
		unix(['rm ', sift_filename]);

		% remove temporary file
		unix(['rm ', kp_filename]);
	end
	% save to .mat file
	save(shot_feature_file, 'clip_kp', 'clip_desc', 'clip_frame', '-v7.3'); 
end
%matlabpool('close');
%unix(['rm -r ', sift_dir]);
quit;
end
