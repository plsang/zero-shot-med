function rerank_using_DPM_and_RANSAC(start_video_id, end_video_id)
start_query_id = 9069;
end_query_id = 9098;
is_compute_RANSAC = 0;
debug_mode = false;
is_using_RANSAC = 1;
if nargin == 0
	start_video_id = 511;
	end_video_id = 511;
	start_query_id = 9071;
	end_query_id = 9071;
	debug_shot = 'shot121_200';	% truong hop nay sai, check tai sao fg_locs tang nhanh sau 2 vong lap
	
	%start_video_id = 676;
	%end_video_id = 676;
	%start_query_id = 9071;
	%end_query_id = 9071;
	%debug_shot = 'shot163_376';
	
	%start_video_id = 450;
	%end_video_id = 450;
	%start_query_id = 9069;
	%end_query_id = 9069;
	%debug_shot = 'shot105_834';
	
	%start_video_id = 437;
	%end_video_id = 437;
	%start_query_id = 9069;
	%end_query_id = 9069;
	%debug_shot = 'shot101_1124';
	
	%start_video_id = 156;
	%end_video_id = 156;
	%start_query_id = 9069;
	%end_query_id = 9069;
	%debug_shot = 'shot37_262';
	
	% Push UP
	start_video_id = 876;
	end_video_id = 876;
	start_query_id = 9069;
	end_query_id = 9069;
	debug_shot = 'shot214_1709';
	
	%start_video_id = 446;
	%end_video_id = 446;
	%start_query_id = 9069;
	%end_query_id = 9069;
	%debug_shot = 'shot104_1334';
	
	%start_video_id = 251;
	%end_video_id = 251;
	%start_query_id = 9069;
	%end_query_id = 9069;
	%debug_shot = 'shot57_1251';

	%start_video_id = 284;
	%end_video_id = 284;
	%start_query_id = 9069;
	%end_query_id = 9069;
	%debug_shot = 'shot65_618';
	
	%start_video_id = 615;
	%end_video_id = 615;
	%start_query_id = 9069;
	%end_query_id = 9069;
	%debug_shot = 'shot148_1058';
	
	% Push DOWN
	start_video_id = 804;
	end_video_id = 804;
	start_query_id = 9069;
	end_query_id = 9069;
	debug_shot = 'shot196_1078';
	
	%start_video_id = 917;
	%end_video_id = 917;
	%start_query_id = 9069;
	%end_query_id = 9069;
	%debug_shot = 'shot224_1609';
	
	%start_video_id = 708;
	%end_video_id = 708;
	%start_query_id = 9069;
	%end_query_id = 9069;
	%debug_shot = 'shot171_181';
	
	%start_video_id = 745;
	%end_video_id = 745;
	%start_query_id = 9069;
	%end_query_id = 9069;
	%debug_shot = 'shot181_1784';
	
	debug_mode = true;
end

%% base level path configuration
ex_bounding_box = 0;
LOOK_UP_DIR = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/2.1.run_query2013-new_test2013-new_TiepBoW_No1_10K/tv2013/test2013-new';
BASE_RESULT_DIR = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/2.7.run_query2013-new_test2013-new_TiepBoW_No1_10K_combine_DPM_RANSAC_max_remove_epsilon_remove_accumulation_nd_nfg_4/tv2013/test2013-new';
BASE_TMP_RANSAC_DIR = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/2.7.run_query2013-new_test2013-new_TiepBoW_No1_10K_combine_DPM_RANSAC/tv2013/test2013-new';
LOG_FILE = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/log/rerank_using_DMP.txt';
BASE_DIR = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/';
BASE_CONFIG_DIR = [BASE_DIR 'metadata/keyframe-5/tv2013/']; 
BASE_LOOKUP_PATH = [BASE_CONFIG_DIR 'test2013-new/'];
LOCAL_DIR = '/tmp/dpm/';

% Change when using different features BoW
qr_raw_bow = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/bow/fg+bg_0.1_hesaff_rootsift_noangle_akmeans_1000000_100000000_50_kdtree_8_800_kdtree_3_0.0125/raw_bow.mat';
db_quant_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/hesaff_rootsift_noangle_cluster/akmeans_1000000_100000000_50/kdtree_8_800/v1_f1_1_sub_quant';
db_frame_info_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/hesaff_rootsift_noangle_mat';
addpath('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/code');
addpath('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/code/web');
run('/net/per610a/export/das11f/ledduy/plsang/nvtiep/libs/vlfeat-0.9.18/toolbox/vl_setup.m');

% Create result folder
if ~exist(BASE_RESULT_DIR, 'dir')
	mkdir(BASE_RESULT_DIR);
	fileattrib(BASE_RESULT_DIR, '+w', 'a');
end

% Run RANSAC for fg only, bg only, both bg and fg
if is_compute_RANSAC==1
	for q_id = start_query_id:end_query_id
		qr_shotID = num2str(q_id);
		
		% create folder containing result file for this video id
		final_result_dir = fullfile(BASE_TMP_RANSAC_DIR, qr_shotID);
		if ~exist(final_result_dir, 'dir')
			try
				mkdir(final_result_dir);
				% make folder writable by all users
				fileattrib(final_result_dir, '+w', 'a');
			catch
				error('error creating final_result_dir');
			end
		end

		% get list of query images
		re = 'frames_png/(.*)/(.*.png)';		
		load(qr_raw_bow); 			% Dung de lay thong tin query_filenames va frame_quant_info
		nquery = length(query_filenames);
		query_set = cell(0);
		count = 0;
		query_id = -1;
		for i = 1:nquery	% find all frames of given query
			for topic_id = 1:length(query_filenames{i})
				[rematch, retok] = regexp(query_filenames{i}{topic_id}, re, 'match', 'tokens');
				if strcmp(qr_shotID, retok{1}{1})
					count = count+1;
					query_set{count} = fullfile('/net/per610a/export/das11g/caizhizhu/ins/ins2013/query/frames_png', retok{1}{1}, retok{1}{2});
					query_id = i;
				end
			end
			if query_id ~= -1
				break;
			end
		end

		for id = start_video_id:end_video_id
			lookup_fname = [qr_shotID,'/TRECVID2013_', num2str(id),'.res']
			% Log
			logfile=fopen(LOG_FILE,'a');
			fprintf(logfile, 'Query: %d. VidId: %d - (%d - %d)\n', q_id, id, start_video_id, end_video_id);
			fclose(logfile);
			fileattrib(LOG_FILE, '+w', 'a');
			% get list of keyframes to perform detection
			lookup_path = fullfile(LOOK_UP_DIR, lookup_fname);
			lookup_file = fopen(lookup_path, 'r');
			db_shot_list = [];
			if lookup_file ~= -1
				db_shot_list = textscan(lookup_file, '%s #$# %s #$# %s');
				db_shot_list = db_shot_list{1};
			else
				continue;
			end
			fclose(lookup_file);
			
			% create result file at temporary directory
			video_id = ['TRECVID2013_' int2str(id)];
			result_file = fullfile(final_result_dir, [video_id '.mat']);
			
			% RANSAC
			if exist(result_file, 'file')
				continue;
			else
				% inliers_struct{ shot_list frame_name fg_inlier_loc bg_inlier_loc fg_bg_inlier{fg_loc bg_loc} } Luu trong file "TRECVID2013_XXX.mat"
				inliers_struct.shot_list = db_shot_list;
				nshot = length(inliers_struct.shot_list);
				inliers_struct.frame_name = cell(1, nshot);
				clear db_shot_list;
				
				for shot_idx=1:nshot	% For all shots in the res file
					db_shotID = inliers_struct.shot_list{shot_idx};
					% get list of visual words of db_image
					db_quant_file = fullfile(db_quant_dir, [db_shotID,'.mat']);
					load(db_quant_file);		% dung de lay thong tin bins
					db_frame_info_file = fullfile(db_frame_info_dir, [db_shotID,'.mat']);
					load(db_frame_info_file);	% dung de lay thong tin clip_frame
					nframe_per_shot = length(clip_frame);
					% Tim frame_id cua db_image trong danh sach frame cua db_shotID
					db_set=clip_frame;
					inliers_struct.frame_name{shot_idx} = db_set;
					inliers_struct.fg_inlier_loc{shot_idx} = cell(1,nframe_per_shot);
					inliers_struct.bg_inlier_loc{shot_idx} = cell(1,nframe_per_shot);
					inliers_struct.fg_bg_inlier{shot_idx} = cell(1,nframe_per_shot);
					
					% Use RANSAC to find no. of inliers
					for db_frame_id=1:nframe_per_shot	% For all frames of a shot
						% db image info
						db_words_id = bins{db_frame_id};
						db_keypoint = round(clip_kp{db_frame_id}(:,:));
						inliers_struct.fg_inlier_loc{shot_idx}{db_frame_id} = cell(1, length(query_set));
						
						for topic_id = 1:length(query_set)	% For all frames of a query
							% parse shot_id + frame name
							[rematch, retok] = regexp(query_set{topic_id}, re, 'match', 'tokens');
							qr_fname = retok{1}{2};
							
							% query image info
							fg_idx = frame_quant_info{query_id}.fg_index{topic_id};
							bg_idx = frame_quant_info{query_id}.bg_index{topic_id};
							qr_words_id_fg = frame_quant_info{query_id}.valid_bins{topic_id}(:,fg_idx);
							qr_words_id_bg = frame_quant_info{query_id}.valid_bins{topic_id}(:,bg_idx);
							qr_keypoint_fg = round(frame_quant_info{query_id}.query_kp{topic_id}(:, fg_idx));
							qr_keypoint_bg = round(frame_quant_info{query_id}.query_kp{topic_id}(:, bg_idx));
							
							[shared_words_fg, iqr_fg, idb_fg] = intersect(qr_words_id_fg(:), db_words_id);
							[shared_words_bg, iqr_bg, idb_bg] = intersect(qr_words_id_bg(:), db_words_id);

							knn=size(qr_words_id_fg,1);
							iqr_fg = floor((iqr_fg+knn-1)/knn);
							iqr_bg = floor((iqr_bg+knn-1)/knn);

							% RANSAC on foreground only
							frame1 = qr_keypoint_fg(:, iqr_fg);
							frame2 = db_keypoint(:,idb_fg);
							matches = [1:size(frame1,2); 1:size(frame1,2)];
							if size(frame1, 2) > 0
								[inliers, H] = geometricVerification(frame1, frame2, matches, 'numRefinementIterations', 10);
								%num_bg_inliers = num_bg_inliers+length(inliers);
								inliers_struct.fg_inlier_loc{shot_idx}{db_frame_id}{topic_id} = frame2(1:2,inliers); % Chi lay nhung points tren db image
							end
							
							% RANSAC on background only
							frame1 = qr_keypoint_bg(:, iqr_bg);
							frame2 = db_keypoint(:,idb_bg);
							matches = [1:size(frame1,2); 1:size(frame1,2)];
							if size(frame1, 2) > 0
								[inliers, H] = geometricVerification(frame1, frame2, matches, 'numRefinementIterations', 10);
								%num_bg_inliers = num_bg_inliers+length(inliers);
								inliers_struct.bg_inlier_loc{shot_idx}{db_frame_id}{topic_id} = frame2(1:2,inliers); % Chi lay nhung points tren db image
							end
							% RANSAC on foreground and background
							frame1 = [qr_keypoint_fg(:, iqr_fg) qr_keypoint_bg(:, iqr_bg)];
							frame2 = [db_keypoint(:,idb_fg) db_keypoint(:,idb_bg)];
							matches = [1:size(frame1,2); 1:size(frame1,2)];
							
							if size(frame1, 2) > 0
								[inliers, H] = geometricVerification(frame1, frame2, matches, 'numRefinementIterations', 10);
								inliers_struct.fg_bg_inlier{shot_idx}{db_frame_id}{topic_id}.fg_loc = frame2(1:2,inliers<=length(iqr_fg)); % Chi lay nhung points tren db image
								inliers_struct.fg_bg_inlier{shot_idx}{db_frame_id}{topic_id}.bg_loc = frame2(1:2,inliers<=length(iqr_bg)); % Chi lay nhung points tren db image
							end
						end
						% write results
						%fprintf(fout, '%s #$# %s #$# %f \n', [db_shotID '_KSC' db_fname], db_shotID, num_bg_inliers);
					end
					% free memory
				end
				save(result_file, 'inliers_struct', '-v6');
				fileattrib(result_file, '+w', 'a');
			end
		end
	end
end

if is_using_RANSAC
	% Combine with DPM
	LOOK_UP_DIR = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/2.2.run_query2013-new_test2013-new_TiepBoW_No1_10K_combine_DPM/tv2013/test2013-new';
	DPM_ORG_FUSION_DIR = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/2.11.run_fusion2013-TiepBoW_No1_10K+DPM[R1=R2.1xw1=2.0-R2=R2.2xw2=1.0-Norm=1]/tv2013/test2013-new';
	re = '.*_KSC(.*)';
	for q_id = start_query_id:end_query_id	% Duyet qua tat ca cac cau query
		qr_shotID = num2str(q_id);
		final_result_dir = fullfile(BASE_RESULT_DIR, qr_shotID);
		final_result_local_dir = fullfile(LOCAL_DIR, qr_shotID);
		% create folder if not existing
		if ~exist(final_result_dir, 'dir')
			mkdir(final_result_dir);
			fileattrib(final_result_dir, '+w', 'a');
		end
		if ~exist(final_result_local_dir, 'dir')
			mkdir(final_result_local_dir);
		end
		
		% Load Original fusion using DPM and BoW: w1 x R_bow + w2 x R_dpm
		fusion_res_file = fullfile(DPM_ORG_FUSION_DIR, qr_shotID , [qr_shotID, '.res']);
		fid = fopen(fusion_res_file);
		dpm_fusion = textscan(fid, '%s #$# %s #$# %f');
		fclose(fid);
		
		% Load DPM scale factor
		scale_factor_file = fullfile(BASE_CONFIG_DIR, [qr_shotID '.cfg']);
		scale_reg = 'Scale : (.*)';
		fid = fopen(scale_factor_file);
		[rematch, retok] = regexp(strtrim(fgetl(fid)), scale_reg, 'match', 'tokens');
		scale_factor = 1.0/str2double(retok{end}{1});
		fclose(fid);

		for id = start_video_id:end_video_id
			fprintf('\rQuery %d, Video: %d - (%d - %d)', q_id, id, start_video_id, end_video_id);
			lookup_fname = [qr_shotID,'/TRECVID2013_', num2str(id),'.res'];
			% Write Log file
			logfile=fopen(LOG_FILE,'a');
			fprintf(logfile, '\r Query: %d. VidId: %d - (%d - %d)\n', q_id, id, start_video_id, end_video_id);
			fclose(logfile);
			fileattrib(LOG_FILE, '+w', 'a');
			
			% Check .res file already existed in data server or not?
			dpm_ransac_res_file = fullfile(final_result_dir, ['/TRECVID2013_', num2str(id),'.res']);
			if exist(dpm_ransac_res_file, 'file')
				continue;
			end
			
			% Load DPM .res files to get bounding box information
			dpm_res_file = fullfile(LOOK_UP_DIR,lookup_fname);
			if ~exist(dpm_res_file, 'file')
				continue;
			end
			fid = fopen(dpm_res_file, 'r');
			C = textscan(fid, '%s #$# %s #$# %f #$# %f #$# %f #$# %f #$# %f #$# %f');
			frame_names = C{1};
			shot_id = C{2};
			%score = C{3};
			left = C{4}.*scale_factor-ex_bounding_box;
			top = C{5}.*scale_factor-ex_bounding_box;
			right = C{6}.*scale_factor+ex_bounding_box;
			bottom = C{7}.*scale_factor+ex_bounding_box;
			fclose(fid);
			clear C;
			
			% Load inlier files file using RANSAC from previous step
			ransac_inlier_file = fullfile(BASE_TMP_RANSAC_DIR, qr_shotID, ['/TRECVID2013_', num2str(id),'.mat']);
			load(ransac_inlier_file);
					
			% Find common shot id and Fuse score
			dpm_ransac_res_local_file = fullfile(final_result_local_dir, ['/TRECVID2013_', num2str(id),'.res']);
			if ~debug_mode
				fid = fopen(dpm_ransac_res_local_file, 'w');
			end

			nshot = length(inliers_struct.shot_list);
			for shot_idx = 1:nshot	% duyet qua tat ca cac shot trong danh sach cua RANSAC
				shot = inliers_struct.shot_list{shot_idx};
				if debug_mode && ~strcmp(shot, debug_shot)
					continue;
				end
				frame_locs = find(ismember(shot_id, shot));	% tim nhung frameID trong DPM .res co shot ID giong voi shotID cua RANSAC
				
				N_fg = 0; % co the nam trong lan ko nam trong DPM region??!!??
				N_bg = 0;
				
				[~, previous_score_id] = ismember(shot, dpm_fusion{1});
				P_score = dpm_fusion{3}(previous_score_id);
				new_scores = [];
				for frame_idx=1:length(frame_locs) % duyet qua tat ca cac frame ma co su dung DPM
					[rematch, frame_name] = regexp(frame_names{frame_locs(frame_idx)}, re, 'match', 'tokens');
					[isa, loc] = ismember(frame_name{end}{1}, inliers_struct.frame_name{shot_idx});
					
					fg_kp = [];
					bg_kp = [];
					if ~isempty(inliers_struct.fg_inlier_loc{shot_idx}{loc})
						% merge all shared words
						fg_kp = [inliers_struct.fg_inlier_loc{shot_idx}{loc}{:}];
						fg_kp = unique(fg_kp', 'rows')';
					end
					if ~isempty(inliers_struct.bg_inlier_loc{shot_idx}{loc})
						% merge all shared words
						nquery = length(inliers_struct.bg_inlier_loc{shot_idx}{loc});
						for ii = 1:nquery
							if size(inliers_struct.bg_inlier_loc{shot_idx}{loc}{ii},2) < 4
								inliers_struct.bg_inlier_loc{shot_idx}{loc}{ii} = [];
							end
						end
						bg_kp = [inliers_struct.bg_inlier_loc{shot_idx}{loc}{:}];
						bg_kp = unique(bg_kp', 'rows')';
					end
					
					%N_fg = N_fg+size(fg_kp,2);
					%N_bg = N_bg+size(bg_kp,2); % version 2.0
					N_bg = size(bg_kp,2);		% version 3.0
					N_fg = size(fg_kp,2);		% version 4.0, bo DPM
					Nd = 0;
					if ~isempty(fg_kp)
						Nd = sum(fg_kp(1,:)>=left(frame_locs(frame_idx))&fg_kp(1,:)<=right(frame_locs(frame_idx))&...
								 fg_kp(2,:)>=top(frame_locs(frame_idx))&fg_kp(2,:)<=bottom(frame_locs(frame_idx)));
					end
					%% compute Nd, N_fg, N_bg, p_score
					% Nd = number of shares words between DPM region and background inliers
					% N_fg = number of shares words in foreground region of query image and frames images
					% N_bg = number of shares words in foreground region of query image and frames images				
					%new_scores(end+1) = exp(Nd)*log2(max(2,N_bg))*P_score;				% Cong thuc goc
					%new_scores(end+1) = (Nd+0.001)*log2(max(2,N_bg))*P_score; 	% binh thuong, co epsilon
					%new_scores(end+1) = max(1,Nd)*log2(max(2,N_bg))*P_score; 	% remove epsilon
					%new_scores(end+1) = max(1,N_fg)*log2(max(2,N_bg))*P_score; 	% remove DPM
					%new_scores(end+1) = max(1,Nd)*max(1,Nd)*max(1,N_fg-Nd)*log2(max(2,N_bg))*P_score; 	% using both Nd and Nfg
					%new_scores(end+1) = max(1,Nd)*max(1,N_fg-Nd)*log2(max(2,N_bg))*P_score; 			% using both Nd and Nfg 1
					%new_scores(end+1) = max(1,Nd)*max(1,Nd)*max(1,N_fg)*log2(max(2,N_bg))*P_score; 	% using both Nd and Nfg 2
					%new_scores(end+1) = max(1,Nd)*max(1,Nd)*max(1,Nd)*max(1,N_fg-Nd)*log2(max(2,N_bg))*P_score; 	% using both Nd and Nfg 3 
					new_scores(end+1) = max(1,Nd)*max(1,N_fg-Nd)*log2(max(2,N_bg))*P_score; 	% using both Nd and Nfg 4
					%new_scores(end+1) = P_score; 								% fuse 2 : 1
					
					if debug_mode
						fullfile('/net/per610a/export/das11g/caizhizhu/ins/ins2013/frames_png/', shot, [frame_name{end}{1} '.png'])
						I = imread(fullfile('/net/per610a/export/das11g/caizhizhu/ins/ins2013/frames_png/', shot, [frame_name{end}{1} '.png']));
						figure; imshow(I); hold on;
						xx = [left(frame_locs(frame_idx))  right(frame_locs(frame_idx)) right(frame_locs(frame_idx)) left(frame_locs(frame_idx)) left(frame_locs(frame_idx))];
						yy = [top(frame_locs(frame_idx)) top(frame_locs(frame_idx)) bottom(frame_locs(frame_idx)) bottom(frame_locs(frame_idx)) top(frame_locs(frame_idx))];
						plot(xx, yy, 'r');
						if ~isempty(fg_kp)
							plot(fg_kp(1,:), fg_kp(2,:), 'g+');
						end
						if ~isempty(bg_kp)
							plot(bg_kp(1,:), bg_kp(2,:), 'b+');
						end
						pause;
					end
				end
				% Write to output file
				if ~debug_mode
					fprintf(fid, '%s #$# %s #$# %f\n', shot, shot, max(new_scores));
				end
			end
			if ~debug_mode
				fclose(fid);
				status = unix(['mv ' dpm_ransac_res_local_file ' ' dpm_ransac_res_file]);
				if status == 0
					try
						delete(dpm_ransac_res_local_file);
					catch
						error('Cannot delete temporary result file');
					end
				end
				fileattrib(dpm_ransac_res_file, '+w', 'a');
			end
			clear inliers_struct fg_locs fg_kp bg_kp frame_names shot_id score left top right bottom frame_locs
		end
	end
else
	% Combine with DPM
	LOOK_UP_DIR = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/2.2.run_query2013-new_test2013-new_TiepBoW_No1_10K_combine_DPM/tv2013/test2013-new';
	DPM_ORG_FUSION_DIR = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/2.11.run_fusion2013-TiepBoW_No1_10K+DPM[R1=R2.1xw1=2.0-R2=R2.2xw2=1.0-Norm=1]/tv2013/test2013-new';
	load(qr_raw_bow); 			% Dung de lay thong tin query_filenames va frame_quant_info
	nquery = length(query_filenames);
	for q_id = start_query_id:end_query_id	% Duyet qua tat ca cac cau query
		qr_shotID = num2str(q_id);
		final_result_dir = fullfile(BASE_RESULT_DIR, qr_shotID);
		final_result_local_dir = fullfile(LOCAL_DIR, qr_shotID);
		% create folder if not existing
		if ~exist(final_result_dir, 'dir')
			mkdir(final_result_dir);
			fileattrib(final_result_dir, '+w', 'a');
		end
		if ~exist(final_result_local_dir, 'dir')
			mkdir(final_result_local_dir);
		end
		
		% Load Original fusion using DPM and BoW: w1 x R_bow + w2 x R_dpm
		fusion_res_file = fullfile(DPM_ORG_FUSION_DIR, qr_shotID , [qr_shotID, '.res']);
		fid = fopen(fusion_res_file);
		dpm_fusion = textscan(fid, '%s #$# %s #$# %f');
		fclose(fid);
		
		% Load DPM scale factor
		scale_factor_file = fullfile(BASE_CONFIG_DIR, [qr_shotID '.cfg']);
		scale_reg = 'Scale : (.*)';
		fid = fopen(scale_factor_file);
		[rematch, retok] = regexp(strtrim(fgetl(fid)), scale_reg, 'match', 'tokens');
		scale_factor = 1.0/str2double(retok{end}{1});
		fclose(fid);
		
		% get list of query images
		re = 'frames_png/(.*)/(.*.png)';		
		query_set = cell(0);
		count = 0;
		query_id = -1;
		for i = 1:nquery	% find all frames of given query
			for topic_id = 1:length(query_filenames{i})
				[rematch, retok] = regexp(query_filenames{i}{topic_id}, re, 'match', 'tokens');
				if strcmp(qr_shotID, retok{1}{1})
					count = count+1;
					query_set{count} = fullfile('/net/per610a/export/das11g/caizhizhu/ins/ins2013/query/frames_png', retok{1}{1}, retok{1}{2});
					query_id = i;
				end
			end
			if query_id ~= -1
				break;
			end
		end

		for id = start_video_id:end_video_id
			fprintf('\rQuery %d, Video: %d - (%d - %d)', q_id, id, start_video_id, end_video_id);
			lookup_fname = [qr_shotID,'/TRECVID2013_', num2str(id),'.res'];
			% Write Log file
			logfile=fopen(LOG_FILE,'a');
			fprintf(logfile, '\r Query: %d. VidId: %d - (%d - %d)\n', q_id, id, start_video_id, end_video_id);
			fclose(logfile);
			fileattrib(LOG_FILE, '+w', 'a');
			
			% Check .res file already existed in data server or not?
			dpm_nonRANSAC_res_file = fullfile(final_result_dir, ['/TRECVID2013_', num2str(id),'.res']);
			if exist(dpm_nonRANSAC_res_file, 'file')
				continue;
			end
			
			% Load DPM .res files to get bounding box information
			dpm_res_file = fullfile(LOOK_UP_DIR,lookup_fname);
			if ~exist(dpm_res_file, 'file')
				continue;
			end
			fid = fopen(dpm_res_file, 'r');
			C = textscan(fid, '%s #$# %s #$# %f #$# %f #$# %f #$# %f #$# %f #$# %f');
			frame_names = C{1};
			shot_id = C{2};
			%score = C{3};
			left = C{4}.*scale_factor-ex_bounding_box;
			top = C{5}.*scale_factor-ex_bounding_box;
			right = C{6}.*scale_factor+ex_bounding_box;
			bottom = C{7}.*scale_factor+ex_bounding_box;
			fclose(fid);
			clear C;
					
			% Find common shot id and Fuse score
			dpm_nonransac_res_local_file = fullfile(final_result_local_dir, ['/TRECVID2013_', num2str(id),'.res']);
			if ~debug_mode
				fid = fopen(dpm_nonransac_res_local_file, 'w');
			end

			shot_list = unique(shot_id);
			nshot = length(shot_list);
			for shot_idx = 1:nshot	% duyet qua tat ca cac shot trong danh sach cua RANSAC
				shot = shot_list{shot_idx};
				if debug_mode && ~strcmp(shot, debug_shot)
					continue;
				end
				frame_locs = find(ismember(shot_id, shot));	% tim nhung frameID trong DPM .res co shot ID giong voi shotID cua RANSAC
				
				N_fg = 0; % co the nam trong lan ko nam trong DPM region??!!??
				N_bg = 0;
				
				% get list of visual words of db_image
				db_quant_file = fullfile(db_quant_dir, [shot,'.mat']);
				load(db_quant_file);		% dung de lay thong tin bins
				db_frame_info_file = fullfile(db_frame_info_dir, [shot,'.mat']);
				load(db_frame_info_file);	% dung de lay thong tin clip_frame
				nframe_per_shot = length(clip_frame);
				% Tim frame_id cua db_image trong danh sach frame cua db_shotID
				db_set=clip_frame;
				
				[~, previous_score_id] = ismember(shot, dpm_fusion{1});
				P_score = dpm_fusion{3}(previous_score_id);
				new_scores = [];
				for frame_idx=1:length(frame_locs) % duyet qua tat ca cac frame ma co su dung DPM
					re = '.*_KSC(.*)';
					[rematch, frame_name] = regexp(frame_names{frame_locs(frame_idx)}, re, 'match', 'tokens');
					[isa, loc] = ismember(frame_name{end}{1}, db_set);
					
					db_words_id = bins{loc};
					db_keypoint = round(clip_kp{loc}(:,:));
					
					fg_kp = [];
					bg_kp = [];
					for topic_id = 1:length(query_set)	% For all frames of a query					
						% query image info
						fg_idx = frame_quant_info{query_id}.fg_index{topic_id};
						bg_idx = frame_quant_info{query_id}.bg_index{topic_id};
						qr_words_id_fg = frame_quant_info{query_id}.valid_bins{topic_id}(:,fg_idx);
						qr_words_id_bg = frame_quant_info{query_id}.valid_bins{topic_id}(:,bg_idx);
						%qr_keypoint_fg = round(frame_quant_info{query_id}.query_kp{topic_id}(:, fg_idx));
						%qr_keypoint_bg = round(frame_quant_info{query_id}.query_kp{topic_id}(:, bg_idx));
						
						[shared_words_fg, iqr_fg, idb_fg] = intersect(qr_words_id_fg(:), db_words_id);
						[shared_words_bg, iqr_bg, idb_bg] = intersect(qr_words_id_bg(:), db_words_id);

						%knn=size(qr_words_id_fg,1);
						%iqr_fg = floor((iqr_fg+knn-1)/knn);
						%iqr_bg = floor((iqr_bg+knn-1)/knn);

						% RANSAC on foreground only
						if ~isempty(db_keypoint)
							frame2 = db_keypoint(1:2,idb_fg);
							fg_kp = [fg_kp frame2];
							
							% RANSAC on background only
							frame2 = db_keypoint(1:2,idb_bg);
							bg_kp = [bg_kp frame2];
						end
					end
					fg_kp = unique(fg_kp', 'rows')';
					bg_kp = unique(bg_kp', 'rows')';
					
					%N_fg = N_fg+size(fg_kp,2);
					%N_bg = N_bg+size(bg_kp,2);%% checkkkkkkkkkkkkkkk lai o tren dung ransac
					N_bg = size(bg_kp,2);
					Nd = 0;
					if ~isempty(fg_kp)
						Nd = sum(fg_kp(1,:)>=left(frame_locs(frame_idx))&fg_kp(1,:)<=right(frame_locs(frame_idx))&...
								 fg_kp(2,:)>=top(frame_locs(frame_idx))&fg_kp(2,:)<=bottom(frame_locs(frame_idx)));
					end
					%% compute Nd, N_fg, N_bg, p_score
					% Nd = number of shares words between DPM region and background inliers
					% N_fg = number of shares words in foreground region of query image and frames images
					% N_bg = number of shares words in foreground region of query image and frames images				
					%new_scores(end+1) = exp(Nd)*log2(max(2,N_bg))*P_score;				% Cong thuc goc
					%new_scores(end+1) = (Nd+0.001)*log2(max(2,N_bg))*P_score; 	% binh thuong, co epsilon
					new_scores(end+1) = max(1,Nd)*log2(max(2,N_bg))*P_score; 	% remove epsilon
					%new_scores(end+1) = P_score; 								% fuse 2 : 1
					
					if debug_mode
						fullfile('/net/per610a/export/das11g/caizhizhu/ins/ins2013/frames_png/', shot, [frame_name{end}{1} '.png'])
						I = imread(fullfile('/net/per610a/export/das11g/caizhizhu/ins/ins2013/frames_png/', shot, [frame_name{end}{1} '.png']));
						figure; imshow(I); hold on;
						xx = [left(frame_locs(frame_idx))  right(frame_locs(frame_idx)) right(frame_locs(frame_idx)) left(frame_locs(frame_idx)) left(frame_locs(frame_idx))];
						yy = [top(frame_locs(frame_idx)) top(frame_locs(frame_idx)) bottom(frame_locs(frame_idx)) bottom(frame_locs(frame_idx)) top(frame_locs(frame_idx))];
						plot(xx, yy, 'r');
						if ~isempty(fg_kp)
							plot(fg_kp(1,:), fg_kp(2,:), 'g+');
						end
						if ~isempty(bg_kp)
							plot(bg_kp(1,:), bg_kp(2,:), 'b+');
						end
						pause;
					end
				end
				% Write to output file
				if ~debug_mode
					fprintf(fid, '%s #$# %s #$# %f\n', shot, shot, max(new_scores));
				end
			end
			if ~debug_mode
				fclose(fid);
				status = unix(['mv ' dpm_nonransac_res_local_file ' ' dpm_nonRANSAC_res_file]);
				if status == 0
					try
						delete(dpm_nonransac_res_local_file);
					catch
						error('Cannot delete temporary result file');
					end
				end
				fileattrib(dpm_nonRANSAC_res_file, '+w', 'a');
			end
			clear inliers_struct fg_locs fg_kp bg_kp frame_names shot_id score left top right bottom frame_locs
		end
	end
end

quit

end
