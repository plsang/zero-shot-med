function rerank_using_RANSAC(start_video_id, end_video_id)
if nargin == 0
	%query_id = '9074';
	start_video_id = 437;
	end_video_id = 437;
	%lookup_fname = [query_id,'/TRECVID2013_', num2str(start_video_id),'.res'];
end

% base level path configuration
LOOK_UP_DIR = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/2.1.run_query2013-new_test2013-new_TiepBoW_No1_10K/tv2013/test2013-new';
BASE_RESULT_DIR = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/2.3.run_query2013-new_test2013-new_TiepBoW_No1_10K_combine_RANSAC/tv2013/test2013-new';
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

for q_id = 9069:9098
	qr_shotID = num2str(q_id);
	
	% create folder containing result file for this video id
	final_result_dir = fullfile(BASE_RESULT_DIR, qr_shotID);
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
	for i = 1:nquery
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
		lookup_fname = [qr_shotID,'/TRECVID2013_', num2str(id),'.res'];
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
		result_file = [LOCAL_DIR qr_shotID '_' video_id '.res']
		fout = fopen(result_file,'w');
		
		for shot_idx=1:length(db_shot_list)
			db_shotID = db_shot_list{shot_idx};
			% get list of visual words of db_image
			db_quant_file = fullfile(db_quant_dir, [db_shotID,'.mat']);
			load(db_quant_file);		% dung de lay thong tin bins
			db_frame_info_file = fullfile(db_frame_info_dir, [db_shotID,'.mat']);
			load(db_frame_info_file);	% dung de lay thong tin clip_frame
			nframe_per_shot = length(clip_frame);
			% Tim frame_id cua db_image trong danh sach frame cua db_shotID
			db_set = cell(1, nframe_per_shot);
			for db_frame_id = 1:nframe_per_shot
				db_set{db_frame_id} = clip_frame{db_frame_id};
			end
			% Use RANSAC to find no. of inliers
			for db_frame_id=1:length(db_set)
				% db image info
				db_words_id = bins{db_frame_id};
				db_keypoint = round(clip_kp{db_frame_id}(:,:));
				db_fname = db_set{db_frame_id};
				num_bg_inliers = 0;
				for topic_id = 1:length(query_set)
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

					% RANSAC
					%frame1 = [qr_keypoint_fg(:, iqr_fg) qr_keypoint_bg(:, iqr_bg)];
					%frame2 = [db_keypoint(:,idb_fg) db_keypoint(:,idb_bg)];
					%matches = [1:size(frame1,2); 1:size(frame1,2)];
					frame1 = qr_keypoint_fg(:, iqr_fg);
					frame2 = db_keypoint(:,idb_fg);
					matches = [1:size(frame1,2); 1:size(frame1,2)];
					if size(frame1, 2) > 0
						[inliers, H] = geometricVerification(frame1, frame2, matches, 'numRefinementIterations', 10);
						num_bg_inliers = num_bg_inliers+length(inliers);
					end
				end
				% write results
				fprintf(fout, '%s #$# %s #$# %f \n', [db_shotID '_KSC' db_fname], db_shotID, num_bg_inliers);
			end

			% free memory 

			% write to log
			%disp(['Finish Detection on ' line ' by model query_' qr_shotID ' with scale factor = ' num2str(scale)]);
		end
		
		fclose(fout);
		% move result file to final directory
		status = unix(sprintf('mv %s %s', result_file, fullfile(final_result_dir, [video_id '.res'])));
		if status == 0
			try
				delete(result_file);
			catch
				error('Cannot delete temporary result file');
			end
		end
		% unix(sprintf('rm -r %s', img_dir));
	end		% end detection on this video_id

end
quit

end
