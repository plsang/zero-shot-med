function comp_dist_query_shot(sID, eID)
query_ids = 9069:9098;

% Load query BoW
load('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/bow/fg+bg_0.1_hesaff_rootsift_noangle_akmeans_1000000_100000000_50_kdtree_8_800_kdtree_3_0.0125/bow_full_notrim_clip_idf_nonorm_-1.mat');
ntopic = size(topic_bows{1}{1},2);
BASE_RESULT_DIR = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/2.1.run_query2013-new_test2013-new_TiepBoW_No1_10K/tv2013/test2013-new/';
BOW_DIR = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/hesaff_rootsift_noangle_cluster/akmeans_1000000_100000000_50/kdtree_8_800/v1_f1_1/sub_bow/';
L1_RAW_RESULT_DIR = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/2.15.run_query2013-new_test2013-new_TiepBoW_No1_10K_recompute_distance_L1/tv2013/test2013-new/';
L2_RAW_RESULT_DIR = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/2.15.run_query2013-new_test2013-new_TiepBoW_No1_10K_recompute_distance_L2/tv2013/test2013-new/';
LOCAL_DIR = '/tmp/dpm/';

for query_index = 1:30
	local_res_dir = fullfile(LOCAL_DIR, num2str(query_ids(query_index)));
	if ~exist(local_res_dir, 'dir')
		mkdir(local_res_dir);
	end
	L1_new_score_dir = fullfile(L1_RAW_RESULT_DIR, num2str(query_ids(query_index)));
	L2_new_score_dir = fullfile(L2_RAW_RESULT_DIR, num2str(query_ids(query_index)));
	if ~exist(L1_new_score_dir, 'dir')
		mkdir(L1_new_score_dir);
		fileattrib(L1_new_score_dir, '+w', 'a');
	end
	if ~exist(L2_new_score_dir, 'dir')
		mkdir(L2_new_score_dir);
		fileattrib(L2_new_score_dir, '+w', 'a');
	end
	
	for vid_idx = sID:eID
		fprintf('\r%d (%d - %d)', vid_idx, sID, eID);
		% Load all shots of a TRECVID_videoID
		old_result_file = fullfile(BASE_RESULT_DIR, num2str(query_ids(query_index)), ['TRECVID2013_' num2str(vid_idx) '.res']);
		fid = fopen(old_result_file, 'r');
		if fid == -1 
			continue
		end
		lst_shots = textscan(fid, '%*s #$# %s #$# %*f');
		lst_shots = lst_shots{1};
		fclose(fid);
		nshot = length(lst_shots);
		
		% Open files to write new score using L1 and L2 distance
		l1_dist_res_local_file = fullfile(local_res_dir, ['TRECVID2013_' num2str(vid_idx) '1.raw']);
		l2_dist_res_local_file = fullfile(local_res_dir, ['TRECVID2013_' num2str(vid_idx) '2.raw']);
		l1_dist_res_file = fullfile(L1_RAW_RESULT_DIR, num2str(query_ids(query_index)), ['TRECVID2013_' num2str(vid_idx) '.raw']);
		l2_dist_res_file = fullfile(L2_RAW_RESULT_DIR, num2str(query_ids(query_index)), ['TRECVID2013_' num2str(vid_idx) '.raw']);
		%if exist(l1_dist_res_file, 'file') && exist(l2_dist_res_file, 'file')
		%		continue;
		%end
		fid_l1 = fopen(l1_dist_res_local_file, 'w');
		fid_l2 = fopen(l2_dist_res_local_file, 'w');
		
		% For each shot in TRECVID2013_VID, recompute score for each query and key frame
		for shot_idx = 1:nshot
			% Load bag of word of current shot
			load([BOW_DIR lst_shots{shot_idx}]);
			frame_bow = full(frame_bow);
			[feat_size, nframe] = size(frame_bow);
			db_bow_l1 = frame_bow./repmat(sum(frame_bow), feat_size, 1);
			db_norm_l2 = sqrt(sum(frame_bow.^2));
			db_bow_l2 = frame_bow./repmat(db_norm_l2, feat_size, 1);
			for topic_idx = 1:ntopic
				qr_bow = full(topic_bows{query_index}{1}(:,topic_idx));
				qr_bow_l1 = repmat(qr_bow/sum(qr_bow), 1, nframe);
				qr_bow_l2 = repmat(qr_bow/norm(qr_bow), 1, nframe);
				fprintf(fid_l1, '%d.%d.src #$# %s #$# ', query_ids(query_index), topic_idx, lst_shots{shot_idx});
				fprintf(fid_l1, '%f ', abs(2 - sum(abs(qr_bow_l1-db_bow_l1))));
				fprintf(fid_l1, '\n');
				fprintf(fid_l2, '%d.%d.src #$# %s #$# ', query_ids(query_index), topic_idx, lst_shots{shot_idx});
				fprintf(fid_l2, '%f ', abs(2 - sum((qr_bow_l2-db_bow_l2).^2)));
				fprintf(fid_l2, '\n');
			end
		end
		fclose(fid_l1);
		fclose(fid_l2);
		unix(['mv ' l1_dist_res_local_file ' ' l1_dist_res_file]);
		unix(['mv ' l2_dist_res_local_file ' ' l2_dist_res_file]);
		fileattrib(l1_dist_res_file, '+w', 'a');
		fileattrib(l2_dist_res_file, '+w', 'a');
	end
end
quit
end