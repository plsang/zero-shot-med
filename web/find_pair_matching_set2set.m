function [nshares_fg, nshares_bg] = find_pair_matching_set2set(qr_shotID, db_shotID, output_dir)

if nargin == 0 % default data used for testing
	qr_shotID = '9069';
	db_shotID = 'shot11_487';
	output_dir= '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/web/';
end
% base dir
db_quant_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/hesaff_rootsift_noangle_cluster/akmeans_1000000_100000000_50/kdtree_8_800/v1_f1_1_sub_quant';
db_frame_info_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/hesaff_rootsift_noangle_mat';
qr_raw_bow = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/bow/fg+bg_0.1_hesaff_rootsift_noangle_akmeans_1000000_100000000_50_kdtree_8_800_kdtree_3_0.0125/raw_bow.mat';

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

% get list of visual words of qr_image
re = 'frames_png/(.*)/(.*.png)';		
load(qr_raw_bow); 			% Dung de lay thong tin query_filenames va frame_quant_info
nquery = length(query_filenames);
query_set = cell(0);
count = 0;
for query_id = 1:nquery
	for topic_id = 1:length(query_filenames{query_id})
		[rematch, retok] = regexp(query_filenames{query_id}{topic_id}, re, 'match', 'tokens');
		if strcmp(qr_shotID, retok{1}{1})
			count = count+1;
			query_set{count} = fullfile('/net/per610a/export/das11g/caizhizhu/ins/ins2013/query/frames_png', retok{1}{1}, retok{1}{2});
		end
	end
end

for i = 1:length(query_set)
	for j=1:3:length(db_set)
		% parse shot_id + frame name

		[rematch, retok] = regexp(query_set{i}, re, 'match', 'tokens');
		qr_shotID = retok{1}{1};
		qr_fname = retok{1}{2};

		db_img = ['/net/per610a/export/das11g/caizhizhu/ins/ins2013/frames_png/',db_shotID, '/', db_set{j},'.png'];
		[rematch, retok] = regexp(db_img, re, 'match', 'tokens');
		db_shotID = retok{1}{1};
		db_fname = retok{1}{2};

		output_image = fullfile(output_dir, [qr_fname,'_',db_shotID,db_fname]);
		
		%find_pair_matching(query_set{i}, db_img, output_image)
		find_pair_matching_RANSAC(query_set{i}, db_img, output_image)
	end
end
quit
end