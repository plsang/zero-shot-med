function [nshares_fg, nshares_bg] = find_pair_matching_set2set_BOW(output_file, runID, qr_shotID, db_shotID, output_dir)

if nargin == 0 % default data used for testing
	qr_shotID = '9069';
	db_shotID = 'shot11_487';
	output_dir= '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/web/';
	runID = '3.1.run_query2013-new_test2013-new_TiepBoW_10K';
end
% base dir
if ~isempty(strfind(runID, 'No1'))
	db_quant_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/hesaff_rootsift_noangle_cluster/akmeans_1000000_100000000_50/kdtree_8_800/v1_f1_1_sub_quant';
	db_frame_info_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/hesaff_rootsift_noangle_mat';
	qr_raw_bow = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/bow/fg+bg_0.1_hesaff_rootsift_noangle_akmeans_1000000_100000000_50_kdtree_8_800_kdtree_3_0.0125/raw_bow.mat';
else
	db_quant_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/perdoch_hesaff_rootsift_cluster/akmeans_1000000_100000000_50/kdtree_8_800/v1_f1_3_sub_quant';
	db_frame_info_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/perdoch_hesaff_rootsift_mat';
	qr_raw_bow = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/bow/fg+bg_0.1_perdoch_hesaff_rootsift_akmeans_1000000_100000000_50_kdtree_8_800_kdtree_3_0.0125/raw_bow.mat';
end
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

% get list of query image
re = 'frames_png/(.*)/(.*.png)';		
load(qr_raw_bow); 			% Dung de lay thong tin query_filenames va frame_quant_info
nquery = length(query_filenames);
query_set = cell(0);
count = 0;
for query_id = 1:nquery
	for topic_id = 1:length(query_filenames{query_id})
		[rematch, retok] = regexp(query_filenames{query_id}{topic_id}, re, 'match', 'tokens');
		if strcmp(qr_shotID, retok{1}{1})
			tmp = query_id;
			count = count+1;
			query_set{count} = fullfile('/net/per610a/export/das11g/caizhizhu/ins/ins2013/query/frames_png', retok{1}{1}, retok{1}{2});
		end
	end
end
query_id = tmp;

sampling_rate = max(floor(length(db_set)/5),1);
list_output = cell(0);
list_score = zeros(0);
list_fg = zeros(0);
list_bg = zeros(0);
num_output = 0;

for i = 1:length(query_set)
	for j=1:sampling_rate:length(db_set)
		% parse shot_id + frame name

		[rematch, retok] = regexp(query_set{i}, re, 'match', 'tokens');
		qr_shotID = retok{1}{1};
		qr_fname = retok{1}{2};

		db_img = ['/net/per610a/export/das11g/caizhizhu/ins/ins2013/frames_png/',db_shotID, '/', db_set{j},'.png'];
		[rematch, retok] = regexp(db_img, re, 'match', 'tokens');
		db_shotID = retok{1}{1};
		db_fname = retok{1}{2};

		output_image = fullfile(output_dir, [qr_fname,'_',db_shotID,db_fname]);
		
		[score, new_output_img, nfg, nbg] = find_pair_matching(query_set{i}, query_id, i, db_img, j, output_image, runID, frame_quant_info, query_filenames, topic_bows, bins, clip_frame, clip_kp);
		%find_pair_matching(query_set{i}, db_img, output_image, runID)
		
		num_output = num_output+1;
		[pathstr,name,ext] = fileparts(new_output_img);
		list_output{num_output} = name;
		list_score(num_output) = score;
		list_fg(num_output) = nfg;
		list_bg(num_output) = nbg;
	end
end

% Write output file including <output image> #$# <score>
[sorted_list, idx] = sort(list_score, 'descend');
fid = fopen(output_file, 'w');
for i=1:length(sorted_list)
	fprintf(fid, '%s #$# %d #$# %d #$# %f \n', list_output{idx(i)}, list_fg(idx(i)), list_bg(idx(i)), sorted_list(i));
end
fclose(fid);
fileattrib(output_file, '+w', 'a');

quit
end