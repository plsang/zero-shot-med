ranked_list_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/result/best_config/txt/fg+bg_0.1_hesaff_rootsift_noangle_akmeans_1000000_100000000_50_kdtree_8_800_v1_f1_1_avg_pooling_full_notrim_clip_idf_nonorm_kdtree_3_0.0125_-1_dist_avg_autoasym_ivf_0.5';

% Read ground truth
fid = fopen('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/ins.search.qrels.tv13', 'r');
lines = textscan(fid, '%d %d %s %d');
fclose(fid);

% Keep only ground truth
gt{1} = lines{1}(lines{4}>0);
gt{2} = lines{3}(lines{4}>0);


query_list_dir = dir(ranked_list_dir);

topK = 1000:1000:10000;
num_shares = zeros(30, length(topK));
nquery = 0;
for k=1:length(query_list_dir)
	if query_list_dir(k).name(1) == '.'
		continue;
	end
	nquery = nquery+1;
	fprintf('\r%s (9069 - 9098)\n', query_list_dir(k).name(1:4));
	% get all shots from ranked list
	rank_file_of_query = fullfile(ranked_list_dir, query_list_dir(k).name);
	fid = fopen(rank_file_of_query, 'r');
	str = strtrim(fgetl(fid));
	
	re = '(shot\d*_\d*)\(dist_(.*)\)';
	C = cell(0);
	while ~feof(fid)
		str = strtrim(fgetl(fid));
		[rematch, retok] = regexp(str, re, 'match', 'tokens');
		C{end+1} = retok{1}{1};
	end
	fclose(fid);
	
	pos_gt = {gt{2}{gt{1}==str2double(query_list_dir(k).name(1:4))}}; % list shot belong to a query
	nline = length(C);
	for i=1:length(topK)
		list_topK = C(1:topK(i));
		num_shares(nquery,i) = num_shares(nquery,i) + sum(ismember(list_topK, pos_gt));
	end
end

figure; plot(topK, sum(num_shares)); hold on; plot(topK, sum(num_shares), 'r*');
for i=1:nquery
	figure; plot(topK, num_shares(i,:)); hold on; plot(topK, num_shares(i,:), 'r*');
end
