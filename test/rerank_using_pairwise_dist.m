function rerank_using_pairwise_dist(query_index)

query_ids = 9069:9098;
RAW_RESULT_DIR = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/2.15.run_query2013-new_test2013-new_TiepBoW_No1_10K_recompute_distance_L1/tv2013/test2013-new/';
LOCAL_DIR = '/tmp/dpm/';

% Ranking shot using pairwise scores
lst_vid_files = dir(fullfile(RAW_RESULT_DIR, num2str(query_ids(query_index)), '*.raw'));
nvid = length(lst_vid_files);
re = '\#\$\# (.*) \#\$\# (.*)';
for i=1:nvid
	fprintf('\r%d', i);
	% Load raw file
	vid_file = fullfile(RAW_RESULT_DIR, num2str(query_ids(query_index)), lst_vid_files(i).name);
	fid = fopen(vid_file, 'r');
	nshot = 0;
	lst_result_shots = cell(0);
	lst_result_score = cell(0);
	while ~feof(fid)
		line = strtrim(fgetl(fid));
		[rematch, retok] = regexp(line, re, 'match', 'tokens');
		[isa loc] = ismember(retok{1}{1}, lst_result_shots);
		if isa
			lst_result_score{loc} = [lst_result_score{loc} str2num(retok{1}{2})];
		else
			nshot = nshot+1;
			lst_result_shots{nshot} = retok{1}{1};
			lst_result_score{nshot} = str2num(retok{1}{2});
		end
	end
	fclose(fid);
	% Write result
	new_res_file = fullfile(RAW_RESULT_DIR, num2str(query_ids(query_index)), [lst_vid_files(i).name(1:end-3) 'res']);
	fid = fopen(new_res_file, 'w');
	for j=1:nshot
		fprintf(fid, '%s #$# %s #$# %f\n', lst_result_shots{j}, lst_result_shots{j}, max(lst_result_score{j}));
	end
	fclose(fid);
end

end