
%if nargin==0
	%res_dir = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/run_query2013-new_test2013-new_TiepBoW_10K_combine_DPM/tv2013/test2013-new';
	%output_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/result/combine_dpm/perf/top10k/';
	%res_dir = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/run_query2013-new_test2013-new_TiepBoW_10K/tv2013/test2013-new';
	%output_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/result/combine_dpm/perf/top10k_org/';
	%res_dir = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/run_query2013-new_test2013-new_dpm/tv2013/test2013-new';
	%output_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/result/combine_dpm/perf/ThayDuy_dpm/';
	res_dir = '/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/run_query2013-new_test2013-new_nsc.bow.dense6mul.sift.Soft-1000.test2013-new.norm1x1_l1norm1/tv2013/test2013-new';
	output_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/result/combine_dpm/perf/dense_dpm/';
%end

rlist_dir = dir(res_dir);

for k=1:length(rlist_dir)
	if rlist_dir(k).name(1) == '.' || ~rlist_dir(k).isdir
		continue;
	end
	
	query_path = fullfile(res_dir, rlist_dir(k).name)
	rfiles_per_query = dir(fullfile(query_path,'*.res'));
	output_file = [output_dir, rlist_dir(k).name ,'.1.list'];
	
	map = containers.Map;
	for i=1:length(rfiles_per_query)
		result_file = fullfile(query_path,rfiles_per_query(i).name);
		file_info = dir(result_file);
		if file_info(1).bytes ==0
			continue;
		end
		fid=fopen(result_file,'r');
		fprintf('\r%d-%d', i, length(rfiles_per_query));
		re = '(shot\d*_\d*)_';
		
		while ~feof(fid)
			line = strtrim(fgetl(fid));
			%C = textscan(line, '%s #$# %s #$# %s');
			C = textscan(line, '%s #$# %s #$# %s #$# %s #$# %s #$# %s #$# %s #$# %s');	
			[rematch, retok] = regexp(C{1}{1}, re, 'match', 'tokens');
			if isKey(map, retok{1}{1})
				map(retok{1}{1}) = max(map(retok{1}{1}), str2double(C{3}{1}));
			else
				map(retok{1}{1}) = str2double(C{3}{1});
			end
		end
		fclose(fid);
	end
	
	fid = fopen(output_file, 'w');
	lstShot = keys(map);
	lstScore = values(map,lstShot);
	lstScore = [lstScore{:}];
	[~, idx] = sort(lstScore,'descend');
	
	for i=1:1000
		fprintf(fid, '%s 0 %s %d %f dense6mul\n', rlist_dir(k).name, lstShot{idx(i)}, i, map(lstShot{idx(i)}));
	end
	fclose(fid);
end