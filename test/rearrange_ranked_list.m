
%type = 0; % best perf
%type = 1; % bad positive first
type = 2; % good negative first
old_runID = '2.2.run_query2013-new_test2013-new_TiepBoW_No1_10K_combine_DPM';
new_runID = '2.2.run_query2013-new_test2013-new_TiepBoW_No1_10K_combine_DPM_good_negative_to_bad_negative';
ranked_list_dir = ['/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/' old_runID '/tv2013/test2013-new'];
output_dir = ['/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/' new_runID '/tv2013/test2013-new'];

% Read ground truth
fid = fopen('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/ins.search.qrels.tv13', 'r');
lines = textscan(fid, '%d %d %s %d');
fclose(fid);

% Keep only ground truth
gt{1} = lines{1}(lines{4}>0);
gt{2} = lines{3}(lines{4}>0);


query_list_dir = dir(ranked_list_dir);

for k=1:length(query_list_dir)
	if ~query_list_dir(k).isdir || query_list_dir(k).name(1) == '.' || ~query_list_dir(k).isdir
		continue;
	end
	
	% get all result files of a query
	query_path = fullfile(ranked_list_dir, query_list_dir(k).name)
	res_files_of_query = dir(fullfile(query_path,'*.res'));
	
	% prepare output dir
	query_output_path = fullfile(output_dir, query_list_dir(k).name );
	if ~exist(query_output_path, 'dir')
		mkdir(query_output_path);
	end
	
	% Load old ranked list (full)
	fid = fopen(fullfile(ranked_list_dir, [query_list_dir(k).name '.rank']), 'r');
	re = '(.*)\#\$\#(.*)';
	stt = 0;
	full_list_shot = cell(1, 10000);
	full_list_score	= zeros(1, 10000);
	while ~feof(fid)
		line = strtrim(fgetl(fid));
		stt = stt+1;
		[rematch, retok] = regexp(line, re, 'match', 'tokens');
		full_list_shot{stt} = retok{1}{1};
		full_list_score(stt) = str2double(retok{1}{2});
	end
	fclose(fid);
	full_list_shot = full_list_shot(1:stt);
	
	for i=1:length(res_files_of_query)
		input_res_file = fullfile(query_path,res_files_of_query(i).name);
		file_info = dir(input_res_file);
		if file_info(1).bytes ==0
			continue;
		end
		% read all lines
		fid=fopen(input_res_file,'r');
		fprintf('\r%d-%d', i, length(res_files_of_query));
		re = '(shot\d*_\d*)';
		% open output file to rerank
		output_res_file = fullfile(output_dir, query_list_dir(k).name, res_files_of_query(i).name);
		ofid = fopen(output_res_file, 'w');
		
		pos_gt = {gt{2}{gt{1}==str2double(query_list_dir(k).name)}};

		while ~feof(fid)
			line = strtrim(fgetl(fid));
			%C = textscan(line, '%s #$# %s #$# %s');
			C = textscan(line, '%s #$# %s #$# %s #$# %s #$# %s #$# %s #$# %s #$# %s');	
			[rematch, retok] = regexp(C{1}{1}, re, 'match', 'tokens');
			% if shot id is in the ground truth file -> write to output file
			%is_mem = cellfun(@(x) strcmp(x, retok{1}{1}), pos_gt, 'UniformOutput', false);
			is_pos = ismember(retok{1}{1}, pos_gt);
			[is_mem loc] = ismember(C{2}{1}, full_list_shot);
			switch type
			case 0
				if is_pos %sum([is_mem{:}]) > 0
					fprintf(ofid, '%s #$# %s #$# %f #$# %f #$# %d\n', C{1}{1}, C{2}{1}, 1.5+rand(), full_list_score(loc), loc);
				else
					fprintf(ofid, '%s #$# %s #$# %f #$# %f #$# %d\n', C{1}{1}, C{2}{1}, 1-rand(), full_list_score(loc), loc);
				end
			case 1
				score = str2double(C{3}{1});
				if is_pos
					fprintf(ofid, '%s #$# %s #$# %f #$# %f #$# %d\n', C{1}{1}, C{2}{1}, 4-score, full_list_score(loc), loc);
				else
					fprintf(ofid, '%s #$# %s #$# %f #$# %f #$# %d\n', C{1}{1}, C{2}{1}, 1-score, full_list_score(loc), loc);
				end
			case 2
				score = str2double(C{3}{1});
				if is_pos
					fprintf(ofid, '%s #$# %s #$# %f #$# %f #$# %d\n', C{1}{1}, C{2}{1}, 1-score, full_list_score(loc), loc);
				else
					fprintf(ofid, '%s #$# %s #$# %f #$# %f #$# %d\n', C{1}{1}, C{2}{1}, 4+score, full_list_score(loc), loc);
				end
			end
		end

		fclose(fid);
		fclose(ofid);
	end
end
