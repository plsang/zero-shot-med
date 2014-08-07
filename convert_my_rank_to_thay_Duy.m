%function convert_my_rank_to_thay_Duy(res_dir)


%if nargin==0
	res_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/result/best_config/txt';
	%res_dir = '/net/per610a/export/das11g/caizhizhu/ins/ins2013/query/result/trecvid_workshop/txt';
	%res_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/result/best_config/perf';
	topK = 10000;
%end

rlist_dir = dir(res_dir);

load('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS2013/metadata/keyframe-5/lstFileName_new.mat')
reg = '(TRECVID2013_\d*)/(shot\d*_\d*)_';
% create dictionary of shotid mapping to trecvid2013_%d
map = containers.Map;
nframe = length(lstFileName);
for i=1:nframe
	[remat, retok] = regexp(lstFileName{i}, reg, 'match', 'tokens');
	map(retok{end}{2}) = retok{end}{1};
end

for k=1:length(rlist_dir)
	if rlist_dir(k).name(1) == '.'
		continue;
	end
	
	query_path = fullfile(res_dir, rlist_dir(k).name) %'/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/result/best_config/txt/fg+bg_0.1_perdoch_hesaff_rootsift_akmeans_1000000_100000000_50_kdtree_8_800_v1_f1_3_0.0125_avg_pooling_full_notrim_clip_idf_nonorm_kdtree_3_0.0125_-1_dist_avg_autoasym_ivf_0.7/';
	output_dir = ['/net/per610a/export/das11f/ledduy/trecvid-ins-2013/result/nvtiep_No1_KNN5_10K', rlist_dir(k).name ,'/tv2013/test2013-new/']
	key_press = input('Xu ly khong? 0-> khong xu ly, 1 xu ly, -1 thoat');
	if key_press==-1
		break;
	end
	if key_press==0
		continue;
	end

	% load tat ca cac file trong top 1000
	query_res_file = dir([query_path,'/*.txt']);
	%query_res_file = dir([query_path,'/*.list']);
	nquery = length(query_res_file);
	reg1 = '(shot\d*_\d*)\(dist_(.*)\): (.*) (.*)';
	reg2 = '(shot\d*_\d*)\(dist_(.*)\): (.*)';


	for i=1:nquery
		res_file = fullfile(query_path, query_res_file(i).name);
		fid = fopen(res_file);
		str=fgetl(fid);
		% read distance list
		dist_map = containers.Map;
		k = 0;
		while ~feof(fid)
			k = k+1;
			if k > topK
				break;
			end
			str=strtrim(fgetl(fid));
			
			C = textscan(str, '%s %d %s %d %f');
			
			[remat, retok] = regexp(str, reg1, 'match', 'tokens');
			if isempty(retok)
				[remat, retok] = regexp(str, reg2, 'match', 'tokens');
			end
			dist_map(retok{end}{1}) = str2double(retok{end}{2});
			
			%dist_map(C{3}{1}) = C{5};
			
		end
		fclose(fid);
		% convert
		keyset = keys(dist_map);
		nkey = length(keyset);
		if ~exist([output_dir, query_res_file(i).name(1:4)], 'dir')
			mkdir([output_dir, query_res_file(i).name(1:4)]);
        end
        unix(['chmod 777 -R ', output_dir, query_res_file(i).name(1:4)])
		for j=1:nkey
			file = [output_dir, query_res_file(i).name(1:4), '/',map(keyset{j}), '.res'];
			fid=fopen(file, 'a');
			tmp = dist_map(keyset{j});
			fprintf(fid,'%s #$# %s #$# %s\n', keyset{j}, keyset{j}, 2-dist_map(keyset{j}));
			%fprintf(fid,'%s #$# %s #$# %s\n', keyset{j}, keyset{j}, dist_map(keyset{j}));
			fclose(fid);
		end
	end
end