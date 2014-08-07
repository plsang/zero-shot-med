clear; clc;
% Add libraries and environmental variable
run('/net/per610a/export/das11f/ledduy/plsang/nvtiep/libs/vlfeat-0.9.18/toolbox/vl_setup.m');
addpath(genpath('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/code/funcs'));
addpath(genpath('/net/per610a/export/das11f/ledduy/plsang/nvtiep/funcs'));

% parameter settings
exp_purpose = 'best_config';
DB = 'INS2013';
disp(exp_purpose)
disp(DB)
eval_topN = 10000;
clobber = false;

work_dir = fullfile('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS', DB);
gt_filename = fullfile(work_dir, 'ins.search.qrels.tv13');

%database.comp_sim= struct('query_obj','fg+bg_0.1','feat_detr','-hesaff', 'feat_desc', '-rootsift -noangle',...
%  'clustering','akmeans','K',1000000,'num_samps',100000000,'iter',50,...
%  'build_params',struct('algorithm', 'kdtree','trees', 8, 'checks', 800, 'cores', 10),...
%  'video_sampling',1,'frame_sampling',1,'knn',1,'delta_sqr',6250,'db_agg','avg_pooling',...
%  'vocab','full','trim','notrim','freq','clip','weight','idf','norm','l1',...
%  'query_knn',3,'query_delta_sqr',6250,'query_num',-1,'query_agg','avg_pooling','dist','l2asym_ivf');
database.comp_sim= struct('query_obj','fg+bg_0.1','feat_detr','-perdoch -hesaff', 'feat_desc', '-rootsift',...
    'clustering','akmeans','K',1000000,'num_samps',100000000,'iter',50,...
    'build_params',struct('algorithm', 'kdtree','trees', 8, 'checks', 800, 'cores', 10),...
    'video_sampling',1,'frame_sampling',1,'knn',3,'delta_sqr',6250,'db_agg','avg_pooling',...
    'vocab','full','trim','notrim','freq','clip','weight','idf','norm','nonorm',...
    'query_knn',3,'query_delta_sqr',6250,'query_num',-1,'query_agg','dist_avg','dist','autoasym_ivf_0.7');
	
database.comp_sim= struct('query_obj','fg+bg_0.1','feat_detr','-hesaff', 'feat_desc', '-rootsift -noangle',...
    'clustering','akmeans','K',1000000,'num_samps',100000000,'iter',50,...
    'build_params',struct('algorithm', 'kdtree','trees', 8, 'checks', 800, 'cores', 10),...
    'video_sampling',1,'frame_sampling',1,'knn',1,'delta_sqr',6250,'db_agg','avg_pooling',...
    'vocab','full','trim','notrim','freq','clip','weight','idf','norm','nonorm',...
    'query_knn',3,'query_delta_sqr',6250,'query_num',-1,'query_agg','avg_pooling','dist','l1_ivf');

database.db_frame_dir = '/net/per610a/export/das11g/caizhizhu/ins/ins2013/frames_png'; 
database.query_dir = fullfile(work_dir,'query');
database.query_mask_dir = fullfile(database.query_dir,'masks');
database.query_mat_dir = fullfile(database.query_dir,'result',exp_purpose,'mat');
database.query_txt_dir = fullfile(database.query_dir,'result',exp_purpose,'txt');
database.query_perf_dir = fullfile(database.query_dir,'result',exp_purpose,'perf');

% Prepare dirs
if ~exist(database.query_mat_dir,'dir')
    mkdir(database.query_mat_dir);
end
if ~exist(database.query_txt_dir,'dir')
    mkdir(database.query_txt_dir);
end
if ~exist(database.query_perf_dir,'dir')
    mkdir(database.query_perf_dir);
end

% set proper norm according to the dist and agg options
if strcmp(database.comp_sim.dist, 'l1_ivf')
	database.comp_sim.norm='l1';
	database.comp_sim.delta_sqr=6250;
	database.comp_sim.query_delta_sqr=6250;
elseif strcmp(database.comp_sim.dist, 'l2_ivf') 
	database.comp_sim.norm='l2';
	database.comp_sim.delta_sqr=6250;
	database.comp_sim.query_delta_sqr=6250;
elseif ~isempty(strfind(database.comp_sim.dist, 'asym'))
	database.comp_sim.norm='nonorm';
	database.comp_sim.delta_sqr=6250;
	database.comp_sim.query_delta_sqr=6250;
end

% rescale delta_sqr
if database.comp_sim.knn>1 && database.comp_sim.delta_sqr~=-1
	if ~isempty(strfind(database.comp_sim.feat_desc, 'root'))
		database.comp_sim.delta_sqr=database.comp_sim.delta_sqr/5e5;
	elseif ~isempty(strfind(database.comp_sim.feat_desc, 'color'))
		database.comp_sim.delta_sqr=database.comp_sim.delta_sqr*2;
	elseif ~isempty(strfind(database.comp_sim.feat_desc, 'mom'))
		database.comp_sim.delta_sqr=database.comp_sim.delta_sqr/1e3;
	end
end
if database.comp_sim.query_knn>1 && database.comp_sim.query_delta_sqr~=-1
	if ~isempty(strfind(database.comp_sim.feat_desc, 'root'))
		database.comp_sim.query_delta_sqr=database.comp_sim.query_delta_sqr/5e5;
	elseif ~isempty(strfind(database.comp_sim.feat_desc, 'color'))
		database.comp_sim.query_delta_sqr=database.comp_sim.query_delta_sqr*2;
	elseif ~isempty(strfind(database.comp_sim.feat_desc, 'mom'))
		database.comp_sim.query_delta_sqr=database.comp_sim.query_delta_sqr/1e3;
	end
end
	
%%  RUN RANKING PROGRAM
disp(database.comp_sim);
feature_detr = database.comp_sim.feat_detr;
feature_desc = database.comp_sim.feat_desc;
feature_config = sprintf('%s %s', feature_detr, feature_desc);
feature_name = strrep(feature_config, '-','');
feature_name = strrep(feature_name, ' ','_');
query_feature_name = sprintf('%s_%s',database.comp_sim.query_obj,feature_name);

database.db_mat_dir = fullfile(work_dir,[feature_name '_mat']);

% clustering_name
if ~isempty(strfind(database.comp_sim.clustering,'akmeans'))
	clustering_name = sprintf('%s_%d_%d_%d',database.comp_sim.clustering,...
		database.comp_sim.K,database.comp_sim.num_samps,database.comp_sim.iter); 
end

%build_name
if strcmp(database.comp_sim.build_params.algorithm,'kdtree')
	build_name = sprintf('kdtree_%d_%d',database.comp_sim.build_params.trees,...
		database.comp_sim.build_params.checks);
	quantize_name = sprintf('%d', database.comp_sim.knn);
	if database.comp_sim.knn>1
		quantize_name = sprintf('%s_%g', quantize_name,database.comp_sim.delta_sqr);
	end
end

db_quantize_name = sprintf('v%d_f%d_%s', database.comp_sim.video_sampling,database.comp_sim.frame_sampling,quantize_name);
db_agg_name = database.comp_sim.db_agg;
%bow_making_name
bow_making_name = sprintf('%s_%s_%s_%s_%s',database.comp_sim.vocab,...
	database.comp_sim.trim,database.comp_sim.freq,...
	database.comp_sim.weight,database.comp_sim.norm);
	   
%query_quantize_name
if strcmp(database.comp_sim.build_params.algorithm,'kdtree')
	query_quantize_name = sprintf('kdtree_%d', database.comp_sim.query_knn);
	if database.comp_sim.query_knn>1
		query_quantize_name = sprintf('%s_%g', query_quantize_name,database.comp_sim.query_delta_sqr);
	end
end
% query_agg_name
query_agg_name = sprintf('%d',database.comp_sim.query_num);
query_bow_making_name = sprintf('%s_%d',bow_making_name,database.comp_sim.query_num);
if database.comp_sim.query_num ~= 1
	query_agg_name = sprintf('%d_%s',database.comp_sim.query_num, database.comp_sim.query_agg);
	if ~isempty(strfind(database.comp_sim.query_agg,'avg_pooling'))
		query_bow_making_name = sprintf('%s_%s',bow_making_name,query_agg_name);
	end
end

% cluster_dir
database.cluster_dir = fullfile(work_dir,[feature_name,'_cluster'],clustering_name);
assert(exist(database.cluster_dir,'dir') == 7);
cluster_filename = dir(fullfile(database.cluster_dir,sprintf('Clustering_l2_%d_%d*.hdf5',database.comp_sim.K,database.comp_sim.num_samps)));
assert(length(cluster_filename) == 1);
database.cluster_filename = cluster_filename(1).name;
database.build_dir = fullfile(database.cluster_dir,build_name);
database.bow_dir = fullfile(database.build_dir,db_quantize_name);

% query_frame_dir
database.query_frame_dir = fullfile(database.query_dir,'frames_png');
database.query_feat_dir = fullfile(database.query_dir,'feature', strrep(feature_name,'rootsift','sift'));

% query_quant_dirname: place contains bow of query
query_quant_dirname = sprintf('%s_%s_%s_%s',query_feature_name,clustering_name,build_name,query_quantize_name);
database.query_bow_dir = fullfile(database.query_dir,'bow',query_quant_dirname);
if ~exist(database.query_bow_dir,'dir')
	mkdir(database.query_bow_dir);
end

%dist_name
dist_name = database.comp_sim.dist;

% res_name
res_name = sprintf('%s_%s_%s_%s_%s_%s_%s_%s_%s',query_feature_name,clustering_name,...
	build_name,db_quantize_name,db_agg_name,bow_making_name,query_quantize_name,query_agg_name,dist_name);
res_filename = fullfile(database.query_mat_dir,['res_' res_name, '.mat']);

if exist(res_filename, 'file')
	disp('Loading result file...')
	load(res_filename);
else
	%% compute result at the beginning
	% Load codebook file
	clustering_file = fullfile(database.cluster_dir, database.cluster_filename);
	if ~exist('prev_clustering_file','var') || ~strcmp(clustering_file,prev_clustering_file)
		if ~exist(clustering_file,'file')
			fprintf('centroid file (%s) doesnot exist!\n',clustering_file)
			return;
		end
		time('centers=hdf5read(clustering_file,''/clusters'');');
		[feat_len,hist_len]=size(centers);
		prev_clustering_file = clustering_file;
	end
	% database bow_file
	if strcmp(db_agg_name,'avg_pooling')
		bow_file = fullfile(database.bow_dir,sprintf('bow_clip_%s_%s.mat', bow_making_name, db_agg_name));
	end
	if ~exist('prev_bow_file','var') || ~strcmp(bow_file,prev_bow_file)
		if exist(bow_file,'file') && ~clobber
			time('load(bow_file)','Loading weighted and normalized database bow file...');
		else
			%% Compute database bow at the beginning
			% Load raw bag of word of the database
			raw_bow_file = fullfile(database.bow_dir,[database.comp_sim.db_agg, '_raw_bow.mat']);
			raw_bow_file = fullfile(database.bow_dir,'raw_bow.mat');
			if ~exist(raw_bow_file,'file')
				fprintf('raw bow file (%s) doesnot exist!\n',raw_bow_file)
				return;
			end;
			time('load(raw_bow_file);','Loading raw database frame bow...');
			
			assert(exist('list_frame_bow', 'var') ~= 0);
			if database.comp_sim.frame_sampling > 1
				disp('sampling database frames');
				list_frame_bow = cellfun(@(x) x(:,1:database.comp_sim.frame_sampling:end), list_frame_bow,'UniformOutput', false);
			end
			
			% Pooling
			switch database.comp_sim.db_agg
			case 'max_pooling'
				list_clip_bow = cellfun(@(x) max(x,[],2), list_frame_bow, 'uniformoutput', false);
				clip_frame_num = cellfun( @(x) size(x,2), list_frame_bow,'uniformoutput',false);
				clip_frame_num = cell2mat(clip_frame_num);
				db_bow = sparse([list_clip_bow{:}]);
			case 'avg_pooling'
				list_clip_bow = cellfun(@(x) mean(x,2), list_frame_bow, 'uniformoutput', false);
				clip_frame_num = cellfun( @(x) size(x,2), list_frame_bow,'uniformoutput',false);
				clip_frame_num = cell2mat(clip_frame_num);
				db_bow = sparse([list_clip_bow{:}]);
			otherwise
				disp('error db_agg option!');
				return;
			end
			big_bow_info_file = fullfile(database.bow_dir,'raw_bow_info.mat');
			time('load(big_bow_info_file);','load raw bow info ...');
			
			% computing tf 
			term_freq = list_term_freq.(database.comp_sim.freq);
			db_lut=list_id2clip_lut;
			%% CLEAR redundant data
			clear list_term_freq list_clip_bow list_frame_bow list_avg_pooling_bow list_max_pooling_bow list_id2clip_lut;
			
			% compute weighting
			weight = get_wei(term_freq,database.comp_sim.weight);
			% trim bow
			if ~strcmp(database.comp_sim.trim,'notrim')
				db_bow = trim_bow(db_bow,database.comp_sim.trim);
			end
			
			disp('weighting and normalizing');
			% matlabpool(8)
			% apply weighting
			% normalize bow
			if ~strcmp(database.comp_sim.weight,'nowei')...
				&& ~strcmp(database.comp_sim.norm,'nonorm')
				norm_id=str2double(database.comp_sim.norm(end));
				assert(norm_id == 1 || norm_id == 2);
				for i=1:size(db_bow,2)
					db_bow(:,i) = db_bow(:,i).*weight;
					bow_norm = norm(db_bow(:,i),norm_id)+eps;
					db_bow(:,i) = db_bow(:,i)./bow_norm;
				end
			elseif ~strcmp(database.comp_sim.weight,'nowei')
				for i=1:size(db_bow,2)
					db_bow(:,i) = db_bow(:,i).*weight;
				end
			elseif ~strcmp(database.comp_sim.norm,'nonorm')
				norm_id=str2double(database.comp_sim.norm(end));
				assert(norm_id == 1 || norm_id == 2);
				for i=1:size(db_bow,2)
					bow_norm = norm(db_bow(:,i),norm_id)+eps;
					db_bow(:,i) = db_bow(:,i)./bow_norm;
				end
			end
			%matlabpool close
			time('save(bow_file,''db_bow'',''db_lut'',''weight'',''clip_frame_num'',''-v7.3'')',...
				 'saving weighted and normalized database bow file ...');
		end
		% Build inverted file
		ivf = [];
		if ~isempty(strfind(database.comp_sim.dist,'ivf'))
			time('ivf = BuildInvFile([],db_bow,0,false);','Building inverted file...');
			db_bow = [];
		end
		prev_bow_file = bow_file;
	end
	
	% query bow file
	query_bow_file = fullfile(database.query_bow_dir,['bow_' query_bow_making_name, '.mat']);
	if exist(query_bow_file,'file') && ~clobber
		time('load(query_bow_file)','load weighted and normalized query bow file ...');
	else
		% raw query bow
		query_raw_bow_file = fullfile(database.query_bow_dir,'raw_bow.mat');
		if exist(query_raw_bow_file,'file') && ~clobber
			time('load(query_raw_bow_file)','load raw query bow file ...');
		else
			if strcmp(database.comp_sim.build_params.algorithm,'kdtree')
				kdtree_filename = fullfile(database.build_dir,'flann_kdtree.bin');
				kdsearch_filename = fullfile(database.build_dir,'flann_kdtree_search.mat');
				assert(exist(kdtree_filename,'file')~=0);
				disp('Loading kdtree ...');
				kdtree = flann_load_index(kdtree_filename,single(centers));
				load(kdsearch_filename);
				search_params.cores = database.comp_sim.build_params.cores;
			end
	
			query_dir = dir(database.query_frame_dir);
			query_dir = {query_dir(:).name};
			valid_ids = cellfun(@(x) ~strcmp(x(1),'.'), query_dir,'UniformOutput',false);
			query_dir = query_dir(cell2mat(valid_ids));
			query_num = length(query_dir);
			topic_bows = cell(1, query_num);
			frame_quant_info = cell(1, query_num);
			query_filenames = cell(1,query_num);
			disp('extract query feature and quantization...');
			for qid = 1:query_num
				tic;
				fprintf('\r%2d(1-%d) ',qid,query_num);
				query_pathname = fullfile(database.query_frame_dir,query_dir{qid});
				query_imgs = dir([query_pathname, '/*.png']);
				query_imgs = {query_imgs(:).name};
				query_filenames{qid} = cellfun(@(x) fullfile(query_pathname,x), query_imgs,'UniformOutput',false);
				num_query_imgs = length(query_imgs);
				if num_query_imgs == 0
					error('query images do not exist in %s!\n',query_pathname);
				end;
		
				sift_dir = fullfile(database.query_feat_dir,query_dir{qid});
				if ~exist(sift_dir,'dir')
					mkdir(sift_dir);
				end
				
				%topic_bows{qid}=zeros(length(vocab_range),num_query_imgs);
				topic_bows{qid}=zeros(hist_len,num_query_imgs);
				frame_quant_info{qid}.fg_index = cell(1, num_query_imgs);
				frame_quant_info{qid}.bg_index = cell(1, num_query_imgs);
				frame_quant_info{qid}.query_kp = cell(1, num_query_imgs);
				frame_quant_info{qid}.query_desc = cell(1, num_query_imgs);
				frame_quant_info{qid}.valid_bins = cell(1, num_query_imgs);
				frame_quant_info{qid}.valid_sqrdists = cell(1, num_query_imgs);
		
				% Tiep: total fg+bg feats
				total_fg_feat = 0;
				total_bg_feat = 0;
				for i=1:num_query_imgs
					clear desc kp
					% query sift extraction
					%disp('query sift extraction ...');
					if ~isempty(strfind(feature_name, 'perdoch'))
						[kp,desc] = mxhesaff(query_filenames{qid}{i},~isempty(strfind(feature_name,'root')),false);
					else
						query_sift_filename=fullfile(sift_dir,strrep(query_imgs{i},'png','txt'));
						if ~isempty(strfind(strrep(feature_config,'root',''),'-sift'))
							if ~exist(query_sift_filename,'file')
								exe = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/code/funcs/compute_descriptors_64bit.ln';
								unix(sprintf('%s %s -i %s -o1 %s>/dev/null 2>&1', exe,...
								strrep(feature_config,'rootsift','sift'), query_filenames{qid}{i}, query_sift_filename));
								if ~exist(query_sift_filename,'file')
									%delete(query_filenames{qid}{i});
									fprintf('query has no sift been detected in query image %s!\n',...
									fullfile(query_pathname, query_imgs{i}));
									continue;
								end
							end
							[kp,desc] = vl_ubcread(query_sift_filename, 'format', 'oxford');
						else
							if ~exist(query_sift_filename,'file')
								if ~isempty(strfind(strrep(feature_config,'root',''),'-sc'))
									feature_id = 'sc';
									exe = './funcs/compute_descriptors.ln';
								elseif ~isempty(strfind(strrep(feature_config,'root',''),'-mom'))
									feature_id = 'mom';
									exe = './funcs/compute_descriptors_linux64';
								end
								pt_file = strrep(query_sift_filename,feature_id,'sift');
								unix(sprintf('%s -p1 %s %s -i %s -o1 %s>/dev/null 2>&1', exe, pt_file, ...
								strrep(feature_desc,'rootsift','sift'), query_filenames{qid}{i}, query_sift_filename));
								if ~exist(query_sift_filename,'file')
									%delete(query_filenames{qid}{i});
									fprintf('query has no %s been detected in query image %s!\n',...
									feature_id, fullfile(query_pathname, query_imgs{i}));
									continue;
								end
							end
							[kp,desc] = ubcread_float(query_sift_filename);
						end
						
						if ~isempty(strfind(feature_name,'rootsift'))
							root_sift = zeros(feat_len,size(desc,2));
							for k = 1:size(desc,2)
								sift = double(desc(1:feat_len,k));
								root_sift(:,k) = sift ./ norm(sift,1);
							end
							desc = sqrt(root_sift);
						end
					end
					
					if ~strcmp(database.comp_sim.query_obj,'crop_fg')  %'crop_fg' is a obsolete option
						query_mask_filename=fullfile(database.query_mask_dir,strrep(query_imgs{i},'src','mask'));
						if exist(query_mask_filename,'file')
							mask = imread(query_mask_filename);
							mask = mask(:,:,1)>128;
							SE = strel('square', 5);
							mask=imdilate(mask,SE);
							% Tiep: blur fg and bg
							%G = fspecial('gaussian',[31 31],2);
							%mask = rgb2gray(mask);
							%mask_blur = imfilter(mask,G,'same');
							xy = floor(kp(1:2,:));
							fg_index=find(mask(sub2ind(size(mask),xy(2,:),xy(1,:))));
							%fg_weis =mask_blur(xy(2,fg_index), xy(1,fg_index));
							bg_index=find(mask(sub2ind(size(mask),xy(2,:),xy(1,:)))==0);
							total_fg_feat = total_fg_feat+length(fg_index);
							total_bg_feat = total_bg_feat+length(bg_index);
						else
							disp('Mask files doesnot exist, please check!');
						end
					end
					
					if strcmp(database.comp_sim.query_obj,'fg_bg') ||...
							strcmp(database.comp_sim.query_obj,'crop_fg') ||...
							~isempty(strfind(database.comp_sim.query_obj,'fg+bg')) ||...
							~exist('mask','var')
						query_kp=kp;
						query_desc=desc;
					elseif strcmp(database.comp_sim.query_obj,'fg')
						query_kp=kp(:,fg_index);
						query_desc=desc(:,fg_index);
					elseif strcmp(database.comp_sim.query_obj,'bg')
						query_kp=kp(:,bg_index);
						query_desc=desc(:,bg_index);
					end        
					
					% quantize query sift
					% disp('quantize query sift ...');
					if strcmp(database.comp_sim.build_params.algorithm,'kdtree')
						[bins,dist_sqr] = flann_search(kdtree,single(query_desc),database.comp_sim.query_knn, search_params);
						% fg+bg_0.1: weight 0.1 for background and 1 for foreground
						bg_index=find(mask(sub2ind(size(mask),xy(2,:),xy(1,:)))==0);
						pos=strfind(database.comp_sim.query_obj,'prop');
						if isempty(pos)
							fgbg_base=1;
						else
							fgbg_base=length(fg_index)/length(bg_index);
						end
						pos=strfind(database.comp_sim.query_obj,'_');
						if isempty(pos(end))
							fgbg_rate=1;
						else
							fgbg_rate=str2double(database.comp_sim.query_obj(pos(end)+1:end));
						end
						re_bins = reshape(bins(:,fg_index),1,[]);
						if database.comp_sim.query_delta_sqr ~= -1 
							weis = exp(-dist_sqr(:,fg_index)./(2*database.comp_sim.query_delta_sqr));
							weis = weis./repmat(sum(weis,1),size(weis,1),1);  % philbin, Lost in Quantization
							weis = reshape(weis,1,[]);
						else
							weis = 1;
						end
						% soft assignment on fg region 
						topic_bows{qid}(:,i) = vl_binsum(topic_bows{qid}(:,i),double(weis),double(re_bins));

						pos=strfind(database.comp_sim.query_obj,'bgsoft');
						if isempty(pos)
							weis = fgbg_rate*fgbg_base;
							% hard assignment on bg region
							topic_bows{qid}(:,i) = vl_binsum(topic_bows{qid}(:,i),double(weis),double(bins(1,bg_index)));
						else
							re_bins = reshape(bins(:,bg_index),1,[]);
							if database.comp_sim.query_knn>1 && database.comp_sim.query_delta_sqr ~= -1
								weis = exp(-dist_sqr(:,bg_index)./(2*database.comp_sim.query_delta_sqr));
								weis = weis./repmat(sum(weis,1),size(weis,1),1);  % philbin, Lost in Quantization
								weis = reshape(weis,1,[]);
							else
								weis = 1;
							end
							weis = weis*fgbg_rate*fgbg_base;
							% hard assignment on bg region
							topic_bows{qid}(:,i) = vl_binsum(topic_bows{qid}(:,i),double(weis),double(re_bins));
						end

						frame_quant_info{qid}.fg_index{i} = fg_index;
						frame_quant_info{qid}.bg_index{i} = bg_index;
						frame_quant_info{qid}.query_kp{i} = query_kp;
						frame_quant_info{qid}.query_desc{i} = query_desc;
						frame_quant_info{qid}.valid_bins{i} = bins;
						frame_quant_info{qid}.valid_sqrdists{i} = dist_sqr;
					end
				end
				% Tiep: save #feature of fr and bg
				nfeat_filename=fullfile(database.query_feat_dir,'num.feat.txt');
				fid = fopen(nfeat_filename, 'a');
				fprintf(fid,'%s %d %d\n', query_dir{qid}, total_fg_feat, total_bg_feat);
				fclose(fid);
				
				topic_bows{qid} = sparse(topic_bows{qid});
				fprintf(' %.0f',toc);
			end
			fprintf('\n');
			time('save(query_raw_bow_file,''topic_bows'',''frame_quant_info'',''query_filenames'',''-v7.3'')',...
				'save raw query bow file ...');
		end

		disp('Normalizing raw query bow ...');
		query_num = length(query_filenames);
		for qid = 1:query_num
			tic;
			fprintf('\r%d(1-%d)',qid,query_num);
			nn_set = 1:size(topic_bows{qid},2);
			kk = database.comp_sim.query_num;
			if kk==-1
				kk=10000;
			end
			if length(nn_set)>kk
				query_subsets = nchoosek(nn_set,kk);
			else
				query_subsets = nn_set;
			end
			new_quant_info.fg_index=cell(1,size(query_subsets,1));
			new_quant_info.query_kp=cell(1,size(query_subsets,1));
			new_quant_info.query_desc=cell(1,size(query_subsets,1));
			new_quant_info.quant_bins=cell(1,size(query_subsets,1));
			new_topic_bow = cell(1,size(query_subsets,1));
			new_query_comb = cell(1,size(query_subsets,1));
			query_feat_corr_bridges{qid} = cell(1,size(query_subsets,1));
			for i=1:size(query_subsets,1)
				new_query_comb{i} = query_filenames{qid}(query_subsets(i,:));
				query_bows = topic_bows{qid}(:,query_subsets(i,:));
				new_quant_info.fg_index{i} = frame_quant_info{qid}.fg_index(query_subsets(i,:));
				new_quant_info.quant_bins{i} = frame_quant_info{qid}.valid_bins(query_subsets(i,:));
				new_quant_info.query_kp{i} = frame_quant_info{qid}.query_kp(query_subsets(i,:));
				new_quant_info.query_desc{i} = frame_quant_info{qid}.query_desc(query_subsets(i,:));
				if ~isempty(strfind(database.comp_sim.query_agg,'auto_'))
					pos = strfind(database.comp_sim.query_agg,'_');
					com_word_thre = str2double(database.comp_sim.query_agg(pos(end-1)+1:pos(end)-1));   
					com_word_add_wei = str2double(database.comp_sim.query_agg(pos(end)+1:end));                             
					com_word_num = get_com_word_num(query_bows);
					adj_mat = com_word_num>=com_word_thre | (eye(size(com_word_num))>0);
					com_words = get_com_word_id(query_bows);
					com_word_add_weis = zeros(hist_len,1);
					for m=1:length(com_words)
						for n=m+1:length(com_words)
							if adj_mat(m,n)
								adj_mat(:,m) = adj_mat(:,m) | adj_mat(:,n);
								adj_mat(:,n) = adj_mat(:,m) | adj_mat(:,n);
							end
							lia_m=ismember(com_words{m,n}, new_quant_info.quant_bins{i}{m}(:,new_quant_info.fg_index{i}{m}));
							lia_n=ismember(com_words{m,n}, new_quant_info.quant_bins{i}{n}(:,new_quant_info.fg_index{i}{n}));
							com_word_add_weis(com_words{m,n}(lia_m & lia_n)) = com_word_add_weis(com_words{m,n}(lia_m & lia_n))+com_word_add_wei;
						end
					end   
					if com_word_add_wei ~= 0
						query_bows = query_bows+query_bows.*repmat(com_word_add_weis,1,size(query_bows,2));
					end
					% for this moment only avg_pooling
					full_adj_mat = unique(adj_mat,'rows');  
					new_query_bows = zeros(size(query_bows,1),size(full_adj_mat,1));
					for j = 1:size(full_adj_mat,1)
						new_query_bows(:,j) = mean(query_bows(:,full_adj_mat(j,:)),2);
					end
					query_bows = new_query_bows;  % now the query_bows is shrinked and not concide with query_filenames, may cause problems afterwards.
				end

				if database.comp_sim.query_num ~= 1
					if ~isempty(strfind(database.comp_sim.query_agg,'max_pooling'))
						query_bows = max(query_bows,[],2);
					elseif ~isempty(strfind(database.comp_sim.query_agg,'avg_pooling'))
						query_bows = mean(query_bows,2);
					end
				end
				new_topic_bow{i} = zeros(size(topic_bows{qid},1),size(query_bows,2));
				for k=1:size(query_bows,2)
					% trim bow
					query_bow = trim_bow(query_bows(:,k),database.comp_sim.trim);

					% apply weighting
					query_bow = wei_bow(query_bow,weight);

					% normalize bow
					query_bow = norm_bow(query_bow,database.comp_sim.norm,...
						length(query_filenames{qid}));
					new_topic_bow{i}(:,k) = query_bow;
				end
				new_topic_bow{i}=sparse(new_topic_bow{i});
			end
			topic_bows{qid} = new_topic_bow;
			query_filenames{qid} = new_query_comb;
			frame_quant_info{qid} = new_quant_info;
			fprintf(' %.0f',toc);
		end
		fprintf('\n');
		if isfield(database.comp_sim, 'query_feat_corr') && database.comp_sim.query_feat_corr.bridge_vq
			time('save(query_bow_file,''topic_bows'',''frame_quant_info'',''query_filenames'',''query_feat_corr_bridges'',''-v7.3'')',...
				'save normalized query bow ...');
		else
			time('save(query_bow_file,''topic_bows'',''frame_quant_info'',''query_filenames'',''-v7.3'')',...
				'save normalized query bow ...');
		end
	end
	
	%% Compute distance and ranking
	query_num = length(query_filenames);
	dists = cell(1,query_num);
	disp('Computing distance ...');
	tic;
	for qid = 1:query_num
		fprintf('\r%2d/%d ',qid,query_num);
		subset_num = length(topic_bows{qid});
		dists{qid} = cell(1,subset_num);
		for sid = 1:subset_num
			% make sure comp_dist output all zero distance with all
			% zero queries
			if isfield(database.comp_sim, 'query_feat_corr') && database.comp_sim.query_feat_corr.bridge_vq
				dists{qid}{sid} = comp_dist(ivf,topic_bows{qid}{sid},db_bow,database.comp_sim.dist,false,{query_feat_corr_bridges{qid}{sid}});
			else
				%topic_bows{qid}{sid} = mean(topic_bows{qid}{sid},2);
				dists{qid}{sid} = comp_dist(ivf,topic_bows{qid}{sid},db_bow,database.comp_sim.dist,false);
			end
		end
	end
	fprintf('\n %.4f\n',toc);
	if database.comp_sim.query_num ~= 1 && isempty(strfind(database.comp_sim.query_agg,'max_pooling'))...
			&& isempty(strfind(database.comp_sim.query_agg,'avg_pooling'))
		disp('query aggregation ...');
		for qid = 1:query_num
			tic;
			fprintf('\r%d(1-%d)',qid,query_num);
			subset_num = length(dists{qid});
			for sid = 1:subset_num
				good_id = find(sum(topic_bows{qid}{sid},1));
				num_good = length(good_id);
				if num_good == 0
					dists{qid}{sid} = dists{qid}{sid}(:,1);
					continue;
				elseif num_good == 1
					dists{qid}{sid} = dists{qid}{sid}(:,good_id);
					continue;
				end
				dist = dists{qid}{sid}(:,good_id);
				if ~isempty(strfind(database.comp_sim.query_agg,'rank'))
					[~,idx]  = sort(dist,1,'ascend');
					[~,dist]  = sort(idx,1,'ascend');
				end
				
				if ~isempty(strfind(database.comp_sim.query_agg,'min'))
					[dist,idx] = min(dist,[],2);
					%print out stat to check if mins are equally
					fprintf('min distribution:');
					true_idx = find(dist~=max(dist));
					idx=idx(true_idx);
					for n = 1:num_good
						fprintf(' %d(%.0f%%);',n,length(find(idx==n))*100/length(idx));
					end
					fprintf('\n');
				elseif ~isempty(strfind(database.comp_sim.query_agg,'avg'))
					if ~isempty(strfind(database.comp_sim.query_agg,'fgwei'))
						fg_pt_num = cell2mat(cellfun(@(x)length(x),frame_quant_info{qid}.fg_index{sid},'UniformOutput', false));
						dist_weight = fg_pt_num/sum(fg_pt_num);
						dist = mean(dist.*repmat(full(dist_weight),size(dist,1),1),2);
					elseif ~isempty(strfind(database.comp_sim.query_agg,'wei'))
						pt_num = cell2mat(cellfun(@(x)size(x,2),frame_quant_info{qid}.quant_bins{sid},'UniformOutput', false));
						dist_weight = pt_num/sum(pt_num);
						dist = mean(dist.*repmat(full(dist_weight),size(dist,1),1),2);
					else
						dist = mean(dist,2);
					end
				elseif ~isempty(strfind(database.comp_sim.query_agg,'max'))
					[dist,idx] = max(dist,[],2);
					%print out stat to check if maxs are equally
					fprintf('max distribution:');
					true_idx = find(dist~=min(dist));
					idx=idx(true_idx);
					for n = 1:num_good
						fprintf(' %d(%.0f%%);',n,length(find(idx==n))*100/length(idx));
					end
					fprintf('\n');
				end
				dists{qid}{sid} = dist;
			end
			fprintf(' %.0f',toc);
		end
		fprintf('\n');
	end
	
	% Sort result
	tic;
	fprintf('sorting ...');
	score = cell(1,query_num);
	ranks = cell(1,query_num);
	for qid = 1:query_num
		subset_num = length(topic_bows{qid});
		score{qid} = cell(1,subset_num);
		ranks{qid} = cell(1,subset_num);
		for sid = 1:subset_num
			[score{qid}{sid},ranks{qid}{sid}]=sort(dists{qid}{sid},1);
		end
	end
	fprintf(' %.0f\n',toc);

	% Save result
	time('save(res_filename,''db_lut'',''score'',''ranks'',''query_filenames'',''-v7.3'')');
end

tic;
disp('write knn result of current run ...');
knn_txt_dir = fullfile(database.query_txt_dir,res_name);
if ~exist(knn_txt_dir,'dir')
	mkdir(knn_txt_dir)
end
num_shown_frames = 4;
%write_knn(query_filenames, db_lut, score, ranks, database.db_frame_dir, ...
%	knn_txt_dir, 1000, 2, false);
write_knn(query_filenames, db_lut, score, ranks, database.db_frame_dir, ...
    knn_txt_dir, eval_topN, num_shown_frames, false);
fprintf(' %.0fs\n',toc);

% Save trec performance
res_names = dir(res_filename);
run_name = res_names(1).name(5:end-4);
trecvid_res_dir = fullfile(database.query_perf_dir, run_name);
if ~exist(trecvid_res_dir, 'dir')
	mkdir(trecvid_res_dir);
end

ex_filename = fullfile(work_dir,'dropped.example.image.shots');
fid = fopen(ex_filename,'r');
ex_list = textscan(fid,'%s');
ex_list = ex_list{1};
fclose(fid);

save_trec_perf_list(query_filenames,score,ranks,db_lut,run_name,trecvid_res_dir,eval_topN, ex_list);

% Evaluation
%[gt,ngt] = read_gt_file(gt_filename);
fgt = fopen(gt_filename,'r');
gt_cell = textscan(fgt,'%d %*d %s %d');
fclose(fgt);
[qid, qid_ind,~]=unique(gt_cell{1});
for i=1:length(qid)
    if i==1
        nz_ids = find(gt_cell{3}(1:qid_ind(i)));
        z_ids = find(gt_cell{3}(1:qid_ind(i))==0);
        check_list = gt_cell{2}(1:qid_ind(i));
    else
        nz_ids = find(gt_cell{3}(qid_ind(i-1)+1:qid_ind(i)));
        z_ids = find(gt_cell{3}(qid_ind(i-1)+1:qid_ind(i))==0);
        check_list = gt_cell{2}(qid_ind(i-1)+1:qid_ind(i));
    end
    gt.(['id_' num2str(qid(i))]) = check_list(nz_ids);
    ngt.(['id_' num2str(qid(i))]) = check_list(z_ids);
end

%[gt,ngt] = read_gt_file(gt_filename);
%trecvid_res_dir = '/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/result/best_config/perf/fg+bg_0.1_hesaff_rootsift_noangle_akmeans_1000000_100000000_50_kdtree_8_800_v1_f1_3_0.0125_avg_pooling_full_notrim_clip_idf_nonorm_kdtree_3_0.0125_-1_dist_avg_autoasym_ivf_0.5';
%eval_topN = 1000;
perf = compute_trec_performance(query_filenames,gt_filename,trecvid_res_dir, eval_topN, gt);
fprintf('Performance of the system is: %f', mean([perf{:}]));
