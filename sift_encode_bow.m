function sift_encode_bow( proj_name, exp_ann, sift_algo, param, start_seg, end_seg )
%ENCODE Summary of this function goes here
%   Detailed explanation goes here
%% kf_dir_name: name of keyframe folder, e.g. keyframe-60 for segment length of 60s   

	% update: Jun 25th, SPM suported
    % setting
    set_env;
	
	if ~exist('version', 'var'),
		version = 'v1.1';  %% using both event video + bg video
	end
	
    % encoding type
    enc_type = 'bow';
	
	if ~exist('codebook_size', 'var'),
		codebook_size = 1000000;
	end
	
	configs = set_global_config();
	logfile = sprintf('%s/%s.log', configs.logdir, mfilename);
	msg = sprintf('Start running %s(%s, %s, %s, %s, %d, %d)', mfilename, proj_name, exp_ann, sift_algo, param, start_seg, end_seg);
	logmsg(logfile, msg);
	change_perm(logfile);
	tic;
	

	feat_pat = sprintf('%s.%s.%s.sift', sift_algo, num2str(param), version);
	feature_ext = sprintf('%s.cb%d.%s', feat_pat, codebook_size, enc_type);
	
	output_dir = sprintf('/net/per610a/export/das11f/plsang/%s/feature/%s/%s', proj_name, exp_ann, feature_ext);
    if ~exist(output_dir, 'file'),
		mkdir(output_dir);
		change_perm(output_dir);
    end
    
    codebook_file = '/net/per610a/export/das11f/plsang/trecvidmed14/feature/bow.codebook.devel/covdet.hessian.v1.1.sift/data/codebook.akm.1000000.50.hdf5';
	
	fprintf('Loading codebook [%s]...\n', codebook_file);
	centers = hdf5read(codebook_file,'/clusters');
	dataset = single(centers);
	[feat_len, hist_len] = size(dataset);
	fprintf('Deduced cluster center info %d %d...\n', feat_len, hist_len);
		
	
	database.comp_sim.build_params = struct('algorithm', 'kdtree','trees', 8, 'checks', 800, 'cores', 8);
	database.comp_sim.knn = 1;
	database.build_dir = '/net/per610a/export/das11f/plsang/trecvidmed14/feature/bow.codebook.devel/covdet.hessian.v1.1.sift/data/';
	
	quant_struct = struct('quantize','kdtree','build_params', database.comp_sim.build_params, 'knn', database.comp_sim.knn);
	
	%% Build flann
	if ~isempty(strfind(quant_struct.quantize, 'kdtree'))
		kdtree_filename = fullfile(database.build_dir,'flann_kdtree.bin');
		kdsearch_filename = fullfile(database.build_dir,'flann_kdtree_search.mat');
		
		if exist(kdtree_filename,'file')
			fprintf('Loading kdtree ...');
			tic;
			kdtree = flann_load_index(kdtree_filename,dataset);
			load(kdsearch_filename);
			search_params.cores = quant_struct.build_params.cores;
			kdtree_time.load = toc;
		else
			tic;
			fprintf('Building kdtree ...');
			[kdtree,search_params,speedup] = flann_build_index(dataset,quant_struct.build_params); 
			kdtree_time.speedup = speedup;
			kdtree_time.build = toc;
			fprintf('%.0f \n',kdtree_time.build);
			fprintf('save kdtree ...');
			tic;
			flann_save_index(kdtree,kdtree_filename);
			save(kdsearch_filename,'search_params');
			kdtree_time.save = toc;
		end
	end
	
	fprintf('Loading metadata...\n');
	medmd_file = '/net/per610a/export/das11f/plsang/trecvidmed14/metadata/medmd_2014.mat';
	load(medmd_file, 'MEDMD'); 
	
	clips = MEDMD.RefTest.KINDREDTEST.clips;
    
    if ~exist('start_seg', 'var') || start_seg < 1,
        start_seg = 1;
    end
    
    if ~exist('end_seg', 'var') || end_seg > length(clips),
        end_seg = length(clips);
    end
    
    %tic
	
    kf_dir = sprintf('/net/per610a/export/das11f/plsang/%s/keyframes', proj_name);
    
    for ii = start_seg:end_seg,
        video_id = clips{ii};                 
        
		if ~isfield(MEDMD.lookup, video_id),
			msg = sprintf('Unknown location of video <%s>\n', video_id);
			logmsg(logfile, msg);
			continue;
		end
		
		output_file = sprintf('%s/%s/%s.mat', output_dir, fileparts(MEDMD.lookup.(video_id)), video_id);
		
        if exist(output_file, 'file'),
            fprintf('File [%s] already exist. Skipped!!\n', output_file);
            continue;
        end
        
		video_kf_dir = fullfile(kf_dir, MEDMD.lookup.(video_id));
		video_kf_dir = video_kf_dir(1:end-4);
		kfs = dir([video_kf_dir, '/*.jpg']);
       
		%% update Jul 5, 2013: support segment-based
		
		fprintf(' [%d --> %d --> %d] Extracting & encoding for [%s - %d kfs]...\n', start_seg, ii, end_seg, video_id, length(kfs));
        
		code = cell(length(kfs), 1);
		
		for jj = 1:length(kfs),
			if ~mod(jj, 10),
				fprintf('%d ', jj);
			end
			img_name = kfs(jj).name;
			img_path = fullfile(video_kf_dir, img_name);
			
			[frames, descrs] = sift_extract_features( img_path, sift_algo, param );
			
            % if more than 50% of points are empty --> possibley empty image
            if isempty(descrs) || sum(all(descrs == 0, 1)) > 0.5*size(descrs, 2),
                %warning('Maybe blank image...[%s]. Skipped!\n', img_name);
                continue;
            end
			
			bin = flann_search(kdtree,single(descrs),quant_struct.knn, search_params);
			
			bin = reshape(bin(1:quant_struct.knn,:),1,[]);
			weis = ones(size(bin));
			code_ = vl_binsum(zeros(hist_len,1),double(weis),double(bin)); 
			code{jj} = code_;	
		end 
        
		code = cat(2, code{:});
		code = mean(code, 2);
		
        par_save(output_file, code, 1); % MATLAB don't allow to save inside parfor loop             
		%change_perm(output_file);
        
    end
    
	elapsed = toc;
	elapsed_str = datestr(datenum(0,0,0,0,0,elapsed),'HH:MM:SS');
	msg = sprintf('Finish running %s(%s, %s, %s, %s, %d, %d). Elapsed time: %s', mfilename, proj_name, exp_ann, sift_algo, param, start_seg, end_seg, elapsed_str);
	logmsg(logfile, msg);
	
    %toc
    quit;

end