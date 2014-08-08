function sift_select_features( sift_algo, param, version )
%SELECT_FEATURES Summary of this function goes here
%   Detailed explanation goes here
	% nSize: step for dense sift
    % parameters
	
	%%
	if ~exist('version', 'var'),
		version = 'v1.0';
	end

	set_env;
	
	proj_dir = '/net/per610a/export/das11f/plsang/trecvidmed14';
    max_features = 100000000; % 100M
	video_sampling_rate = 0.5;
    sample_length = 25; % frames
    ensure_coef = 1.1;
	
	configs = set_global_config();
	logfile = sprintf('%s/%s.log', configs.logdir, mfilename);
	msg = sprintf('Start running %s(%s, %s)', mfilename, sift_algo, param);
	logmsg(logfile, msg);
	tic;
	
    kf_dir = '/net/per610a/export/das11f/plsang/trecvidmed13/keyframes';
		
	fprintf('Loading metadata...\n');
	medmd_file = '/net/per610a/export/das11f/plsang/trecvidmed14/metadata/medmd_2014.mat';
	MEDMD = load(medmd_file, 'MEDMD'); 
	MEDMD = MEDMD.MEDMD;
	
	clips = MEDMD.Research.default.clips;
	list_video = unique(clips);	% 10161 clips, ~314 hours
	
	num_selected_videos = ceil(video_sampling_rate * length( list_video ));
	rand_index = randperm(length(list_video));
	selected_index = rand_index(1:num_selected_videos);
    selected_videos = list_video(selected_index);
	
	max_features_per_video = ceil(ensure_coef * max_features/length(selected_videos));
    
    data = cell(length(selected_videos), 1);
	location = cell(length(selected_videos), 1);

	output_file = sprintf('%s/feature/bow.codebook.devel/%s.%s.%s.sift/data/selected_feats_%d.hdf5', proj_dir, sift_algo, num2str(param), version, max_features);
	output_location_file = sprintf('%s/feature/bow.codebook.devel/%s.%s.%s.sift/data/selected_location_%d.mat', proj_dir, sift_algo, num2str(param), version, max_features);
	if exist(output_file),
		fprintf('File [%s] already exist. Skipped\n', output_file);
		return;
	end
	
	%matlabpool open;
    for ii = 1:length(selected_videos),
        video_name = selected_videos{ii};
        
		if isfield(MEDMD.lookup, video_name),
			video_kf_dir = fullfile(kf_dir, MEDMD.lookup.(video_name));
			video_kf_dir = video_kf_dir(1:end-4);							% remove .mp4
		else	% err: Reference to non-existent field 'HVC257099'.
			warning('File [%s] does not exist\n',  video_name);
			continue;
		end
		
		kfs = dir([video_kf_dir, '/*.jpg']);
        
		selected_idx = [1:length(kfs)];
		if length(kfs) > sample_length,
			rand_idx = randperm(length(kfs));
			selected_idx = selected_idx(rand_idx(1:sample_length));
		end
		
		fprintf('Computing features for: %d - %s %f %% complete\n', ii, video_name, ii/length(selected_videos)*100.00);
		feat = [];
		frame = [];
		for jj = selected_idx,
			img_name = kfs(jj).name;
			img_path = fullfile(video_kf_dir, img_name);
			
			[frames, descrs] = sift_extract_features( img_path, sift_algo, param );
            
            % if more than 50% of points are empty --> possibley empty image
            if isempty(descrs) || sum(all(descrs == 0, 1)) > 0.5*size(descrs, 2),
                %warning('Maybe blank image...[%s]. Skipped!\n', img_name);
                continue;
            end
			feat = [feat descrs];
			frame = [frame frames];
		end
        
		%size(feat, 2)
		
        if size(feat, 2) > max_features_per_video,
            [data{ii}, sel] = vl_colsubset(feat, max_features_per_video);
			location{ii} = frame(:, sel);
        else
            data{ii} = feat;
			location{ii} = frame;
        end
        
    end
    
	%matlabpool close;
	
    % concatenate features into a single matrix
    data = cat(2, data{:});
	location = cat(2, location{:});
    
    if size(data, 2) > max_features,
         data = vl_colsubset(data, max_features);
		 location = vl_colsubset(location, max_features);
    end

	output_dir = fileparts(output_file);
	if ~exist(output_dir, 'file'),
		mkdir(output_dir);
	end
	
	fprintf('Saving selected features to [%s]...\n', output_file);
	hdf5write(output_file, '/data', data);
	fprintf('Saving selected location to [%s]...\n', output_file);
    save(output_location_file, 'location', '-v7.3');
    
	elapsed = toc;
	elapsed_str = datestr(datenum(0,0,0,0,0,elapsed),'HH:MM:SS');
	msg = sprintf('Finish running %s(%s, %s). Elapsed time: %s', mfilename, sift_algo, param, elapsed_str);
	logmsg(logfile, msg);
end
