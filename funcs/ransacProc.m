function [ score, good ] = ransacProc(ranks,db_lut,rep_qid,rep_vid,fg_scoring,fg_wei,esc_thre,quant_knn, query_frame_info,database, query_topics,db_frame_dir,frame_sampling,ransac_method,dup,thre,show,midres_filename)

if ~exist('esc_thre','var')
    esc_thre = inf;
end

if ~exist('quant_knn','var')
    quant_knn = 1;
end

if ~exist('db_frame_dir','var')
    show = false;
end

if ~exist('show','var')
    show = false;
end

if fg_scoring == 2
    % for fg_scoring method 2, only counting fg region, so fg_wei is not necessary
    fg_wei = 1;
end

[num_video,num_query] = size(ranks);
num_feat = length(database);
score = zeros(size(ranks));
good = false(size(ranks));
total_time = 0;
for qid = 1:num_query
    prev_good_idx = 0;
    for rid = 1:num_video
        tic;
        fprintf('\r%d/%d %d/%d ',qid,num_query,rid,num_video);
        if rid-prev_good_idx>esc_thre
            fprintf('abort from topic:%s,image:%d,video:%d\n', topic_name, iid, rid);
            break;
        end
        if rep_vid(rid,qid) == 0
            continue;
        end
        rep_qid_kp=cell(1,num_feat);
        rep_vid_kp=cell(1,num_feat);
        rep_qid_kp_wid=cell(1,num_feat);
        rep_vid_kp_wid=cell(1,num_feat);
        rep_qid_fg_idx=cell(1,num_feat);
        for nid = 1:num_feat
            rep_vid_id = 1+(rep_vid(rid,qid)-1)*frame_sampling;
            video_name = db_lut{ranks(rid,qid)};
            video_feat_fullname = fullfile(database{nid}.db_mat_dir,[video_name,'.mat']);
            video_quant_fullname = fullfile(database{nid}.quant_dir,[video_name,'.mat']);
            if ~exist(video_feat_fullname,'file')
                fprintf('feature mat %s doesnot exist\n',video_feat_fullname);
                continue;
            end
            if ~exist(video_quant_fullname,'file')
                fprintf('quant file %s doesnot exist\n',video_quant_fullname);
                continue;
            end
            load(video_feat_fullname, 'clip_kp');
            if rep_vid_id > length(clip_kp)
                fprintf('warning: rep vid id exceed the length of the video\n');
            end
            rep_kps = clip_kp{min(rep_vid_id,end)};
            clear clip_kp
            load(video_quant_fullname, 'bins');
            rep_bins = bins{min(rep_vid_id,end)}(1:min(quant_knn,end),:);
            
            rep_kps = repmat(rep_kps,1,size(rep_bins,1));
            rep_bins = reshape(rep_bins',1, []);
            rep_vid_kp_wid{nid} = double(rep_bins);
            rep_vid_kp{nid} = double(rep_kps);

            rep_kps = query_frame_info{nid}{qid}.query_kp{rep_qid(rid,qid)};
            rep_bins = query_frame_info{nid}{qid}.valid_bins{rep_qid(rid,qid)}(1:min(quant_knn,end),:);
            rep_fg_idx = false(1,size(rep_bins,2));
            rep_fg_idx(query_frame_info{nid}{qid}.fg_index{rep_qid(rid,qid)}) = true;
            rep_kps = repmat(rep_kps,1,size(rep_bins,1));
            rep_fg_idx = repmat(rep_fg_idx,1,size(rep_bins,1));
            rep_bins = reshape(rep_bins',1, []);
            rep_qid_kp_wid{nid} = double(rep_bins);
            rep_qid_kp{nid} = double(rep_kps);
            rep_qid_fg_idx{nid} = rep_fg_idx;
        end
        if show
            rep_qid_filename = query_topics{qid}{1}{rep_qid(rid,qid)};
            load(strrep(video_feat_fullname,'rootsift','sift_color'), 'clip_frame');
            rep_vid_filename = fullfile(db_frame_dir,db_lut{ranks(rid,qid)},[clip_frame{rep_vid(rid,qid)} '.png']);
            clear clip_frame
        else
            rep_qid_filename = [];
            rep_vid_filename = [];
        end
        
        rep_qid_kp=cell2mat(rep_qid_kp);
        rep_vid_kp=cell2mat(rep_vid_kp);
        rep_qid_kp_wid=cell2mat(rep_qid_kp_wid);
        rep_vid_kp_wid=cell2mat(rep_vid_kp_wid);
        rep_qid_fg_idx=cell2mat(cellfun(@(x) logical(x),rep_qid_fg_idx, 'UniformOutput', false));
        
        com_wid = intersect(rep_qid_kp_wid(:),rep_vid_kp_wid(:));
        if isempty(com_wid)
            continue;
        end            
        matches=cell(1,length(com_wid));
        qid_ind = repmat(1:size(rep_qid_kp_wid,2),size(rep_qid_kp_wid,1),1);
        vid_ind = repmat(1:size(rep_vid_kp_wid,2),size(rep_vid_kp_wid,1),1);
        for i=1:length(com_wid)
            com_qid_ind = qid_ind(rep_qid_kp_wid==com_wid(i));
            com_vid_ind = vid_ind(rep_vid_kp_wid==com_wid(i));
            [p,q]=meshgrid(com_qid_ind,com_vid_ind);
            matches{i} = [p(:) q(:)]'; 
        end
        matches=cell2mat(matches);
        if fg_scoring == 1 || fg_scoring == 3
            % only check fg region and discard other matches from the bg rgion
            full_matches = matches;
            matches = matches(:,rep_qid_fg_idx(int32(matches(1,:))));
        end
        if dup<=0 
            X = [rep_qid_kp(1:2,matches(1,:)); rep_vid_kp(1:2,matches(2,:))];
            X = unique(X','rows');
            X = X';
            rep_qid_kp = X(1:2,:);
            rep_vid_kp = X(3:4,:);
        else            
            X = [rep_qid_kp(1:2,matches(1,:)); rep_vid_kp(1:2,matches(2,:))];
            uniq_flag=sparsifypoints(X,dup);
            rep_qid_kp=X(1:2,uniq_flag);
            rep_vid_kp=X(3:4,uniq_flag);
            %fprintf('%d--->%d\n',size(X,2),size(rep_qid_kp,2));
        end
        matches = repmat(1:size(rep_qid_kp,2),2,1);     
        if length(matches)<=4 || length(unique(rep_qid_kp', 'rows'))<=4 ...
                || length(unique(rep_vid_kp', 'rows'))<=4
            continue;
        end
        [H,inliers] = ransacAlgo(rep_qid_kp, rep_vid_kp, ...
            matches,rep_qid_filename, rep_vid_filename, ransac_method, thre, show);
        if isempty(inliers) || length(unique(matches(1,inliers)))<=4 || length(unique(matches(2,inliers)))<=4
            continue;
        end
        
        switch fg_scoring
        case 0
            weight=ones(1,size(rep_qid_kp,2));
            weight(rep_qid_fg_idx) = fg_wei;
        case 1 
            weight=ones(1,size(rep_qid_kp,2));
        case 2
            % only scoring pts within fg region, differ from 1 in that bg points can help to compute H
            weight=zeros(1,size(rep_qid_kp,2));
            weight(rep_qid_fg_idx) = 1;
        case 3 %% something wrong with this option, tobe checked
            % Calculate, in both directions, the transfered points    
            matches = full_matches;
            X1 = [rep_qid_kp(1:2,matches(1,:)); ones(1,length(matches))];
            X2 = [rep_vid_kp(1:2,matches(2,:)); ones(1,length(matches))];
            X1 = normalise2dpts(X1);
            X2 = normalise2dpts(X2);
            FGX1 = X1(:,rep_qid_fg_idx(matches(1,:)));
            FGX2 = X2(:,rep_qid_fg_idx(matches(1,:)));
            H = homography2d(FGX1(:,inliers),FGX2(:,inliers)); % NOTE, this is necessary for LORANSAC, since the H returned by LORANSAC is sort of different.
            Hx1    = H*X1;
            invHx2 = H\X2;
            
            % Normalise so that the homogeneous scale parameter for all coordinates
            % is 1.
            Hx1    = hnormalise(Hx1);
            invHx2 = hnormalise(invHx2); 
            
            d2 = sum((X1-invHx2).^2)  + sum((X2-Hx1).^2);
            THRESHOLD = 0.001;
            inliers = find(abs(d2) < THRESHOLD);    
            
            weight=ones(1,size(rep_qid_kp,2));
            weight(rep_qid_fg_idx) = fg_wei;
        end
        score(rid,qid) = sum(weight(int32(matches(1,inliers))));
        good(rid,qid) = true;
        if exist('midres_filename','var')~=0 && (mod(rid,10)==0 || rid==num_video)
            save(midres_filename,'score','good','rid','qid','-v7.3');
        end
        prev_good_idx = rid;
        total_time = total_time+toc;
        fprintf('score %d %.0fs',score(rid,qid),total_time);
    end
end
fprintf('\r ransac time %.0fs\n',total_time);
end
