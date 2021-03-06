function [ score, good ] = ransacProcGE(merged_ranks,db_lut,rep_ids,fg_scoring,fg_wei,tune_fg_wei,esc_thre, quant_knn, query_frame_info,database, query_topics,db_frame_dir,frame_sampling,ransac_method,dup,thre,show,midres_filename)

if ~exist('esc_thre','var')
    esc_thre = inf;
end

if quant_knn>0
    VQ4paring = true;
else
    VQ4paring = false;
    %minmax_install();
end

if ~exist('db_frame_dir','var')
    show = false;
end

if ~exist('show','var')
    show = false;
end

[num_video,num_query] = size(merged_ranks);
num_feat = length(database);
if tune_fg_wei
    assert(~(isinf(fg_wei)&&fg_method==1)); %as we will use lazy computing, not correct for weight tuning.
    fg_weis = unique([0,0.01,0.05,0.1,0.2,0.5,0.8,1,2,5,8,10,10000,Inf,fg_wei]);
else
    fg_weis = fg_wei;
end
sum_big_score = cell(length(fg_weis));
for i=1:length(fg_weis)
    sum_big_score{i} = zeros(size(merged_ranks));
end
if rep_ids.rep_qid==-1 || rep_ids.rep_vid==-1
    max_big_score = sum_big_score;
end
good = false(size(merged_ranks));
run_idx = find(fg_weis == fg_wei);

total_time = 0;
for nid = 1:num_feat
    for qid = 1:num_query
        prev_good_idx = 0;
        for rid = 1:num_video
            tic;
            fprintf('\r%d/%d %d/%d %d/%d ',nid, num_feat,qid,num_query,rid,num_video);
            if rid-prev_good_idx>esc_thre
                fprintf('abort from topic:%s,image:%d,video:%d\n', topic_name, iid, rid);
                break;
            end
                
            video_name = db_lut{merged_ranks(rid,qid)};
            video_feat_fullname = fullfile(database{nid}.db_mat_dir,[video_name,'.mat']);
            if ~exist(video_feat_fullname,'file')
                fprintf('feature mat %s doesnot exist\n',video_feat_fullname);
                continue;
            end
            load(video_feat_fullname, 'clip_kp');
            if length(rep_ids.rep_vid)>1
                if rep_ids.rep_vid(rid,qid) == 0
                    continue;
                end
                rep_vids = 1+(rep_ids.rep_vid(rid,qid)-1)*frame_sampling;
                if rep_vids > length(clip_kp)
                    fprintf('warning: rep vid id exceed the length of the video\n');
                end
                rep_kps = {clip_kp{min(rep_vids,end)}};
            elseif rep_ids.rep_vid == 0
                rep_vids = 1:frame_sampling:length(clip_kp);
                empty_flag = cellfun(@(x) ~isempty(x), clip_kp(rep_vids), 'UniformOutput', false);
                rep_vids = rep_vids(cell2mat(empty_flag));
                rand_id = randperm(length(rep_vids));
                rep_vids = rep_vids(rand_id(1));
                rep_kps = {clip_kp{rep_vids}};
            elseif rep_ids.rep_vid == -1 
                rep_vids = 1:frame_sampling:length(clip_kp);
                empty_flag = cellfun(@(x) ~isempty(x), clip_kp(rep_vids), 'UniformOutput', false);
                rep_vids = rep_vids(cell2mat(empty_flag));
                rep_kps = clip_kp(rep_vids);
            end
            clear clip_kp
            if VQ4paring
                video_quant_fullname = fullfile(database{nid}.quant_dir,[video_name,'.mat']);
                if ~exist(video_quant_fullname,'file')
                    fprintf('quant file %s doesnot exist\n',video_quant_fullname);
                    continue;
                end
                load(video_quant_fullname, 'bins');
                rep_bins = cellfun(@(x) double(reshape(x(1:quant_knn,:)',1,[])), bins(rep_vids), 'UniformOutput', false);
                rep_kps = cellfun(@(x) double(repmat(x,1,quant_knn)), rep_kps, 'UniformOutput', false);
                rep_vid_kp_wid = rep_bins;
                rep_vid_kp = rep_kps;
            else
                load(video_feat_fullname, 'clip_desc');
                video_kp = rep_kps;
                video_desc = clip_desc(rep_vids);
            end
            
            if length(rep_ids.rep_qid)>1
                rep_qids = rep_ids.rep_qid(rid,qid);
            elseif rep_ids.rep_qid == 0 
                rep_qids = randperm(length(query_frame_info{nid}{qid}.query_kp));
                rep_qids = rep_qids(1);
            elseif rep_ids.rep_qid == -1
                rep_qids = 1:length(query_frame_info{nid}{qid}.query_kp);
            end
            rep_qid_fg_idx = cell(size(rep_qids));
            for j=1:length(rep_qids)
                i = rep_qids(j);
                rep_fg_idx = false(1,size(query_frame_info{nid}{qid}.valid_bins{i},2));
                rep_fg_idx(query_frame_info{nid}{qid}.fg_index{i}) = true;
                rep_qid_fg_idx{j} = rep_fg_idx;
            end
            if VQ4paring
                rep_qid_kp_wid = cell(size(rep_qids));
                rep_qid_kp = cell(size(rep_qids));
                for j=1:length(rep_qids)
                    i = rep_qids(j);
                    rep_kps = query_frame_info{nid}{qid}.query_kp{i};
                    rep_bins = query_frame_info{nid}{qid}.valid_bins{i}(1:quant_knn,:);
                    rep_kps = repmat(rep_kps,1,size(rep_bins,1));
                    rep_qid_fg_idx{j} = repmat(rep_qid_fg_idx{j},1,size(rep_bins,1));
                    rep_bins = reshape(rep_bins',1, []);
                    rep_qid_kp_wid{j} = double(rep_bins);
                    rep_qid_kp{j} = double(rep_kps);
                end
            else
                query_kp = query_frame_info{nid}{qid}.query_kp(rep_qids);
                query_desc = query_frame_info{nid}{qid}.query_desc(rep_qids);
            end
            
            num_sel_query_img = length(rep_qids);
            num_sel_video_img = length(rep_vids);
            for q=1:num_sel_query_img
                for v=1:num_sel_video_img
                    if VQ4paring
                        com_wid = intersect(rep_qid_kp_wid{q}(:),rep_vid_kp_wid{v}(:));
                        if isempty(com_wid)
                            continue;
                        end
                        matches=cell(1,length(com_wid));
                        qid_ind = repmat(1:size(rep_qid_kp_wid{q},2),size(rep_qid_kp_wid{q},1),1);
                        vid_ind = repmat(1:size(rep_vid_kp_wid{v},2),size(rep_vid_kp_wid{v},1),1);
                        for i=1:length(com_wid)
                            com_qid_ind = qid_ind(rep_qid_kp_wid{q}==com_wid(i));
                            com_vid_ind = vid_ind(rep_vid_kp_wid{v}==com_wid(i));
                            [m,n]=meshgrid(com_qid_ind,com_vid_ind);
                            matches{i} = [m(:) n(:)]';
                        end
                        matches=cell2mat(matches);
                        matched_fg_idx = rep_qid_fg_idx{q}(int32(matches(1,:)));
                        if fg_scoring == 1
                            % in this case, for ransac H computing, only use fg region and discard other matches from the bg rgion
                            if ~isinf(fg_wei)
                                full_X = [rep_qid_kp{q}(1:2,matches(1,:)); rep_vid_kp{v}(1:2,matches(2,:))];
                            end
                            matches = matches(:,matched_fg_idx);
                        end
                        X = [rep_qid_kp{q}(1:2,matches(1,:)); rep_vid_kp{v}(1:2,matches(2,:))];
                    else
                        if fg_scoring == 1 && isinf(fg_wei)
                            sel_query_desc = query_desc{q}(:,rep_qid_fg_idx{q});
                            sel_query_kp = query_kp{q}(:,rep_qid_fg_idx{q});
                        else
                            sel_query_desc = query_desc{q};
                            sel_query_kp = query_kp{q};
                        end
                        dist=vl_alldist2(double(video_desc{v}),sel_query_desc,'l2');
                        [res,loc]=mink(dist,2);
                        % ratio testing
                        matches = res(1,:)<res(2,:)*0.9;
                        matched_fg_idx = false(size(matches));
                        matched_fg_idx(rep_qid_fg_idx{q})=matches(rep_qid_fg_idx{q});
                        if fg_scoring == 1 
                            if ~isinf(fg_wei)
                                full_X = [sel_query_kp(1:2,matches); video_kp{v}(1:2,loc(matches))];
                            end
                            matches = matched_fg_idx;
                        end
                        X = [sel_query_kp(1:2,matches); video_kp{v}(1:2,loc(matches))];
                    end
                    if fg_scoring == 0 || isinf(fg_wei)
                        full_X = X;
                        if fg_scoring == 1
                            matched_fg_idx = true(1,size(full_X,2));%as full_X are all fg pts
                        end
                    end
                    % remove duplicates
                    if dup<=0
                        X = unique(X','rows');
                        X = X';
                    else
                        uniq_flag=sparsifypoints(X,dup);
                        %fprintf('%d--->%d\n',size(X,2),sum(true(uniq_flag)));
                        X=X(:,uniq_flag);
                    end
                    matches = repmat(1:size(X,2),2,1);
                    if size(matches,2)<=4 || length(unique(X(1:2,:)', 'rows'))<=4 ...
                            || length(unique(X(3:4,:)', 'rows'))<=4
                        continue;
                    end
                    
                    if show
                        rep_qid_filename = query_topics{qid}{1}{rep_qids(q)};
                        %load(strrep(video_feat_fullname,'rootsift','sift_color'), 'clip_frame');
                        load(video_feat_fullname, 'clip_frame');
                        rep_vid_filename = fullfile(db_frame_dir,db_lut{merged_ranks(rid,qid)},[clip_frame{rep_vids(v)} '.png']);
                        clear clip_frame
                    else
                        rep_qid_filename = [];
                        rep_vid_filename = [];
                    end
                    [H,inliers] = ransacAlgo(X(1:2,:), X(3:4,:), ...
                        matches,rep_qid_filename, rep_vid_filename, ransac_method, thre, show);
                    if sum(sum(abs(H))) < 0.00001 || isempty(inliers) ||...
                            length(unique(X(1:2,unique(matches(1,inliers)))', 'rows'))<=4 || ...
                            length(unique(X(3:4,unique(matches(2,inliers)))', 'rows'))<=4
                        continue;
                    end
                    
                    matches = repmat(1:size(full_X,2),2,1);
                    inliers_old = inliers;
                    inliers = homogdist2d(H, full_X, thre);
                    if show
                        fprintf('inliers added: %d\n',length(find(inliers))-length(find(inliers_old)));
                    end
                            
                    for i=1:length(fg_weis)
                        this_fg_wei = fg_weis(i);
                        if isinf(this_fg_wei)
                            % only scoring pts within fg region
                            weight=zeros(1,size(full_X,2));
                            weight(matched_fg_idx) = 1;
                        else
                            weight=ones(1,size(full_X,2));
                            weight(matched_fg_idx) = this_fg_wei;
                        end
                        if rep_ids.rep_qid==-1 || rep_ids.rep_vid==-1
                            max_big_score{i}(rid,qid) = max(max_big_score{i}(rid,qid),sum(weight(int32(matches(1,inliers)))));
                        end
                        sum_big_score{i}(rid,qid) = sum_big_score{i}(rid,qid)+sum(weight(int32(matches(1,inliers)))); %should be better than max
                        good(rid,qid) = true;
%                        if ~tune_fg_wei && exist('midres_filename','var')~=0 && (mod(rid,100)==0 || rid==num_video)
%                            this_midres_filename = strrep(midres_filename,'.mat','_tmp.mat');
%                            save(this_midres_filename,'score','good','rid','qid','-v7.3');
%                        end
                    end
                    prev_good_idx = rid;
                end
            end
            total_time = total_time+toc;
            fprintf('sum score %d %.0fs', sum_big_score{run_idx}(rid,qid),total_time);
        end
    end
end
fprintf('\r ransac time %.0fs\n',total_time);
for i=1:length(fg_weis)
    this_fg_wei = fg_weis(i);
    if this_fg_wei == fg_wei
        score = sum_big_score{i};
    else
        new_score = sum_big_score{i};
        good_record = good;
        this_midres_filename = strrep(midres_filename,sprintf('_w%g',fg_wei),sprintf('_w%g',this_fg_wei));
        save(this_midres_filename,'db_lut','query_topics','new_score','good_record','merged_ranks','-v7.3');
    end
end
if rep_ids.rep_qid==-1 || rep_ids.rep_vid==-1
    for i=1:length(fg_weis)
        this_fg_wei = fg_weis(i);
        new_score = max_big_score{i};
        good_record = good;
        if length(rep_ids.rep_qid)>1
            subfix = '_q1_'; 
        elseif rep_ids.rep_qid == 0 
            subfix = '_q0_';
        elseif rep_ids.rep_qid == -1
            subfix = '_q-1_';
        end
        this_midres_filename = strrep(strrep(midres_filename, sprintf('_w%g',fg_wei),...
           sprintf('_w%g',this_fg_wei)),subfix,[subfix 'max_']);
        save(this_midres_filename,'db_lut','query_topics','new_score','good_record','merged_ranks','-v7.3');
    end
end
end


function [inliers, H] = homogdist2d(H, x, t)
    
    x1 = [x(1:2,:); ones(1,size(x,2))];   % Extract x1 and x2 from x
    x2 = [x(3:4,:); ones(1,size(x,2))];    
    
    % Calculate, in both directions, the transfered points    
    Hx1    = H*x1;
    invHx2 = H\x2;
    
    % Normalise so that the homogeneous scale parameter for all coordinates
    % is 1.
    
    x1     = hnormalise(x1);
    x2     = hnormalise(x2);     
    Hx1    = hnormalise(Hx1);
    invHx2 = hnormalise(invHx2); 
    
    inliers = find((sqrt(sum((x1-invHx2).^2)) < t) & (sqrt(sum((x2-Hx1).^2)) < t));  
end
