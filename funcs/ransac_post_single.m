function [ score, good ] = ransac_post(ranks,query_topics,db_lut,rep_qid,rep_vid,query_feat_info,db_feat_dir,...
    db_quant_dir,method,rerank_topk,esc_thre,quant_knn,query_image_dir,db_frame_dir,frame_sampling, show)

if ~exist('rerank_topk','var')
    rerank_topk = length(ranks{1}{1});
end

if ~exist('esc_thre','var')
    esc_thre = inf;
end

if ~exist('quant_knn','var')
    quant_knn = 1;
end

if ~exist('query_image_dir','var') || ~exist('db_frame_dir','var')
    show = false;
end

if ~exist('frame_sampling','var')
    frame_sampling = 1;
end

if ~exist('show','var')
    show = false;
end

num_query = length(query_topics);
score = zeros(rerank_topk, num_query);
good = false(rerank_topk, num_query);
total_time = 0;
for qid = 1:num_query
    prev_good_idx = 0;
    for rid = 1:rerank_topk
        tic;
        if rep_vid(rid,qid) > 0 && rep_qid(rid,qid) > 0
            rep_vid_id = 1+(rep_vid(rid,qid)-1)*frame_sampling;
            video_name = db_lut{ranks{qid}{1}(rid)};
            video_feat_fullname = fullfile(db_feat_dir,video_name);
            load(video_feat_fullname, 'clip_kp');
            rep_vid_kp = clip_kp{rep_vid_id};
            clear clip_kp
            video_quant_fullname = fullfile(db_quant_dir,video_name);
            load(video_quant_fullname, 'bins');

            rep_vid_kp_wid = bins{rep_vid_id}(1:min(quant_knn,end),:);
            clear clip_kp bins

            rep_qid_kp = query_feat_info{qid}.query_kp{rep_qid(rid,qid)};
            rep_qid_kp_wid = query_feat_info{qid}.valid_bins{rep_qid(rid,qid)}(1:min(quant_knn,end),:);
            if show
                rep_qid_filename = fullfile(query_image_dir,query_topics{qid}{1}{rep_qid(rid,qid)});
                load(video_feat_fullname, 'clip_frame');
                rep_vid_filename = fullfile(db_frame_dir,db_lut{ranks{qid}{1}(rid)},[clip_frame{rep_vid(rid,qid)} '.png']);
                clear clip_frame
            else
                rep_qid_filename = [];
                rep_vid_filename = [];
            end

            com_wid = intersect(rep_qid_kp_wid(:),rep_vid_kp_wid(:));
            matches=[];
            qid_ind = repmat(1:size(rep_qid_kp_wid,2),size(rep_qid_kp_wid,1),1);
            vid_ind = repmat(1:size(rep_vid_kp_wid,2),size(rep_vid_kp_wid,1),1);
            for i=1:length(com_wid)
                com_qid_ind = qid_ind(rep_qid_kp_wid==com_wid(i));
                com_vid_ind = vid_ind(rep_vid_kp_wid==com_wid(i));
                for j=1:length(com_vid_ind)
                    matches=[matches, [reshape(com_qid_ind,1,[]);ones(1,length(com_qid_ind))*com_vid_ind(j)]];
                end
            end

            [score(rid,qid),passed] = ransac_verify(rep_qid_kp, rep_vid_kp, matches, method, rep_qid_filename, rep_vid_filename, show);
        else
            passed = false;
        end
        
        if passed
            prev_good_idx = rid;
            good(rid,qid) = true;
        elseif rid-prev_good_idx>esc_thre
            fprintf('abort from topic:%s,image:%d,video:%d\n', topic_name, iid, rid);
            break;
        end
        total_time = total_time+toc;
        fprintf('\r%d/%d %d/%d score %d %.0fs',qid,num_query,rid,rerank_topk,score(rid,qid),total_time);
    end
end
fprintf('\r ransac time %.0fs\n',total_time);
end
