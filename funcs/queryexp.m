function [ new_score, new_ranks ] = queryexp(topic_bows, db_bow, score, ranks, avg_topn, norm_type, dist_type, exp_type,is_ransac_res,good_dist_thre)  
% a simple implementation of the average query expansion described in
% Total Recall: Automatic Query Expansion with a Generative Feature Model
% for Object Retrieval

query_num = length(topic_bows);
avg_topic_bows = topic_bows;

assert(~isempty(strfind(dist_type,'ivf')));
fprintf('build invfile...');
tic;
ivf = BuildInvFile([],db_bow,0,false);
fprintf('%.0fs\n',toc);
dists = cell(1,query_num);
new_score = cell(1,query_num);
new_ranks = cell(1,query_num);
total_time = 0;
disp('query expansion and invfile based search...');
for qid = 1:query_num
    tic;
    subset_num = length(topic_bows{qid});
    dists{qid} = cell(1,subset_num);
    new_score{qid} = cell(1,subset_num);
    new_ranks{qid} = cell(1,subset_num);
    for sid = 1:subset_num
        if exist('good_dist_thre','var')
            good_dist_idx = score{qid}{sid}<good_dist_thre;
            avg_idx = good_dist_idx;
        else
            avg_idx = true(size(score{qid}{sid}));
        end
        avg_idx(min(avg_topn+1,end):end) = false;
        sel_topic_bows = [topic_bows{qid}{sid},db_bow(:,ranks{qid}{sid}(avg_idx))]; 
        if strcmp(exp_type,'avg_dist')
            multi_dists = comp_dist(ivf,sel_topic_bows,db_bow,dist_type,false);
            dists{qid}{sid} = mean(multi_dists,2);
        elseif ~isempty(strfind(exp_type,'avg_pooling'))
            if ~isempty(strfind(exp_type,'eq'))
                pos = strfind(exp_type,'_');
                wei = str2num(exp_type(pos(end)+1:end));
                sel_topic_bows = sel_topic_bows .* repmat([ones(1,size(topic_bows{qid}{sid},2))*wei/size(topic_bows{qid}{sid},2)*sum(avg_idx),ones(1,sum(avg_idx))],size(sel_topic_bows,1),1);  % weighting with wei,1,1,1...
            elseif ~isempty(strfind(exp_type,'power'))
                pos = strfind(exp_type,'_');
                wei = str2num(exp_type(pos(end)+1:end));
                sel_topic_bows = sel_topic_bows .* repmat([ones(1,size(topic_bows{qid}{sid},2))*wei/size(topic_bows{qid}{sid},2),0.5.^[1:sum(avg_idx)]],size(sel_topic_bows,1),1);  % weighting with wei,1/2,1/4,1/8...
            end
            avg_topic_bow = sparse(sum(sel_topic_bows,2));
            % normalize bow
            avg_topic_bow = norm_bow(avg_topic_bow,norm_type,size(sel_topic_bows,2));
            dists{qid}{sid} = comp_dist(ivf,avg_topic_bow,db_bow,dist_type,false);
        else
            fprintf('Error: unrecognized expansion type: %s', exp_type);
            exit;
        end
        if is_ransac_res
            % we always trust ransac results
            % that means we keep ransac results first

            % re-ranked scores are set to above 0
            min_dist = min(dists{qid}{sid});
            if min_dist<0
                dists{qid}{sid} = dists{qid}{sid} - min_dist;
            end
            % inject the minus scores of previous ransac results
            dists{qid}{sid}(ranks{qid}{sid}(good_dist_idx))...
                = score{qid}{sid}(good_dist_idx);
        end 
        [new_score{qid}{sid},new_ranks{qid}{sid}]=sort(dists{qid}{sid},1);
    end 
    total_time = total_time + toc;
    fprintf('\r%d/%d %.0fs',qid,query_num,total_time);
end
fprintf('\n%.0fs\n',toc);

fprintf('clean invfile...');
tic;
mxCleanInvFile(ivf);
fprintf('%.0fs\n',toc);
clear ivf
