function save_trec_perf_list( query_filenames,score,ranks,lut,run_name,trecvid_res_dir,eval_topN,ex_list )
for qid=1:length(query_filenames)
    subset_num = length(query_filenames{qid});
    for sid = 1:subset_num
        pathstr=fileparts(query_filenames{qid}{sid}{1});
        [~, topic_name]=fileparts(pathstr);
        trec_list_file = fullfile(trecvid_res_dir,sprintf('%s.%d.list',topic_name,sid));
        if exist(trec_list_file,'file')
            continue;
        end
        trec_fid = fopen(trec_list_file,'w');
        score{qid}{sid} = max(score{qid}{sid})-score{qid}{sid};
        assert(size(ranks{qid}{sid},2) == 1);
        for j = 1:eval_topN
            ret = lut{ranks{qid}{sid}(j)};
            if exist('ex_list', 'var') && ismember(ret,ex_list)
                continue;
            end
            fprintf(trec_fid, '%s 0 %s %d %G %s\n', topic_name,ret,j,...
                score{qid}{sid}(j),run_name);
        end
        fclose(trec_fid);
    end
end
end

