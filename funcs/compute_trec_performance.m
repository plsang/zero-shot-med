function perf = compute_trec_performance(query_filenames,gt_filename, trecvid_res_dir, eval_topN, gt)
query_num = length(query_filenames);
perf = cell(1,query_num);
for qid=1:query_num
    subset_num = length(query_filenames{qid});
    perf{qid} = zeros(1,subset_num);
    for sid = 1:subset_num
        pathstr=fileparts(query_filenames{qid}{sid}{1});
        [~, topic_name]=fileparts(pathstr);
        trec_list_file = fullfile(trecvid_res_dir,sprintf('%s.%d.list',topic_name,sid));
        assert(exist(trec_list_file,'file')~=0);
        
        trecvid_perf_filename = strrep(trec_list_file,'list','perf');
		trecvid_perf_filename
        if ~exist(trecvid_perf_filename,'file')
            [status, performance]=unix(sprintf('/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/code/funcs/trec_eval -q -a -c %s %s %d', gt_filename, trec_list_file, eval_topN));
            assert(status == 0);
            perf_fid = fopen(trecvid_perf_filename,'wt');
            fwrite(perf_fid,performance);
            fclose(perf_fid);
        end
        map = read_trecvid_perf(trecvid_perf_filename);
        perf{qid}(sid) = map(1);
        if length(gt.(['id_' topic_name]))>eval_topN
            perf{qid}(sid) = perf{qid}(sid)*length(gt.(['id_' topic_name]))/eval_topN;
        end
    end
end
end
