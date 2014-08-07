function [ score, rank ] = imageprop(query_details, db_frame_dir, ...
    rank_list_dir,  video_rank, db_lut)

feature_comb = fieldnames(query_details);
rep_video_idx=query_details.(feature_comb{1}).rep_video_idx;
query_names=query_details.(feature_comb{1}).query_names;

[rerank_topk,num_query] = size(rep_video_idx);

old_rank_list_dir = fullfile(rank_list_dir,'old');
save_old_rank_list(old_rank_list_dir, query_names, db_frame_dir, ...
     rep_video_idx, video_rank, db_lut, rerank_topk);

new_rank_list_dir = fullfile(rank_list_dir,'new');
new_rank_detail_dir = fullfile(rank_list_dir,'match');
image_propogate(old_rank_list_dir, new_rank_list_dir, new_rank_detail_dir);

verbose = true;
score = load_new_rank_list(new_rank_list_dir, query_names, rerank_topk, verbose);
end

function save_old_rank_list(rank_list_dir, query_names, db_frame_dir, ...
     rep_video_idx, video_rank, db_lut, rerank_topk)
if ~exist(rank_list_dir,'dir')
    mkdir(rank_list_dir);
end
num_query = length(query_names);
for i=1:num_query
    % get query topic name
    topic_path = fileparts(query_names{i}{1});
    [~,topic_name,~] = fileparts(topic_path);
    ftopic_rank_list = fopen(fullfile(rank_list_dir,[topic_name '.txt']),'w');
    num_query_img = length(query_names{i});
    fprintf(ftopic_rank_list,'%d %d\n',num_query_img,rerank_topk);
    for j=1:num_query_img
        fprintf(ftopic_rank_list,'%s\n',query_names{i}{j});
    end
    for j=1:rerank_topk
        fprintf(ftopic_rank_list,'%s\n',fullfile(db_frame_dir,db_lut{video_rank(i,j)},...
            sprintf('%s_%04d.jpg',db_lut{video_rank(i,j)},rep_video_idx(j,i))));
    end
    fclose(ftopic_rank_list);
end
end

function image_propogate(old_rank_list_dir, new_rank_list_dir, new_rank_detail_dir)
if ~exist(new_rank_list_dir,'dir')
    mkdir(new_rank_list_dir);
end
if ~exist(new_rank_detail_dir,'dir')
    mkdir(new_rank_detail_dir);
end
old_rank_list = dir(fullfile(old_rank_list_dir,'*.txt'));
old_rank_list = {old_rank_list(:).name};
cmd = './ImageExploration/bin/learnMultiViewMosaic';
for i=1:length(old_rank_list)
    unix(sprintf('%s %s >/dev/null 2>&1 &', cmd,fullfile(old_rank_list_dir,old_rank_list{i})));
end
end

function score = load_new_rank_list(new_rank_list_dir, query_names, rerank_topk, verbose)
num_query = length(query_names);
score = zeros(num_query, rerank_topk);
new_rank_list = dir(fullfile(new_rank_list_dir,'*.txt'));
new_rank_list = {new_rank_list(:).name};
assert(num_query == length(new_rank_list));
for i=1:num_query
    ftopic_rank_list = fopen(fullfile(new_rank_list_dir,new_rank_list{i}));
    num_query_img = length(query_names{i});
    components = textscan(ftopic_rank_list,repmat('%d ',1,num_query_img));
    components = cell2mat(components);
    total_rows = size(components,1);
    assert(rerank_topk+num_query_img == total_rows);    
    if verbose
        for j=1:num_query_img
            fprintf('ranking of query itself: %s \n',query_names{i}{j});
            ind = true(1,total_rows);
            ind(j) = false;
            [~,order] = sort(components(ind,j),'descend');
            for n=1:min(length(order),10)
                fprintf('%d ',order(n));
            end
            fprintf('\n');
        end
    end
    score(i,:) = sum(components(num_query_img+1:end,:),2);
    
    fclose(ftopic_rank_list);
end
end