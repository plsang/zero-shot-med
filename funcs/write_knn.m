function write_knn(query_filenames, db_lut, score, ranks, db_frame_dir, knn_dir, topk, num_shown_frames,clobber)
if ~exist('clobber','var')
    clobber = false;
end
assert(length(query_filenames) == length(score));
if ~exist(knn_dir,'dir')
    mkdir(knn_dir);
end
pos = strfind(db_frame_dir,'ins20');
DB = db_frame_dir(pos:pos+length('ins2012')-1);
query_num = length(query_filenames);
for qid = 1:query_num
    fprintf('\r%d(1-%d)',qid,query_num);
    subset_num = length(query_filenames{qid});
    for sid = 1:subset_num
        pathstr=fileparts(query_filenames{qid}{sid}{1});
        [~, topic_name]=fileparts(pathstr);
        
        knn_filename = fullfile(knn_dir,sprintf('%s.%d.txt',topic_name,sid));
        assert(size(ranks{qid}{sid},2) == 1);
        if ~exist(knn_filename,'file') || clobber
            % print knn list
            fid = fopen(knn_filename,'w');
            img_num = length(query_filenames{qid}{sid});
            for iid = 1:img_num
                pos = strfind(query_filenames{qid}{sid}{iid},'/');
                %topic_name/query_name.jpg
                fprintf(fid, '%s ',query_filenames{qid}{sid}{iid}(pos(end-1)+1:end));
            end
            fprintf(fid,'\n');
            for i=1:topk
                clip_name = db_lut{ranks{qid}{sid}(i)};
                fprintf(fid, '%s(dist_%.4f): ',clip_name,score{qid}{sid}(i));
                frame_list_file = fullfile('/tmpfs', DB, 'frame_list',[clip_name,'.txt']);
                if ~exist(frame_list_file,'file')  
                    % slow but in local frame dir---- for other members, later can be deleted
                    clip_dir = fullfile(db_frame_dir,clip_name);
                    frame_list_file = fullfile(clip_dir,'frames.txt');
                    if ~exist(frame_list_file,'file')
                        if ~isempty(strfind(clip_dir,'ins2012')) %should be remove later
                            frames = dir(fullfile(clip_dir,'*.jpg'));
                        else
                            frames = dir(fullfile(clip_dir,'*.png'));
                        end
                        frames = {frames(:).name};
                        num_frames = length(frames);
                        if num_frames == 0
                            fprintf(fid, '\n');
                            continue
                        end
                        
                        list_fid=fopen(frame_list_file,'w');
                        for j=1:num_frames
                            fprintf(list_fid,[frames{j} '\n']);
                        end
                        fclose(list_fid);
                    end
                end
                if exist(frame_list_file,'file')~=0  
                    list_fid=fopen(frame_list_file,'r');
                    frames = textscan(list_fid,'%s');
                    frames = frames{1};
                    fclose(list_fid);
                    num_frames = length(frames);
                    if num_frames == 0
                        fprintf(fid, '\n');
                        continue
                    end
                end
                rep_frame_ids = linspace(1,num_frames,min(num_frames,num_shown_frames));
                rep_frame_ids = floor(rep_frame_ids);
                for j=1:length(rep_frame_ids)
                    fprintf(fid, '%s ',frames{rep_frame_ids(j)}(1:end-4));
                end
                fprintf(fid, '\n');
            end
            fclose(fid);
        end
    end
end
