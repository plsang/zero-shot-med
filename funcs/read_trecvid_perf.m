function perf = read_trecvid_perf(trecvid_perf_filename)
perf_fid = fopen(trecvid_perf_filename,'r');
perf_cell = textscan(perf_fid,'%s %s %f');
count = 0;
perf=[];
for n = 1:length(perf_cell{1})
    if strcmp(perf_cell{1}(n),'infAP')
        perf = [perf, perf_cell{3}(n)];
    end
end
fclose(perf_fid);
end
