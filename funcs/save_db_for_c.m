function save_db_for_c(database,lut,db_pathname,lut_pathname)
num_item = size(database,2);
assert(num_item==length(lut));
fdb = fopen(db_pathname,'w');
flut = fopen(lut_pathname,'w');
for i=1:num_item
    fprintf('\r%d/%d',i,num_item);
    fprintf(flut,'%d:%s\n',i,lut(i,:));
    fwrite(fdb,i,'int32');
    nid = find(database(:,i));
    num_nz = length(nid);
    fwrite(fdb,num_nz,'int32');
    for j=1:num_nz
        fwrite(fdb,nid(j),'int32');
        fwrite(fdb,database(nid(j),i),'float');        
    end
end
fclose(fdb);
fclose(flut);
