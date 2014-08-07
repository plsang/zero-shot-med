function save_qinfo_for_c(tree,weight,norm,qinfo_pathname)
% should be consistent with load_qinfo in the online_retrieval_func.c
global M;
global K;
global num_centers;
method_type = 0; %VL_IKM_LLOYD    VL_IKM_ELKAN
dist_metric = 0; % VL_IKM_L2       VL_IKM_CHI
[M,K] = size(tree.centers);
D = tree.depth;
fqinfo = fopen(qinfo_pathname,'w');
fwrite(fqinfo,int32(str2double(norm(end))),'int');
fwrite(fqinfo,length(weight),'int');
fwrite(fqinfo,weight,'float');
fwrite(fqinfo,method_type,'int');
fwrite(fqinfo,dist_metric,'int');
fwrite(fqinfo,M,'int');
fwrite(fqinfo,K,'int');
fwrite(fqinfo,D,'int');

num_centers=0;

flog = fopen('/home/caizhizhu/cflog.txt','w');
save_tree_centers(fqinfo,flog,tree,D);
fclose(flog);
assert(num_centers == K*(K^D-1)/(K-1));
fclose(fqinfo);
end

function save_tree_centers(fqinfo,flog,tree,height)
global M;
global K;
global num_centers;
centers = tree.centers;
KK = size(tree.centers,2);
if KK<K
    centers = [centers,zeros(M,K-KK)];   %padded with zeros
    if height > 1
        for k = KK+1:K
            tree.sub(k).centers = [];
        end
    end
end
fwrite(fqinfo, centers,'int');
fprintf(flog,'off: %d hei: %d 1st:%d 2rd:%d end:%d\n',num_centers/K, height,centers(1,1),centers(2,1),centers(end,end));
num_centers = num_centers+K;
if height > 1
    for k=1:K
        save_tree_centers(fqinfo,flog,tree.sub(k),height-1);
    end
end
end