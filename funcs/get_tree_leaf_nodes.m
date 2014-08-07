function leaf_nodes = get_tree_leaf_nodes(tree)
D = tree.depth;
[M,K] = size(tree.centers);

if D <= 1
    leaf_nodes = tree.centers;
elseif D==2
    n = 1;
    leaf_nodes=zeros(M,K^D);
    for a = 1:K
        len=size(tree.sub(a).centers,2);
        leaf_nodes(:,n:n+len-1) = tree.sub(a).centers;
        n = n+K;
    end
elseif D==3
    n = 1;
    leaf_nodes=zeros(M,K^D);
    for a = 1:K
        for b = 1:K
            len=size(tree.sub(a).sub(b).centers,2);
            leaf_nodes(:,n:n+len-1) = tree.sub(a).sub(b).centers;
            n = n+K;
        end
    end
elseif D==6    
    n = 1;
    leaf_nodes=zeros(M,K^D);
    for a = 1:K
        for b = 1:K
            for c = 1:K
                for d = 1:K
                    for e = 1:K
                        len=size(tree.sub(a).sub(b).sub(c).sub(d).sub(e).centers,2);
                        leaf_nodes(:,n:n+len-1) = tree.sub(a).sub(b).sub(c).sub(d).sub(e).centers;
                        n = n+K;
                    end
                end
            end
        end
    end
else
    fprintf('unsupported tree depth (1,2,3,6 only) %d',D);
end
end

