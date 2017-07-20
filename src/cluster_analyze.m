function [ ciM, clMerge, cl, clust, cli ] = cluster_analyze( f_r, f_F , tol)

cl = cluster_find(f_r,f_F);

[r,~] = size(cl);
%tol = 4;  % cluster size cutoff-- 4 usually

clust = cell(r,1);
cli = cell(length(r),1);

for i = 1:r
    if length(cl{i,1}) > tol
            %disp(cl{i,1})
            mtrx = zeros(length(cl{i,1})/2,2);
            mtrx1 = zeros(length(cl{i,1})/2,2);
            n = 1;
            p = 1;
            for m = 1:length(cl{i,1})
                if rem(m,2)
                    mtrx(n,p) = str2double(cl{i,1}(m));
                    mtrx1(n,p) = cl{i,2}(m);
                    p = p+1;
                else
                    mtrx(n,p) = str2double(cl{i,1}(m));
                    mtrx1(n,p) = cl{i,2}(m);
                    n = n+1;
                    p = 1;
                end
            end
            % assign, and elimate any possible duplicates
            cli{i,1} = unique(mtrx,'rows');
            clust{i,1} = unique(mtrx1,'rows');
            clearvars mtrx mtrx1
    end
end

%clearvars cl

clust = clust(~cellfun('isempty', clust));
cli = cli(~cellfun('isempty', cli));

%clearvars a b i j m n p r
[ciM,clMerge] = cl_merge(cli,f_F);