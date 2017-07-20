function [cl_indx_M, clusters_M] = cl_merge( cl_I, freq )

    j = 1;

    cl_indx_M = {};
    clusters_M = {};

    for i = 1:length(cl_I)
        if i < length(cl_I)
            a = cl_I{i,1}(length(cl_I{i,1}(:,1)),1);
            b = cl_I{i+1,1}(1,1);

            cc1 = min(cl_I{i,1}(:,2));
            cc2 = max(cl_I{i,1}(:,2));
            dd1 = min(cl_I{i+1,1}(:,2));
            dd2 = max(cl_I{i+1,1}(:,2));
            
            c2_flag = 0;
            
            if (dd1 >= (cc1 - 1) && dd1 <= (cc2 + 1)) || ... 
                (dd2 >= (cc1 - 1) && dd2 <= (cc2 + 1))
                c2_flag = 1;
            end
            
            if (b == a + 1) && c2_flag
                cl_indx_M{j} = [cl_I{i,1};cl_I{i+1,1}];
                clusters_M{j} = [freq(cl_I{i,1});freq(cl_I{i+1,1})];
                j = j + 1;
            end
        end
    end

    cl_indx_M = cl_indx_M';
    clusters_M = clusters_M';

end