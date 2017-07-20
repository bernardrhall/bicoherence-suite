function [ clusters_f ] = cluster_find( fr, frq )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
j = 1;
len = length(fr);
clusters_f = cell(len,1);

%-------------------Make sure only one of mirrored bifrequencies are taken
%frSort = fr(:,1:2);
%[~,a,~]=unique(sort(frSort,2),'rows','stable');
%fr = fr(a,:);
%-------------------------

if len > 0
	fr = sortrows(fr,1);

%test = '';
for i = 1:length(fr(:,7))
    a = fr(i,7);
    b = fr(i,8);
    nxt = 0;
    if (i+1 <=length(fr(:,7)))
        if ((fr(i+1,7)==a+1)&&(fr(i+1,8)==b+1)) || ((fr(i+1,7)==a)&&(fr(i+1,8)==b+1)) ...
                || ((fr(i+1,7)==a+1)&&(fr(i+1,8)==b))
            %disp(strcat(test,',',num2str(a),',',num2str(b),',',num2str(fr(i+1,7)),',',num2str(fr(i+1,8))))
            %clusters_f{j} = strcat(clusters_f{j},',',num2str(a),',',num2str(b),',',num2str(fr(i+1,7)),',',num2str(fr(i+1,8)));
            clusters_f{j} = strcat(clusters_f{j},',',num2str(a),',',num2str(b));
            nxt = 1;
        else
            if nxt
                clusters_f{j} = strcat(clusters_f{j},',',num2str(fr(i+1,7)),',',num2str(fr(i+1,8)));
            end
            if j <= length(fr)
                j = j+1;
            end
        end
    end
end
clusters_f = clusters_f(~cellfun('isempty', clusters_f));
len = length(clusters_f);
for i = 1:length(clusters_f)
    br1 = strsplit(clusters_f{i,1},',');
    br1 = br1(~cellfun('isempty', br1));
    %clusters_f{i,1} = unique(br1);
    clusters_f{i,1} = br1;
    
end

col = cell(len,1);
clusters_f = [clusters_f col];

for i = 1:len
    for k = 1:length(clusters_f{i,1})
        clusters_f{i,2}(k) = frq(str2double(clusters_f{i,1}(k)));
    end
end
% -------------------old endpoint
else
	col = cell(len,1);
	clusters_f = [clusters_f col]; % just return empty cell in this case
end % end check fr length condition
end % fuction end
