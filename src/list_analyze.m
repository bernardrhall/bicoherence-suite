function list_analyze( path, tol )

dirString = strcat(path,'/','*.mat');
directory = dir(dirString);
nfiles = length(directory);
event = cell(nfiles,1);
sEvent = cell(nfiles,1);
for i = 1 : nfiles
    clearvars -except directory i nfiles event path tol sEvent
    list = directory(i).name;
    disp('Loading variables...')
    q = matfile(strcat(path,'/',list));
    %load(strcat(path,'/',list));
    disp(strcat('Variables loaded--analyzing...',num2str(i),'...of...',num2str(nfiles)))
    [cl_indx_Merge,clusters_Merge,cl,clusters,cl_indx] = cluster_analyze(q.final_results,q.freqBsp_F,str2num(tol));
    [a,b] = size(clusters);
    [c,d] = size(clusters_Merge);
    if (a > 0 && b > 0) && (c == 0 || d == 0)
        disp('Minor cluster...recording...')
        disp(list)
	event{i} = list;
    elseif (c > 0 && d > 0)
        disp('Major cluster...recording...')
        disp(list)
	sEvent{i} = list;
	disp('Saving cluster variables...')
	nm = strsplit(list,'.');
	nm = nm(1:(end-1));
	nm = strjoin(nm,'.');
    	save(strcat(path,'/',nm,'_clusters.mat'), 'cl_indx_Merge', 'clusters_Merge', 'cl', 'clusters', 'cl_indx', '-v7.3');
    	disp('Done!')
    end
end

event = event(~cellfun('isempty',event));
sEvent = sEvent(~cellfun('isempty',sEvent));

file = fopen(strcat(path,'/','event.dat'),'w');
for j = 1:length(event)
	fprintf(file,'%s\n',event{j});
end
fclose(file);

file1 = fopen(strcat(path,'/','sEvent.dat'),'w');
for j = 1:length(sEvent)
        fprintf(file1,'%s\n',sEvent{j});
end
fclose(file1);

save(strcat(path,'/','event.mat'), 'event', '-v7.3');
save(strcat(path,'/','sEvent.mat'), 'sEvent', '-v7.3');

end
