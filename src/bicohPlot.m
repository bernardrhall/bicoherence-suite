function bicohPlot(ipath,opath,mFile,mName,fontsize,time_window,fs,findCL,tol,mleg,overload,fpng)

%mName = genvarname(mName);
disp('Loading...')
load(strcat(ipath,'/',mFile))
disp('done...')
fs = str2num(fs);
%disp(mName)
if ischar(mName)  % check to see if mName is interpreted as string; if so, evaluate so it is matrix variable name
	mName = eval(mName);
end
%mName = genvarname(mName);
findCL = str2num(findCL);
GW_channel = mFile;
GWCN = strsplit(GW_channel,'.');
GWCN = GWCN(1:(end-1));
GWCN = strjoin(GWCN,'.');
flag = 0;
if findCL
	disp('Retrieving clusters...')
	clfile = matfile(strcat(ipath,'/',GWCN,'_clusters.mat'));
	clMerge = clfile.clusters_Merge;
        %[ clindexMerge, clMerge, cl, cluster, clusterIndex ] = cluster_analyze( final_results, freqBsp_F, str2num(tol));
        [x, ~] = size(clMerge);
	GWCN1 = strrep(GWCN,'.','_d_'); % need to replace '.' with '_d_' (or other valid character) so that the .fig/.png filename will work
	if x <= overload
		%sct = zeros(x,1);
		figure; colormap lines
		cmap = colormap;
		[lmap, ~] = size(cmap);
		%smap = randi([0,255],x,3)
		smap = rand([x 3]);
		%hld = 1;
		hold all
	        for i = 1:x
        	        %fileID = fopen(strcat('clusters_',GWCN{1},'_',num2str(i),'.dat'),'w');
			disp(char(strcat('Writing',{' '},strcat(opath,'/','clusters_',GWCN,'_',num2str(i),'.dat'))))
                	dlmwrite(strcat(opath,'/','clusters_',GWCN,'_',num2str(i),'.dat'),clMerge{i,1},' ');
			disp('Plotting clusters...')
			if i <= lmap
				scatter(clMerge{i,1}(:,1),clMerge{i,1}(:,2),[],cmap(i,:));
			else
				scatter(clMerge{i,1}(:,1),clMerge{i,1}(:,2),[],smap(i,:))
			end
			lgd{i} = num2str(i);
        	end
        	set(gca,'YDir','normal');
	        set(gca,'FontSize',str2double(fontsize));
        	xlabel('f1 (Hz)'), ylabel('f2 (Hz)')
		title_string_A = strcat('Clusters in',{' '},strrep(GWCN,'_','\_'));
		title_string_A = strrep(title_string_A,'+',' ');
		title_string_A = textwrap(title_string_A,40);
		title(title_string_A)
		%title(char(strcat('Clusters in',{' '},strrep(GWCN,'_','\_'))))
		grid on
		if str2num(mleg)
			legend(lgd);
		end
		hold off
		saveas(gca,strcat(opath,'/',GWCN1,'_','clusters'),'png')
	else
		disp('Cluster overload!...Skipping...')
		dlmwrite(strcat(opath,'/','cl_overload_',GWCN,'_','.dat'),x);
		flag = 1;
	end
end
%  plot figure
disp('Plotting matrix...')
figure; colormap jet
plt = imagesc(freqBsp_F,freqBsp_F,abs(mName));
colorbar
GW_channel_n = strrep(GW_channel,'_','\_');
%GW_channel_n = regexprep(GW_channel_n,'\d','');
GW_channel_n = strsplit(GW_channel_n,'.');
title_string = strcat('(Auto) Bicoherence of&',GW_channel_n{1},' (@ f_s =&',num2str(fs),'),\newline GPS:&',num2str(GPS_central),', Dur. =&', ...
time_window,' sec.');
title_string = strrep(title_string,'&',' ');
title(title_string)
set(gca,'YDir','normal');
set(gca,'FontSize',str2double(fontsize));
xlabel('f1 (Hz)'), ylabel('f2 (Hz)')
disp(char(strcat('Saving',{' '},strcat(opath,'/',GWCN),'.fig')))
%plot2svg(strcat(GW_channel_n{1},'.svg'), plt, 'png')
%saveas(plt,strcat(opath,'/',GWCN),'epsc')
saveas(plt,strcat(opath,'/',GWCN1),'fig')
if flag
	saveas(plt,strcat(opath,'/',GWCN1,'_','clusters.png'))
end
if fpng && ~flag
	saveas(plt,strcat(opath,'/imsc_',GWCN1),'png')
end
%saveas(plt,strcat(opath,'/',GWCN),'m')
%savefig(plt,strcat(opath,'/',GWCN))
%close(plt)
