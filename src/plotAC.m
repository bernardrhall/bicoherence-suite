function plotAC(matrix,fs,fontsize,opath,number,gps)

mat = matfile(matrix);

figure; 
%h = plot(fs*mat.freq,abs(mat.acoh));
fs = str2num(fs);
%number = str2num(number);

h = semilogx(fs*mat.freq,abs(mat.acoh));
%title_string = char(strcat('Periodogram for',{' '},C));
title_string = 'Autocoherence';
title(title_string)
set(gca,'YDir','normal');
set(gca,'FontSize',str2double(fontsize));
xlabel('frequency (Hz)'), ylabel('')
if strcmp(opath,'.')
	saveString = strcat(gps,'_autocoherence_',number);
else
	saveString = strcat(opath,'/',gps,'_autocoherence_',number);
end

disp(char(strcat('Saving',{' '},saveString,'.png')))
%plot2svg(strcat(GW_channel_n{1},'.svg'), plt, 'png')
%saveas(plt,strcat(opath,'/',GWCN),'epsc')
saveas(h,saveString,'png')
