function makePSD(opath,ipath,ts,fs,fontsize,lg,fg)

disp('Loading time series data...')
timeSeries = load(strcat(ipath,'/',ts));
disp('Done!')

timeSeries = timeSeries(:,2);

[pxx,w] = periodogram(timeSeries,[],[],str2double(fs));

if str2num(lg)
	figure
	h = semilogx(w,10*log10(pxx));
	% Set the axis limits and turn on the grid
	%axis([min(w) max(w) min(10*log10(pxx)) max(10*log10(pxx))])
	grid on
	grid minor
else
	figure
	h = plot(w,10*log10(pxx));
	grid on
	grid minor
end
C = strsplit(ts,'.');
C = C(1:(end-1));
C = strjoin(C,'.');
C1 = strrep(C,'.','_d_');

title_string = strcat('Periodogram for',{' '},strrep(C,'_','\_'));
title_string = strrep(title_string,'+',' ');
title_string = textwrap(title_string,40);

%title_string = char(strcat('Periodogram for',{' '},C));
%grid on
%grid minor
ax = gca;
%set(ax,'XMinorTick','on')  % sets Minor X Ticks to display

ax.XTick = [10^-2 10^0 10^1 10^2 10^3 10^4];
title(title_string)
set(ax,'YDir','normal');
set(ax,'FontSize',str2double(fontsize));
xlabel('frequency (Hz)'), ylabel('dB')

saveString = strcat(opath,'/','PSD_',C1);
disp(char(strcat('Saving',{' '},saveString,'.png')))
%plot2svg(strcat(GW_channel_n{1},'.svg'), plt, 'png')
%saveas(plt,strcat(opath,'/',GWCN),'epsc')
saveas(h,saveString,'png')

if str2num(fg)
	disp(char(strcat('Saving',{' '},saveString,'.fig')))
	saveas(h,saveString,'fig')
end
