function makeSpectrogram(opath,ipath,ts,fs,fontsize)

disp('Loading time series data...')
timeSeries = load(strcat(ipath,'/',ts));
disp('Done!')

timeSeries = timeSeries(:,2);

[~,f,t,s] = spectrogram(timeSeries,[],2000,2048,str2double(fs),'yaxis');
%[~,~,~,s] = spectrogram(timeSeries,'yaxis');

%[pxx,w] = periodogram(timeSeries,[],[],str2double(fs));

h = surf(t,f,10*log10(s))
view(0, 90)
axis tight
set(h, 'LineStyle','none')
%h = imagesc(t,f,10*log10(s))
C = strsplit(ts,'.');
C = C{1};

colormap jet
colorbar
title_string = char(strcat('Spectrogram for',{' '},C));
title(title_string)
set(gca,'YScale','log');
set(gca,'YDir','normal');
set(gca,'FontSize',str2double(fontsize));
%xlabel('frequency (Hz)'), ylabel('dB')

saveString = strcat(opath,'/','specGm_',C);
disp(char(strcat('Saving',{' '},saveString,'.png')))
%plot2svg(strcat(GW_channel_n{1},'.svg'), plt, 'png')
%saveas(plt,strcat(opath,'/',GWCN),'epsc')
saveas(h,saveString,'png')

end
