function [data, fs, time] = getTs(framePath,frName,debugLevel,allowRedundantFlag,channelName,outpath,write_text)  

%framePath = '1126254778.499023_1126254790.499023.lcf';
%frameType = 'H1_R';
%debugLevel = 1;
%allowRedundantFlag = 0;
%startTime = 1126254778.499023;
%stopTime = 1126254790.499023;
%channelName = 'H1:CAL-DELTAL_EXTERNAL_DQ';

debugLevel = str2num(debugLevel);
allowRedundantFlag =str2num(allowRedundantFlag);

frameCache = loadframecache(framePath)

frameType = frameCache.frameTypes;
frameType = frameType{1};

div = strsplit(frName,'.');
div = div(1:end-1);
if length(div) > 1
        %jn = '.';
        div = strjoin(div,'.');
else
        div = div{1};
end
div = strsplit(div,'_');
startTime = str2double(div{1});
stopTime = str2double(div{2});
duration = round(stopTime - startTime);

[data, fs, time] = ...
           readframedata(frameCache, channelName, frameType, ...
                         startTime, stopTime, allowRedundantFlag, ...
                         debugLevel);

timeSeries = [time' data'];

x = length(data);

if write_text
	if (fs > 0) && (x > 0)
		disp('Writing time-series text file...')
	        nameString = strcat(outpath,'/',num2str(startTime),'+',channelName,'+',num2str(fs),'+',num2str(duration),'.txt');
	        tfile = fopen(nameString,'w');
        	for i = 1:x
                	ln = char(strcat(num2str(timeSeries(i,1),16),{' '},num2str(timeSeries(i,2),16)));
	                fprintf(tfile,'%s\r\n',ln);
        	end
	        fclose(tfile);
		disp('Done!')
	        %dlmwrite(nameString,timeSeries,' ','precision', 16);
	else
		disp('No data!...noting')
		missing = fopen(strcat(outpath,'/',num2str(startTime),'_missing'),'w');
                fclose(missing);
	end
end
