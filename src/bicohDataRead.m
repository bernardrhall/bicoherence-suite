function [ timeSeries fs ] = bicohDataRead(framePath, frameFile, channelName, ...
                         startTime, duration, write_text, outpath);

%writes a time series text/data file from a frame file

%fullchName = strcat(channelName,ifo,':');

%addpath(framePath);
startTime = str2double(startTime);
duration = str2double(duration);
write_text = str2num(write_text);

frameFile = fullfile(framePath,frameFile)

[data times] = frgetvect(frameFile,channelName,startTime,duration,1);

if times(1) < 1000
	times = times + startTime;
end

timeSeries = [times data];

[x y] = size(data);

fs = round(x/duration);

if write_text
	nameString = strcat(outpath,'/',num2str(startTime),'+',channelName,'+',num2str(fs),'+',num2str(duration),'.txt');
	tfile = fopen(nameString,'w');
	for i = 1:x
		ln = char(strcat(num2str(timeSeries(i,1),16),{' '},num2str(timeSeries(i,2),16)));
		fprintf(tfile,'%s\r\n',ln);
	end
	fclose(tfile);
	%dlmwrite(nameString,timeSeries,' ','precision', 16);
end
