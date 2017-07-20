function [ timeSeries fs ] = cohDataRead(framePath, frameFile, channelName, ...
                         startTime, duration, write_text);

%wtites a time series text/data file from a frame file

%fullchName = strcat(channelName,ifo,':');

%addpath(framePath);

frameFile = fullfile(framePath,frameFile)

[data times] = frgetvect(frameFile,channelName,startTime,duration,1);

timeSeries = [data times];

[x y] = size(data);

fs = round(x/duration);

if write_text
	nameString = strcat(num2str(startTime),'+',channelName,'+',num2str(fs),'+',num2str(duration),'.txt');

	dlmwrite(nameString,timeSeries,' ');
end
