function makeGaussian(seconds, GPS, fs, amplitude, startCount, count, directory)

rng('shuffle')
%---------------------------
%seconds = 8;
%GPS = 1134628650;
%fs = 4096;
%%amplitude = 0.01;
%amplitude = 10^(-19);
%directory = 'gaussian_1000_source';
%count = 1000;
%------------------------------
count = str2num(count);
startCount = str2num(startCount);
amplitude = str2double(amplitude);
fs = str2double(fs);
GPS = str2double(GPS);
seconds = str2double(seconds);
%------------------------------
N = seconds * fs;
col1 = ((1:N)/fs) + GPS;

for i = startCount:count
    disp(char(strcat('Creating time series',{' '},num2str(i),{' '},'of',{' '},num2str(count))))
    GtimeSeries = amplitude * randn(1,N);
    
    GtimeSeries = [col1;GtimeSeries];
    GtimeSeries = GtimeSeries';
    
    disp('Writing to file...')
    outfile = fopen(strcat(directory,'/GTS',num2str(i),'.txt'),'w');
    for j = 1:N
        fprintf(outfile,'%s\r\n',char(strcat(num2str(GtimeSeries(j,1),16),{' '},num2str(GtimeSeries(j,2),16))));
    end
    fclose(outfile);
end

end
