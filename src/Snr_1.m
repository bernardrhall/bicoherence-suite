function [ S , sDb, noise, signal] = Snr_1( ts,fs,mn,verbose )
% This function attempts to estimate the SNR of a signal in a method closer
% to to that described by Taekhyun Kim.  Currently, the "variance" of the
% noise is found as the 2 times the standard deviation of the signal minus
% the peaks identified in the power spectrum.  The "signal" is found as the
% square root of the sum of the peaks in the power spectrum.  Then the SNR
% is found as the ratio of these signal/noise_variance.  This seems to work
% best when used in determining a bicoherence threshold.  This gives values
% for the critical angle that are sufficiently small when noise is minimal
% or absent, and increases the critical angle as expected when the noise
% level rises.  A final formula may need to be determined with more
% experimentation.

ts_trim = 0.25;
le = 5;% 3
lc = 4; % 2
uc = 4; % 2
ue = 5; % 3

% le = 20;
% lc = 18;
% uc = 18;
% ue = 20;

%test = var(ts)

ts = ts - mean(ts);
%[Pxx, F] = periodogram(ts,kaiser(numel(ts),2),numel(ts),fs);
[Pxx, F] = periodogram(ts,[],numel(ts),fs);
test = var(Pxx);
%Y = peak2rms(Pxx);
%figure; plot(F,10*log10(Pxx));
multiplier = max(Pxx);
Pxx = Pxx/multiplier;
%std(Pxx)
%hn = 0;
ln = 0;
prom = max(Pxx)/10;
if test < 0.0001
    dvs = 10;
    ln = 1;
elseif test < 4
    dvs = 10;  % default 4?
elseif test < 35
    dvs = 3;
elseif test < 85
    dvs = 1.5;
else
    dvs = 1;
end

%dvs = 10;
normal = 1;
%prom = max(Pxx)/20;
%prom = max(Pxx)/1.5;

% if mean(Pxx) < 1
%     prom = max(Pxx)/20;
% else
%     prom = max(Pxx)/1.5;
% end
%findpeaks(Pxx,fs,'MinPeakProminence',0.15)
%peaks = findpeaks(Pxx,fs,'MinPeakProminence',0.15);
%peaks = findpeaks(Pxx,fs,'MinPeakProminence',prom)
if ln
    %[peaks,indx] = findpeaks(log(Pxx),'MinPeakProminence',10.55);% *********
    %[peaks,indx] = findpeaks(log(Pxx),'MinPeakProminence',12.55);
    [indx,peaks] = peakfinder(log(Pxx),(max(Pxx)-min(Pxx))/dvs,10.55,1,1);
    peaks = exp(peaks);
else
    [indx,peaks] = peakfinder(Pxx,(max(Pxx)-min(Pxx))/dvs,prom,1,1);
end

%indx
%numel(peaks)
%Freqs = F((indx+1));
Freqs = F((indx));

ts1 = ts;

if (length(ts)/fs) < 10
    [xx,~] = size(ts1);
    if xx > 1
        ts1 = ts1';
    end
    %size(hann(length(ts1),'symmetric')')
    ts1 = ts1 .* hann(length(ts1),'symmetric')';
    ts_trim = 0.5;
    %ts1 = horzcat(ts1,fliplr(ts1));
    %ts1 = horzcat(ts1,zeros(10*length(ts1),1)');
end
% for m = 1:15
%     ts1 = horzcat(ts1,ts);
% end
% le1 = le;
% lc1 = lc;
% uc1 = uc;
% ue1 = ue;

goflag = 1;

for k = 1:numel(Freqs)
    disp('Checking parameters...')
    if ((Freqs(k) - le)) < 0 && ((Freqs(k) - lc) > 0)
        %le = Freqs(k)/4
        goflag = 0;
        disp('Using highpass (lstop okay)...')
        disp(strcat('Filtering frequency...',num2str(k),'...of...',num2str(numel(Freqs))))
        highpass_filter(Freqs(k),Freqs(k)+ue,fs);
    elseif (Freqs(k) - lc) < 0
        %lc = Freqs(k)/2
        goflag = 0;
        disp('Using highpass (lstop small)...')
        disp(strcat('Filtering frequency...',num2str(k),'...of...',num2str(numel(Freqs))))
        highpass_filter(Freqs(k),Freqs(k)+ue,fs);
    end
    if (Freqs(k) + uc) > (fs/2)
        %uc = (fs/2 - Freqs(k))/2;
        goflag = 0;
    end
    if (Freqs(k) - ue) > (fs/2)
        %ue = fs/2 - Freqs(k);
        goflag = 0;
    end
    
    if goflag
        disp(strcat('Filtering frequency...',num2str(k),'...of...',num2str(numel(Freqs))))
        ts1 = filter(notch_filter(Freqs(k)-le,Freqs(k)-lc,Freqs(k)+uc,Freqs(k)+ue,fs),ts1);
    end
    
    goflag = 1;
    % notch is slightly wider to properly eliminate signal frequencies
end

tindx = floor(length(ts1)*ts_trim); % trim beginning of time series where filter is not yet engaged
ts1 = ts1(tindx:length(ts1));

%figure; periodogram(ts1,[],numel(ts1),fs)  % use for filter test plot

%figure; plot(ts1)
%figure; plot(ts)
%peaks = findpeaks(Pxx,fs,'Threshold',0.0001)
%sig = mean(peaks);
%pwr_sig_sum = (sum(peaks)/2) * multiplier
if isempty(peaks)
    normal = 0;
end

if normal
    %sig = var(sqrt(peaks))
    if mn
        sig = mean(sqrt(peaks));
        signal = mean(sqrt(peaks*multiplier));
    else
        sig = (sum(peaks)/2);
        signal = sig * multiplier;
    end
    %sig = mean(sqrt(peaks));
    
    %signal = var(sqrt(peaks*multiplier))
    %cutoff = min(peaks);
else
    sig = max(Pxx)/2;
    signal = sig * multiplier;
    %cutoff = sig;
end
%sig * multiplier

% for i=1:numel(Pxx)
%     if Pxx(i) >= cutoff
%         Pxx(i) = 0;
%     end
% end

%mean(Pxx)
%if normal
    %varNoise = var(sqrt(Pxx))
    %varNoise = max(Pxx)

%varNoise = std(Pxx) * 10;% Pxx is already squared, so std(Pxx) ~ var(sqrt(Pxx))?
%varNoise = std(Pxx);
%sqrt(peak2rms(Pxx))
%varNoise = sqrt(evar(Pxx));
%varNoise = evar(ts1)/multiplier;
varNoise = var(ts1/sqrt(multiplier));
%varNoise = evar(ts/sqrt(multiplier));
%varNoise1 = evar(ts);
varNoise1 = var(ts1);

% noise_std = std(Pxx*multiplier);
% noise_rms = rms(Pxx*multiplier);

noise_std = std(ts1);
noise_rms = rms(ts1);

if verbose
    disp(strcat('Number of peaks found:...',num2str(numel(peaks))))
    disp(strcat('Noise amplititude RMS:...',num2str(noise_rms)))
    disp(strcat('(Normalized) Signal power:...',num2str(sig)))
    disp(strcat('(Normalized) Noise variance:...',num2str(varNoise)))
    disp(strcat('Noise variance:...',num2str(varNoise1)))
    disp(strcat('Noise amplititude std:...',num2str(noise_std)))
    disp(strcat('Signal Power:...',num2str(signal)))
end
%noise_mean = mean(Pxx*multiplier)
%multiplier
noise = noise_std;
    %noise = var(sqrt(Pxx*multiplier))
% else
%     varNoise = mean(sqrt(Pxx))
%     noise = varNoise * sqrt(multiplier)
%end
%noise = var(Pxx)
%S = sig/noise;
%S = sig/varNoise;  % ~ mean -> ~sqrt(SNR)

S = sqrt(signal)/(2*noise);  % this formula seems to work best
%S = signal/varNoise1;
sDb = 10*log10(S);
end

