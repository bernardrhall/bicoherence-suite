function [ T ] = bspcThresh( segments, samples, sig, pr, fs, verbose, SNRdb )
%This function calculates T_min, the bicoherence magnitude squared
%threshold

%   Detailed explanation goes here
M = segments;  % number of segments
N = samples;  % number of (time series) samples per segment
%SNR = snr(sig); % gives snr in decibels
if ~SNRdb
    [~,SNR,~,~] = Snr_1(sig,fs,0,verbose); % customized SNR calculation function
else
    SNR = SNRdb;
end

varTheta = 1/(M*N*10^(SNR/10)); % expression in denominator converts snr
%back into a ratio of sigma_sig/sigma_noise

varR = var(pr);
%varR = 0.5;

sigmaTheta = sqrt(varTheta);
%sigmaR = sqrt(varR);

%T = vpa((sigmaTheta/sqrt(2*pi))^(2*varR*varTheta),32);
T = (sigmaTheta/sqrt(2*pi))^(2*varR*varTheta);
%T = char(T);  % output as string to prevent roundoff errors
T = num2str(T);
%T = T(1:18);
%T = str2double(T);
end

