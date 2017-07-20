function [results] = bpCritTheor( segments, nsamples, SNR, bcMatrix, fVec )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
M = segments;  % number of segments
N = nsamples;  % number of samples per segment (time series)
results = [];

[a,b] = size(SNR);

if (a > 1) || (b > 1)
    return;
end

varTheta = 1/(M*N*10^(SNR/10));

%k = 1;
verbose = 1;
% ~ 99.999
reqPercentile = 99.996;

%reqPercentile = 99.9;
maxFreq = max(fVec);
%sigmaTheta = (1/M)*((1/bc)-1);
%sigmaTheta = sqrt(abs(var(angle(bcMatrix(:)))));
%bPmean = mean(angle(bcMatrix(:)));

[m,n] = size(bcMatrix);

%Zc = 1.65;  % Number of standard deviations from the mean (theoretically 0)

%Zc = 2.05;
% alpha = 1 - 0.9798 = 0.0202
% increases the likelihood that an angle will fit

Zc = 2.10;
%alpha = 1 - 0.9821 = 0.0179 

c = Zc * sqrt(varTheta);
%ccalc = Zc * std(angle(bcMatrix(:)))
disp(strcat('Critical angle:...',num2str(c),'...radians'))
%Pfa = c/pi;
%cC = 0;
% ch1 = 0;
% ch2 = 0;
for k = 1:1
    mx = max(abs(bcMatrix(:)));
    %magThresh = ((mx)^2) * sqrt(9.2/(2*M));
    magThresh = prctile(abs(bcMatrix(:)),reqPercentile);  % require magnitude be in the 99th percentile
    %mx = 1;

    for i = 1:n
        for j = 1:m
            %cC = abs(Z(j,i))*sqrt((1/(2*M))*((1/bc)-1));
            cC = angle(bcMatrix(j,i));
            %pdiff = 100 * (abs(fVec(ceil(j/2)) - fVec(ceil(i/2)))/((fVec(ceil(j/2))+fVec(ceil(i/2)))/2));
            pdiff = 100 * (abs(fVec(j) - fVec(i))/((fVec(j)+fVec(i))/2));
            if ((abs(bcMatrix(j,i)) >= (magThresh)) && verbose)
                disp(strcat('Found peak!...',num2str(fVec(j)),'...',num2str(fVec(i)), ...
                    '...Magnitude:...',num2str(abs(bcMatrix(j,i)))))
                %ch1 = 1;        
            end
            if (abs(cC) < c) && (abs(bcMatrix(j,i)) >= (magThresh)) && (pdiff > 0.5)

                if ((fVec(j) + fVec(i)) > maxFreq)
                    if verbose
                        disp('WARNING: Frequency combination exceeds maximum allowed...')
                    end
                    %disp('Flattening...')
                    %bcMatrix(j,i) = 0;
                else            
                    disp('Found QPC candidate!')
                    disp(strcat(num2str(fVec(j)),'...',num2str(fVec(i)),'...Angle=...', ...
                    num2str(cC),'...Maximum magnitude=...',num2str(mx), ...
                    '...Magnitude threshold:...',num2str(magThresh),'...Magnitude:...',num2str(abs(bcMatrix(j,i)))))
                    %ch2 = 1;
                    if k == 1
                        results = vertcat(results,[fVec(j),fVec(i), ...
                        cC,abs(bcMatrix(j,i)),mx,magThresh,j,i]);
                    end
                    %k = k + 1;
                end
            end
%             if (ch1 == 1) && (ch2 == 0)
%                 disp('Flattening...')
%                 bcMatrix(j,i) = 0;
%             end
%             ch1 = 0;
%             ch2 = 0;
        end
    end
end
Pfa = c/pi;
disp(strcat('(Angular) False Alarm Probability:...',num2str(Pfa)))

end

