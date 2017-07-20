function [combined] = bpCrit( segments, bcMatrix, fVec, T )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
M = segments;  % number of segments
%N = samples;  % number of samples per segment
%bc = 0.9;
%digits(32);
% Snr = vpa(10^(S/10),32);
% %correction = 0.99999999999999999999999999867551;
% bc = vpa(1 - (1/Snr),32);
% bc = char(bc);
bc = T;
if length(bc) > 17
	bc = bc(1:18);
end
bc = str2double(bc);
%bc = 0.99999999999;
bc1 = 0.99;
bc2 = 0.9;
bc3 = 0.8;
%bc = 0.7;
results = [];
results1 = [];
results2 = [];
results3 = [];
%k = 1;
verbose = 1;
% ~ 99.999
%reqPercentile = 99.999;
reqPercentile = 99.996;
maxFreq = max(fVec);
%sigmaTheta = (1/M)*((1/bc)-1);
%sigmaTheta = sqrt(abs(var(angle(bcMatrix(:)))));
%bPmean = mean(angle(bcMatrix(:)));

[m,n] = size(bcMatrix);

% if length(fVec) < n
%    if m == n
%        tmpVec = zeros(n,1);
%        for i = 1:n
%            if rem(i,2) == 0 
%               tmpVec(i) = (fVec((i/2) + 1)+fVec(i/2))/2;
%               %tmpVec(i + 1) = fVec((i/2) + 1);
%            else
%               tmpVec(i) = fVec(ceil(i/2));
%            end
%        end
%        fVec = tmpVec
%        clearvars tmpVec
%    else
%        disp('Not a square matrix!!')
%    end
% end


%Z = (1/sigmaTheta) * (angle(bcMatrix) - bPmean);
%disp('Calculating angular Z-scores')
%Z = zscore(angle(bcMatrix));

%bcMatrix = angle(bcMatrix);

Zc = 1.65;

%magThresh = sqrt(9.2/(2*M));
c = Zc * sqrt((1/(2*M))*((1/bc)-1));
c1 = Zc * sqrt((1/(2*M))*((1/bc1)-1));
c2 = Zc * sqrt((1/(2*M))*((1/bc2)-1));
c3 = Zc * sqrt((1/(2*M))*((1/bc3)-1));
%ccalc = Zc * std(angle(bcMatrix(:)))
disp(strcat('Critical angle 1:...',num2str(c),'...radians'))
disp(strcat('Critical angle 2:...',num2str(c1),'...radians'))
disp(strcat('Critical angle 3:...',num2str(c2),'...radians'))
disp(strcat('Critical angle 4:...',num2str(c3),'...radians'))
%Pfa = c/pi;
%cC = 0;

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
            %disp(num2str(pdiff))
        end
        if (abs(cC) < c) && (abs(bcMatrix(j,i)) >= (magThresh)) && (pdiff > 0)  % 1.8 ?
        %pdiff = 0.5 normally    
            if ((fVec(j) + fVec(i)) > maxFreq)
                if verbose
                    disp('WARNING: Frequency combination exceeds maximum allowed...')
                end
            else            
                disp('Found primary QPC candidate!')
                disp(strcat(num2str(fVec(j)),'...',num2str(fVec(i)),'...Angle=...', ...
                num2str(cC),'...Maximum magnitude=...',num2str(mx), ...
                '...Magnitude threshold:...',num2str(magThresh),'...Magnitude:...',num2str(abs(bcMatrix(j,i)))))

                results = vertcat(results,[fVec(j),fVec(i), ...
                        cC,abs(bcMatrix(j,i)),mx,magThresh,j,i]);
                %k = k + 1;
            end
        elseif (abs(cC) < c1) && (abs(bcMatrix(j,i)) >= (magThresh)) && (pdiff > 0)
            if ((fVec(j) + fVec(i)) > maxFreq)
                if verbose
                    disp('WARNING: Frequency combination exceeds maximum allowed...')
                end
            else            
                disp('Found secondary QPC candidate!')
                disp(strcat(num2str(fVec(j)),'...',num2str(fVec(i)),'...Angle=...', ...
                num2str(cC),'...Maximum magnitude=...',num2str(mx), ...
                '...Magnitude threshold:...',num2str(magThresh),'...Magnitude:...',num2str(abs(bcMatrix(j,i)))))

                results1 = vertcat(results1,[fVec(j),fVec(i), ...
                        cC,abs(bcMatrix(j,i)),mx,magThresh,j,i]);
            end
        elseif (abs(cC) < c2) && (abs(bcMatrix(j,i)) >= (magThresh)) && (pdiff > 0)
            if ((fVec(j) + fVec(i)) > maxFreq)
                if verbose
                    disp('WARNING: Frequency combination exceeds maximum allowed...')
                end
            else            
                disp('Found tertiary QPC candidate!')
                disp(strcat(num2str(fVec(j)),'...',num2str(fVec(i)),'...Angle=...', ...
                num2str(cC),'...Maximum magnitude=...',num2str(mx), ...
                '...Magnitude threshold:...',num2str(magThresh),'...Magnitude:...',num2str(abs(bcMatrix(j,i)))))

                results2 = vertcat(results2,[fVec(j),fVec(i), ...
                        cC,abs(bcMatrix(j,i)),mx,magThresh,j,i]);
            end
        elseif (abs(cC) < c3) && (abs(bcMatrix(j,i)) >= (magThresh)) && (pdiff > 0)
            if ((fVec(j) + fVec(i)) > maxFreq)
                if verbose
                    disp('WARNING: Frequency combination exceeds maximum allowed...')
                end
            else            
                disp('Found outlying QPC candidate...')
                disp(strcat(num2str(fVec(j)),'...',num2str(fVec(i)),'...Angle=...', ...
                num2str(cC),'...Maximum magnitude=...',num2str(mx), ...
                '...Magnitude threshold:...',num2str(magThresh),'...Magnitude:...',num2str(abs(bcMatrix(j,i)))))

                results3 = vertcat(results3,[fVec(j),fVec(i), ...
                        cC,abs(bcMatrix(j,i)),mx,magThresh,j,i]);
            end
        end
    end
end

combined = {results results1 results2 results3};

Pfa = c/pi;
Pfa1 = c1/pi;
Pfa2 = c2/pi;
Pfa3 = c3/pi;

disp(strcat('(Angular) False Alarm Probability (primary):...',num2str(Pfa)))
disp(strcat('(Angular) False Alarm Probability (secondary):...',num2str(Pfa1)))
disp(strcat('(Angular) False Alarm Probability (tertiary):...',num2str(Pfa2)))
disp(strcat('(Angular) False Alarm Probability (outlier):...',num2str(Pfa3)))

end

