function [results] = bcohCrit( magThresh, bcMatrix, fVec, T )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

results = [];

%k = 1;
verbose = 1;

maxFreq = max(fVec);

[m,n] = size(bcMatrix);

for k = 1:1
    mx = max(abs(bcMatrix(:)));
    mx = mx^2;
    for i = 1:n
        for j = 1:m
            %pdiff = 100 * (abs(fVec(ceil(j/2)) - fVec(ceil(i/2)))/((fVec(ceil(j/2))+fVec(ceil(i/2)))/2));
            pdiff = 100 * (abs(fVec(j) - fVec(i))/((fVec(j)+fVec(i))/2));
            %pdiff = 0.1;
            if ((abs(bcMatrix(j,i))^2 >= (magThresh)) && verbose)
                disp(strcat('Found peak!...',num2str(fVec(j)),'...',num2str(fVec(i)), ...
                    '...Magnitude:...',num2str(abs(bcMatrix(j,i)))))       
            end
            if (abs(bcMatrix(j,i))^2 >= (magThresh)) && (pdiff > 0.5)
                if ((fVec(j) + fVec(i)) > maxFreq)
                    if verbose
                        disp('WARNING: Frequency combination exceeds maximum allowed...')
                    end
                else            
                    disp('Found final QPC candidate!')
                    disp(strcat(num2str(fVec(j)),'...',num2str(fVec(i)), ...
                    '...Maximum magnitude=...',num2str(mx), ...
                    '...Magnitude threshold:...',num2str(magThresh),'...Magnitude:...',num2str(abs(bcMatrix(j,i)))))
                    if k == 1
                        results = vertcat(results,[fVec(j),fVec(i), ...
                        abs(bcMatrix(j,i)),mx,T,magThresh,j,i]);
                    end
                end
            end
        end
    end
end

end

