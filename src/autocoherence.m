function [ acoh,freq ] = autocoherence( y, nfft )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
nsamp = nfft;

if rem(nfft,2)
    divnfft = (nfft + 1)/2;
else
    divnfft = (nfft/2) + 1;
end

acoh    = zeros(divnfft,1);
yk    = 0;
ykj    = zeros(divnfft,1);

xseg = y;

% ----------------------------------------------------------------------
disp('Creating time-domain window...')

%wind = hanning(nsamp);

%wind = hann(nsamp,'periodic');

wind = chebwin(nsamp,200);
% ------------------apply window---------------------------------

xseg = (xseg(:) - mean(xseg)) .* wind;

% ---------------- accumulate double products----------------------
Xf = fft(xseg, nfft)/nsamp;
%Xf = fft(xseg-mean(xseg), nfft)/nsamp;  % original
Xf = Xf(1:divnfft);
Xf(2:end) = 2 * Xf(2:end);
        
CXf    = conj(Xf);

disp('Calculating double product:...')
            progBar = progressBar((divnfft)^2);
            for row = 1:divnfft
                yk = yk + (abs(Xf(row)))^2;
                for index = 1:divnfft
                    indx = index - row;
                    if indx < 1
                        indx = (-1) * indx;
                    end
                    indx = indx + 1;
                    %acoh(row) = acoh(row) + (Xf(index) * CXf(indx));  % individual double product for this experiment
		    acoh(row) = acoh(row) + (Xf(indx) * CXf(index));  % individual double product for this experiment
                    ykj(row) = ykj(row) + (abs(Xf(indx)))^2;
                    progBar(index + ((row - 1)*divnfft));
                end
            end

            acoh = abs(acoh);
            %acoh = acoh.^2;
	    yk
            ykj = yk * ykj;
            disp('Normalizing...')
            for index = 1:divnfft
                acoh(index) = acoh(index)/sqrt(ykj(index));
            end
            %den = yk * ykj;
            
            %acoh = acoh * (1/den);
            
            a = length(acoh);
            freq = 0:1/(a-1):1;
            freq = 0.5 * freq;
            
end
