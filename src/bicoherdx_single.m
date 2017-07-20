function [Bspec,waxis,randVec] = bicoherdx_single (w, x, y,  plot, randomize, rvector, tile, nfft, wind, nsamp, overlap)
%BISPECD Bispectrum estimation using the direct (fft-based) approach.
%	[Bspec,waxis] = bispecd (y,  nfft, wind, segsamp, overlap)
%	y    - data vector or time-series
%	nfft - fft length [default = power of two > segsamp]
%	wind - specifies the time-domain window to be applied to each
%	       data segment; should be of length 'segsamp' (see below);
%		   otherwise, the default Hanning window is used.
%	segsamp - samples per segment [default: such that we have 8 segments]
%	        - if y is a matrix, segsamp is set to the number of rows
%	overlap - percentage overlap [default = 50]
%	        - if y is a matrix, overlap is set to 0.
%
%	Bspec   - estimated bispectrum: an nfft x nfft array, with origin
%	          at the center, and axes pointing down and to the right.
%	waxis   - vector of frequencies associated with the rows and columns
%	          of Bspec;  sampling frequentiledcy is assumed to be 1.
%   plot    - boolean value to set whether or not to create plot (0 = no,
%             1 = yes)  %added
% randomize - boolean to designate whether or not to use phase
%             randomization (0 = no, 1 = yes)  %added

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc.
%       $Revision: 1.8 $
%  A. Swami   January 20, 1993.

%     RESTRICTED RIGHTS LEGEND
% Use, duplication, or disclosure by the Government is subject to
% restrictions as set forth in subparagraph (c) (1) (ii) of the
% Rights in Technical Data and Computer Software clause of DFARS
% 252.227-7013.
% Manufacturer: United Signals & Systems, Inc., P.O. Box 2374,
% Culver City, California 90231.
%
%  This material may be reproduced by or for the U.S. Government pursuant
%  to the copyright license under the clause at DFARS 252.227-7013.

% --------------------- parameter checks -----------------------------

    [ly, nrecs] = size(y);
    if (ly == 1) y = y(:);  ly = nrecs; nrecs = 1; end

    if (exist('nfft') ~= 1)            nfft = 128; end
    if (exist('overlap') ~= 1)      overlap = 50;  end
    overlap = min(99.9,max(overlap,0));
    if (nrecs > 1)                  overlap =  0;  end
    if (exist('nsamp') ~= 1)          nsamp = 0;   end
    if (nrecs > 1)                    nsamp = ly;  end

    if (nrecs == 1 & nsamp <= 0)
       nsamp = fix(ly/ (8 - 7 * overlap/100));
    end
    if (nfft  < nsamp)   nfft = 2^nextpow2(nsamp); end

    overlap  = fix(nsamp * overlap / 100);             % added 2/14
    nadvance = nsamp - overlap;
    nrecs    = fix ( (ly*nrecs - overlap) / nadvance);

% ----------------------------------------------------------------------

    if rem(nfft,2)
        divnfft = (nfft + 1)/2;
    else
        divnfft = (nfft/2) + 1;
    end

% ----------------------------------------------------------------------
    disp('Creating time-domain window...')

    if (exist('wind') ~= 1) 
        wind = hanning(nsamp); 
    end
    
    [rw,cw] = size(wind);
    if (min(rw,cw) ~= 1 || max(rw,cw) ~= nsamp)
	   disp(['Segment size  is ',int2str(nsamp)])
	   disp(['"wind" array  is ',int2str(rw),' by ',int2str(cw)])
	   disp('Using default Hanning window')
	   wind = hanning(nsamp);
    end
    wind = wind(:);
% ---------------- accumulate triple products ----------------------

    rng('shuffle')  %reseed based on time

    Bspec    = zeros(divnfft,divnfft);
    prspec    = zeros(divnfft,divnfft);
    %normalI = zeros((nfft/2)+1,(nfft/2)+1);
    %normalS = zeros((nfft/2)+1,(nfft/2)+1);
    normalization = zeros(divnfft,divnfft);
    %mask = hankel([1:nfft],[nfft,1:nfft-1] );   % the hankel mask (faster)
    
    locseg = [1:nsamp]';
    
    %tileSeg = 1;  %0 to disable, 1 to enable
    tileSeg = tile;  %0 to disable, 1 to enable
    tileFlag = 1;
    
    if tileSeg
        prspecT    = zeros(divnfft,divnfft);
    end
    
%    progBar = progressBar(2*nfft);
%      for tmp = 1:n
%        progBar(tmp);

     if randomize == 1
        randVec = zeros(nrecs,1);
     else
        randVec = 0;
     end
    
    for krec = 1:nrecs
        
        yseg   = y(locseg);
        yseg = (yseg(:) - mean(yseg)) .* wind;% original "detrending" method/syntax
        %xseg = detrend(xseg,0);  %not sure if we want to detrend or not;
        %                          may lose low-frequency features?
        wseg   = w(locseg);
        wseg = (wseg(:) - mean(wseg)) .* wind;
        
        xseg   = x(locseg);
        xseg = (xseg(:) - mean(xseg)) .* wind;
        
        Wf     = fft(wseg, nfft)/nsamp;
        Xf     = fft(xseg, nfft)/nsamp;
        Yf     = fft(yseg, nfft)/nsamp;
        
        %Xf     = fft(xseg, nfft)/nsamp;  % no detrending
        
        Wf = Wf(1:divnfft);
        Wf(2:end) = 2 * Wf(2:end);
        
        Xf = Xf(1:divnfft);
        Xf(2:end) = 2 * Xf(2:end);
        
        Yf = Yf(1:divnfft);
        Yf(2:end) = 2 * Yf(2:end);
        %------------------------------------  break up into smaller chunks
        %disp('Creating Hankel Mask...')
        %mask = hankel( [1:nfft],[nfft,1:(nfft)-1] );   % the hankel mask (faster)
        %CWf    = conj(Wf);
        %CXf    = conj(Xf);
        CYf    = conj(Yf);
        %CXf    = [flipud(CXf(2:end));CXf];
        %size(mask)
        %[m,n] = size(prspec);
        %disp('Creating conjugate matrix')
        %CXf = CXf(mask);
        
        %clearvars mask  % clear some space before next operations
        
        %disp('Reshaping conjugate matrix')
        %prconj = reshape(CXf, nfft, nfft)  % is this necessary?  The martix sould already be the right size
        
        %clearvars CXf
        
        %disp('Taking outer product of Xf...')
        %prspec = (Xf * Xf.');  % outer product
        
        if tileFlag
            disp(strcat('Calculating triple product:...',num2str(krec),'...of...',num2str(nrecs)))
            progBar = progressBar((divnfft)^2);
            for col = 1:divnfft
                for row = 1:divnfft
                    indx = (row + col) - 1;
                    if indx > divnfft  % lower "anti-quadrant" of Hankel Matrix
                        indx = indx - divnfft;
                    end
                    %XfC = XfC/abs(XfC);
                    prspec(row,col) = Wf(col) * Xf(row) * CYf(indx);  % individual triple product for this experiment
                    
                    %normalization(row,col) = normalization(row,col) + abs(prspec(row,col)); 
                    
                    % Not taking magnitude squared here, since we undo
                    % that operation later in the algorithm.  It should
                    % be fine, as long as we remember to account for this.
                    %progBar(row*col);
                    progBar(row + ((col - 1)*(divnfft)));
                end
            end
        else
            disp(strcat('Tiled, cycle:...',num2str(krec),'...of...',num2str(nrecs)))
            prspec = prspecT;
            %normalization = normalization + normT;
        end
        
        if tileSeg
            if tileFlag
                prspecT = prspec;
                %normT = normalization;
            end
            tileFlag = 0;
        end
        
        %disp('Multiplying the matrices...')
        %prspec = prspec .*prconj; % element-wise multiplication of two arrays
        
        %clearvars prconj
        %-----------------------------------
        
        %CXf    = conj(Xf);
        
        %prspec = (Xf * Xf.') .*reshape(CXf(mask), nfft, nfft);
        
        if randomize == 1
            disp('Randomizing the biphase for this experiment')
            
            m=abs(prspec);
            % Angles
            p=angle(prspec);
            % The imaginary unit
            i=sqrt(-1);
            
            pr = -pi;  % initialize and make sure 1 iteration of the loop is completed
            
            while pr == -pi % the value of -pi is not included in the interval
                            % since exp(i*pi) = exp(-i*pi).
                            % The expectation value of a uniform
                            % distribution of complex numbers with
                            % phase/argument (-pi, pi] is 0.
                if ~rvector
                    disp('(Using MATLAB rand)')
                    pr = -pi + (pi+pi)*rand;  %multiplicative constant; uniform distribution between -pi and pi
                else
                    disp('(Using user-provided random numbers)')
                    pr = -pi + (pi+pi)*rvector(krec);
                    if pr == -pi
                        pr = pr + 0.0001;
                    end
                end
            end
            
            randVec(krec) = pr;
            disp(strcat('R =...',num2str(pr)))
            
            %pr = -pi + (pi+pi)*randn;  %multiplicative constant; gaussian distribution between -pi and pi
        
            %pr = 1 + ((2*pi)-1)*rand; %multiplicative constant; uniform distribution between 1 and 2pi
                                      % * We may not want to include zero
                                      % in the range of this variable.  If
                                      % the biphase is zero, it will not
                                      % matter what we multiply it by.  If
                                      % it is not, then we do not want to
                                      % artificially make it so.
            %----------------------------
%             progBar = progressBar(nfft);
%             for col = 1:nfft/2
%                 for row = 1:nfft/2
%                     m=abs(prspec(row,col));
%                     % Angles
%                     p=angle(prspec(row,col));
%                     % Back to the complex numbers
%                     prspec(row,col) = m.*exp(i*(p*pr));
%                     progBar(row+col);
%                 end
%             end
            
            prspec = m.*exp(i*(p*pr));
            clearvars m p
        end
        
        disp('Adding to total (bispectrum)...')
%         progBar = progressBar(nfft);
%         
%         for col = 1:nfft/2
%             for row = 1:nfft/2
%                 Bspec(row,col) = Bspec(row,col) + prspec(row,col);
%                 progBar(row+col);
%             end
%         end
        % see https://en.wikipedia.org/wiki/Bicoherence
        normalization = normalization + abs(prspec);  % ~ perfect biphase
        Bspec = Bspec + prspec; % ~ calculated biphase
        
        %Bspec  = Bspec + (Xf * Xf.') .* ...
	    %     reshape(CXf(mask), nfft, nfft);
         
        locseg = locseg + nadvance;
    end
    
    %clearvars prspec normT
    clearvars prspec
    %disp('Normalizing bispectrum...')
    %Bspec = fftshift(Bspec)/(nrecs); %Normalization: 'M'
    %Bspec = Bspec/(nrecs); %Normalization: 'M'
    
%     disp('Normalizing first term of denominator...')
%     
%     normalI = ((1/nrecs) * normalI);
%     
%     disp('Normalizing second term of denominator...')
%     
%     normalS = ((1/nrecs) * normalS);
%     
    
    disp('Calculating bicoherence (normalizing bispectrum)...')
    
    %Bspec = Bspec ./ normalization;
    
    progBar = progressBar((divnfft)^2);
    for col = 1:divnfft
        for row = 1:divnfft            
            Bspec(row,col) = Bspec(row,col)/normalization(row,col);
            progBar(row + ((col - 1)*(divnfft)));
        end
    end
    
%     m = abs(Bspec);
%     n = angle(Bspec);
%     i = sqrt(-1);
%     
%     progBar = progressBar((nfft/2 + 1)^2);
%     for col = 1: ((nfft/2) + 1)
%         for row = 1: ((nfft/2)+1)            
%             m(row,col) = m(row,col)/normalization(row,col);
%             progBar(row*col);
%         end
%     end
%     
%     Bspec = m.*exp(i*n);
%     
%     clearvars m n
    
% ------------ contour plot of magnitude bispectum --------------------

   [a,~] = size(Bspec);
   %waxis = 0:1/a:1;
   waxis = 0:1/(a-1):1;
   waxis = 0.5 * waxis;
   
   if plot == 1
%         hold off, clf
% %  contour(abs(Bspec),4,waxis,waxis),grid
%         contour(waxis,waxis,abs(Bspec),4),grid on 
%         title('Bispectrum estimated via the direct (FFT) method')
%         xlabel('f1'), ylabel('f2')
%         set(gcf,'Name','Hosa BISPECD')
   end
return
