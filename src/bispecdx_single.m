function [Bspec,waxis] = bispecdx_single (x, y, z,  plot, randomize, nfft, wind, nsamp, overlap)
%BISPECD Bispectrum estimation using the direct (fft-based) approach.
%	[Bspec,waxis] = bispecd (y,  nfft, wind, segsamp, overlap)
%	y    - data vector or time-series
%	nfft - fft length [default = power of two > segsamp]
%	wind - window specification for frequency-domain smoothing
%	       if 'wind' is a scalar, it specifies the length of the side
%	          of the square for the Rao-Gabr optimal window  [default=5]
%	       if 'wind' is a vector, a 2D window will be calculated via
%	          w2(i,j) = wind(i) * wind(j) * wind(i+j)
%	       if 'wind' is a matrix, it specifies the 2-D filter directly
%	segsamp - samples per segment [default: such that we have 8 segments]
%	        - if y is a matrix, segsamp is set to the number of rows
%	overlap - percentage overlap [default = 50]
%	        - if y is a matrix, overlap is set to 0.
%
%	Bspec   - estimated bispectrum: an nfft x nfft array, with origin
%	          at the center, and axes pointing down and to the right.
%	waxis   - vector of frequencies associated with the rows and columns
%	          of Bspec;  sampling frequency is assumed to be 1.
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
    overlap = min(99,max(overlap,0));
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


% ------------------- create the 2-D window -------------------------
  if (exist('wind') ~= 1) wind = 5; end
  [m,n] = size(wind);
  window = wind;
  if (max(m,n) == 1)     % scalar: wind is size of Rao-Gabr window
     winsize = wind;
     if (winsize < 0) winsize = 5; end        % the window length L
     winsize = winsize - rem(winsize,2) + 1;  % make it odd
     if (winsize > 1)
        mwind   = fix (((nfft/2)+1)/winsize);            % the scale parameter M
        lby2    = (winsize - 1)/2;

        theta  = -lby2:lby2;
        opwind = ones(winsize,1) * (theta .^2);       % w(m,n)=m^2
        opwind = opwind + opwind' + theta' * theta;   % m^2 + n^2 + mn
        opwind = 1 - (2*mwind/((nfft/2)+1))^2 * opwind;       %
        hex    = ones(winsize,1) * theta;             % m
        hex    = abs(hex) + abs(hex') + abs(hex+hex');
        hex    = (hex < winsize);
        opwind = opwind .* hex;
        opwind = opwind * (4 * mwind^2) / (7 * pi^2) ;
     else
        opwind = 1;
     end

  elseif (min(m,n) == 1)  % 1-D window passed: convert to 2-D
     window = window(:);
     if (any(imag(window) ~= 0))
        disp(['1-D window has imaginary components: window ignored'])
        window = 1;
     end
     if (any(window < 0))
        disp(['1-D window has negative components: window ignored'])
        window = 1;
     end
     lwind  = length(window);
     windf  = [window(lwind:-1:2); window];    % the full symmetric 1-D
     window = [window; zeros(lwind-1,1)];
     opwind = (windf * windf')      ...
              .* hankel(flipud(window), window); % w(m)w(n)w(m+n)
     winsize = length(window);

  else                    % 2-D window passed: use directly
    winsize = m;
    if (m ~= n)
       disp('2-D window is not square: window ignored')
       window = 1;
       winsize = m;
    end
    if (rem(m,2) == 0)
       disp('2-D window does not have odd length: window ignored')
       window = 1;
       winsize = m;
    end
    opwind  = window;
  end

% ---------------- accumulate triple products ----------------------

    rng('shuffle')  %reseed based on time

    Bspec    = zeros((nfft/2)+1,(nfft/2)+1);
    prspec    = zeros((nfft/2)+1,(nfft/2)+1);
    %mask = hankel([1:nfft],[nfft,1:nfft-1] );   % the hankel mask (faster)
    
    locseg = [1:nsamp]';
    
    tileSeg = 1;  %0 to disable, 1 to enable
    tileFlag = 1;
    
    if tileSeg
        prspecT    = zeros((nfft/2)+1,(nfft/2)+1);
    end
    
%    progBar = progressBar(2*nfft);
%      for tmp = 1:n
%        progBar(tmp);
    
    nrecs = 8;  %temp

    for krec = 1:nrecs
        
        xseg   = x(locseg);
        yseg   = y(locseg);
        zseg   = z(locseg);
        %xseg = detrend(xseg,0);  %not sure if we want to detrend or not;
        %                          may lose low-frequency features?
        Xf     = fft(xseg-mean(xseg), nfft)/nsamp;  % original
        %                                              "detrending" method/syntax
        Yf     = fft(xseg-mean(yseg), nfft)/nsamp;  % original
        Zf     = fft(xseg-mean(zseg), nfft)/nsamp;  % original
        
        %Xf     = fft(xseg, nfft)/nsamp;  % no detrending
        
        Xf = Xf(1:(nfft/2)+1);
        Xf(2:end) = 2 * Xf(2:end);
        
        Yf = Yf(1:(nfft/2)+1);
        Yf(2:end) = 2 * Yf(2:end);
        
        Zf = Zf(1:(nfft/2)+1);
        Zf(2:end) = 2 * Zf(2:end);
        
        %------------------------------------  break up into smaller chunks
        %disp('Creating Hankel Mask...')
        %mask = hankel( [1:nfft],[nfft,1:(nfft)-1] );   % the hankel mask (faster)
        CZf    = conj(Zf);
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
            progBar = progressBar(nfft);
            for col = 1:(nfft/2)+1
                for row = 1:(nfft/2)+1
                    indx = (row + col) - 1;
                    if indx > (nfft/2) + 1  % lower "anti-quadrant" of Hankel Matrix
                        indx = indx - ((nfft/2)+1);
                    end
                    prspec(row,col) = Xf(col) * Yf(row) * CZf(indx);  % individual triple product for this experiment
                    progBar(row+col);
                end
            end
        else
            disp(strcat('Tiled, cycle:...',num2str(krec),'...of...',num2str(nrecs)))
            prspec = prspecT;
        end
        
        if tileSeg
            if tileFlag
                prspecT = prspec;
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
            
            %pr = -pi + (pi+pi)*rand;  %multiplicative constant; uniform distribution between -pi and pi
            %pr = -pi + (pi+pi)*randn;  %multiplicative constant; gaussian distribution between -pi and pi
        
            pr = 1 + ((2*pi)-1)*rand; %multiplicative constant; uniform distribution between 1 and 2pi
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
        end
        
        clearvars m p
        
        disp('Adding to total...')
%         progBar = progressBar(nfft);
%         
%         for col = 1:nfft/2
%             for row = 1:nfft/2
%                 Bspec(row,col) = Bspec(row,col) + prspec(row,col);
%                 progBar(row+col);
%             end
%         end

        Bspec = Bspec + prspec;
        
        %Bspec  = Bspec + (Xf * Xf.') .* ...
	    %     reshape(CXf(mask), nfft, nfft);
         if tileFlag
            locseg = locseg + nadvance;
         end
    end
    
    clearvars prspec
    disp('Normalizing...')
    %Bspec = fftshift(Bspec)/(nrecs); %Normalization: 'M'
    Bspec = Bspec/(nrecs); %Normalization: 'M'
% ----------------- frequency-domain smoothing ------------------------
    disp('Applying window...')
  if (winsize > 1)
      lby2 = (winsize-1)/2;
      Bspec = conv2(Bspec,opwind);
      Bspec = Bspec(lby2+1:lby2+((nfft/2)+1),lby2+1:lby2+((nfft/2)+1));
  end
% ------------ contour plot of magnitude bispectum --------------------

   if (rem(nfft,2) == 0)
       waxis = [0:((nfft/2)+1)/2]'/((nfft/2)+1);
   else
       waxis = [0:((nfft+1)/2)/2]'/((nfft+1)/2);
   end
   
   if plot == 1
        hold off, clf
%  contour(abs(Bspec),4,waxis,waxis),grid
        contour(waxis,waxis,abs(Bspec),4),grid on 
        title('Bispectrum estimated via the direct (FFT) method')
        xlabel('f1'), ylabel('f2')
        set(gcf,'Name','Hosa BISPECD')
   end
return
