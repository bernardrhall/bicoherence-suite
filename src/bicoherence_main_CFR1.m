function bicoherence_main_CF( GW_channel, fs, threshMan, HbicohThresh,HbicohBspcThresh, fnct, mode, filt, ChannelA, ChannelB, segslice, segStart, segEnd, full_SNR, autom, time_window, offset_multiplier, ol, seg, tile, GPS_central, check, dnsmple, decm, plot_img, randomize, preSeq, sequenceN, fsd, pass_band_fl, pass_band_fh, uSnrV, upr, verbose, chStr, path, ifoCh )

% -------------------------Control Panel----------------------------------
% ------------------------------------------------------------------------
fs = str2double(fs);
threshMan = str2num(threshMan);
HbicohThresh = str2double(HbicohThresh);
HbicohBspcThresh = str2double(HbicohBspcThresh);
fnct = str2num(fnct);
mode = str2num(mode);
filt = str2num(filt);
segslice = str2num(segslice);
segStart = str2double(segStart);
segEnd = str2double(segEnd);
full_SNR = str2num(full_SNR);
autom = str2num(autom);
time_window = str2double(time_window);
offset_multiplier  = str2double(offset_multiplier);
ol = str2double(ol);
seg = str2num(seg);
tile = str2num(tile);
GPS_central = str2double(GPS_central);
check = str2num(check);
dnsmple = str2num(dnsmple);
decm = str2num(decm);
plot_img = str2num(plot_img);
randomize = str2num(randomize);
preSeq = str2num(preSeq);
fsd = str2double(fsd);
pass_band_fl = str2double(pass_band_fl);
pass_band_fh = str2double(pass_band_fh);
uSnrV = str2double(uSnrV);
upr = str2double(upr);
%verbose
%chStr

% GW input channel ************************

%GW_channel = 'H1-CAL-DELTAL_EXTERNAL_DQ_112625944630.csv';

wvlt = 0; % experimental (0 is off)
%fs = 16384; %(original) sample frequency; set to 0 will auto-calculate

% ----------------------------------------------------------
%fnct = 6; % 0: cross bispectrum, 1: cross bicoherence 2: auto bicoherence, 
%3: auto bispectrum, 
%4: magnitude squared coherence (between wigh band and low band),
%5:autocoherence, 6: bispectrogram, 7: time-slide coherence

% manual upper limit set for bispectrogram
amax = 0;  % 0 for manual, 1 for auto (normal)
%upr = 17*10^(-42);
% ----------------------------------------------------------
if fnct ~= 7
%    mode = 0;
    % With the auto bicoherence (fnct = 2): 0 for periodic Hann window;
    % 1 for Chebyshev; with 'tiled': 2 for default window with nfft =
    % seg_length; 3 for default window with nfft = 16384 (will zero-pad if
    % necessary)
    % With autobispectrum, 1 uses seg-length for FFT; 0 uses 16384 for FFT
    % length;
    % 
    % With function 3, 0 uses static fft length of 16384 (zero-padding, if necessary);
    %   1 uses time segment length for fft
    % ----------------------------------------------------------
    % optional
    %filt = 1; %set to 1 to enable.  Enables filter for "auto" cross bicoherence
    slow = 0; % experimental, "Slow motion" setting (reduces the reported
    % sampling rate, i.e., to 'fsd')

    % X, Y Channel input **********************
    %if (fnct == 0 || fnct == 1) && filt == 0
    %    ChannelA = 'TCS-ITMY_CO2_ISS_IN_AC_OUT_DQ.csv';
    %    ChannelB = 'ASC-AS_A_RF45_Q_PIT_OUT_DQ.csv';
    %end
    %----------------------------------------------

    % if tiled, this is the length of an individual experiment
    % if not tiled, this is the entire length of the segment to be analyzed.
    % This segment will be divided into 'seg' experiments, the length of which
    % will be calculated based on the length of the segment divided by the
    % number of experiments.  (No overlap is default.)  A value of '-1' uses
    % the entire loaded time series.

    %segment slice-----------------
    % Reduce the input GW time series to time between
    % 'segStart' and 'segEnd'
    %segslice = 0; % 1 to enable, 0 to disable
    %segStart = 2.0; %seconds
    %segEnd = 6.0;

    %full_SNR = 1; % calculate SNR before slicing; 1 recommended
    uSnr = 1; % use manual (user) SNR input: 0 -> off, 1 -> on
    %uSnrV = 3.507599869825159; % user provided SNR value 
	
    %autom = 0; % 1 for autoset, 0 for manual set (below)

    %----------------------------------manual set
    if ~autom % use with manual set
        if fnct == 1 || fnct == 2
            disp(strcat('Time window for bicoherence:...',num2str(time_window),'...s'))
        elseif fnct == 0 || fnct == 3
            disp(strcat('Time window for bispectrum:...',num2str(time_window),'...s'))
        else
            disp(strcat('Time window for study:...',num2str(time_window),'...s'))
        end
    else % use with auto set
        if fnct == 1 || fnct == 2
            disp(strcat('Time window for bicoherence:...',num2str(time_window),'...s'))
        else
            disp(strcat('Time window for bispectrum:...',num2str(time_window),'...s'))
        end
    end
    %----------------------------------

    %------------------------------
    % Select time around GW time series central time on which to
    % perform analysis
    %time_window = 8; %seconds; "-1" to use entire time series
    if segslice
        if time_window > (segEnd - segStart)
                time_window = segEnd - segStart;
        end
    end

    %ol = 99.9;  %overlap amount
    %ol = 98.0;
    %------------------------manual segmentation (use no tile)
    % Note that the actual number of segments may be greater than
    % 'seg' when using 'no tile' option in conjunction with an overlap > 0.
    % Use 'offset_multiplier' to compensate for this, if desired.
    %seg = 30;  % number of time segments (experiments) to use in analysis
    %tile = 1;  % 0 to disable, 1 to enable
    %---------------------------------------

    % Expected GPS Central
    %GPS_central = 1114677781;

    %------------------------------------
    % sanity check to make sure the central time you think you have
    % is actually the time you have

    %check = 1;   % 1 to enable, 0 to disable
    %-------------------------

    %dnsmple = 0; %set to 1 to enable; use with caution (may cause aliasing)
    %decm = 1; % set to 1 to enable
    %plot_img = 1; % Create plot?  (1 for yes, 0 for no)
    %randomize = 1; %1 enables phase randomization, 0 disables
    % set preSeq to 0 to use MATLAB's rand() for generation -----
    %preSeq = 1;

    if preSeq
        %sequenceN = 'seq_552015_2.txt'; % 
        %sequenceN = 'seq_552015_6.txt';
        rSequence = load(sequenceN);
    else
        rSequence = 0;
    end
    %----------------------------------------------------

    %fsd = 8192; % reduced sampling frequency. Used if downsample or decimate is selected
    %if (fnct == 0 || fnct == 1 || fnct == 4) && filt == 1
    %    pass_band_fl = 20;  %frequency in Hz
    %    pass_band_fh = 20;
    %end

    %---------------------------------------------------------------------
    %-----------------End Control Panel------------------------------------

    % ----------------Process Input----------------------------------------
    GWc = load(GW_channel);

    disp('Channel loaded...')

    if full_SNR && ~(fnct==6)% calculate before slicing
	if ~uSnr
        	disp('Calculating SNR (full)...')
        	[~,Snr,~,~] = Snr_1(GWc(:,2),fs,0,1); % calculate SNR on raw signal
        	Sr = 10^(Snr/10);
	end
    end
    
    if (fnct == 0 || fnct == 1) && filt == 0
        chA = load(ChannelA);
        chB = load(ChannelB);
        GW_high = chA(:,2);
        GW_low = chB(:,2);

        % make sure all the sampling rates are the same
	fsCheckM = floor(length(GWc(:,1))/(GWc(length(GWc(:,1))-1,1) - GWc(1,1)));
	if fsCheckM ~= fs
		disp('Warning: Resetting GW supplied sampling rate...')
		fs = fsCheckM;
	end
        fsCheck1 = round(length(chA(:,1))/(chA(length(chA(:,1))-1,1) - chA(1,1)));
        fsCheck2 = round(length(chB(:,1))/(chB(length(chB(:,1))-1,1) - chB(1,1)));

        %fsCheck1 = fs; %%%%% temp quick fix
        %fsCheck2 = fs;%%%%
	if fsCheckM < fsCheck1
            GW_high = decimate(GW_high,(round(fsCheck1/fsCheckM)));
            fsCheckM = fsCheck1;
	    fs = fsCheckM;
        end

	if fsCheckM < fsCheck2
            GW_low = decimate(GW_low,(round(fsCheck2/fsCheckM)));
            fsCheckM = fsCheck2;
            fs = fsCheckM;
        end

        if fsCheck1 > fsCheck2
            GW_high = decimate(GW_high,(round(fsCheck1/fsCheck2)));
            fsCheck1 = fsCheck2;
        end
        if fsCheck2 > fsCheck1
            GW_low = decimate(GW_low,(round(fsCheck2/fsCheck1)));
            fsCheck2 = fsCheck1;
        end
        if fs > fsCheck1
            A = GWc(1,1)
            B = GWc(end,1)
            %GW_1c = decimate(GW_1c,(round(fs/fsCheck1)));
            GW1 = decimate(GWc(:,2),(round(fs/fsCheck1)));
            %GW2 = decimate(GWc(:,1),(round(fs/fsCheck1)));
            fs = fsCheck1;
            GW2 = A:1/fs:B;
            GW2 = GW2';
	    fs
	    size(GW1)
	    size(GW2)
            GWc = [GW2 GW1];
            %clearvars GW1 GW2
        end
    end
	
    if segslice
        fss = round(length(GWc(:,1))/(GWc(length(GWc(:,1))-1,1) - GWc(1,1)));
        %GWc = GWc(round(segStart*fss) + 1:round(segEnd*fss) + 1,:);
        if length(GWc(:,1)) >= (round((segEnd*fss)) + 1)
        %if length(GWc(:,1)) >= (round((segEnd*fss)))
            if segStart
		%GWc = GWc(round(segStart*fss):round(segEnd*fss) + 1,:);
                GWc = GWc(round(segStart*fss)+ 1:round(segEnd*fss) + 1,:);
                if (fnct == 0 || fnct == 1)
                    GW_high = GW_high(round(segStart*fss)+ 1:round(segEnd*fss) + 1);
                    GW_low = GW_low(round(segStart*fss)+ 1:round(segEnd*fss) + 1);
                end
            else
		%GWc = GWc(1:round(segEnd*fss) + 1,:);
                GWc = GWc(1:round(segEnd*fss) + 1,:);  % possibly redundant ?
                if (fnct == 0 || fnct == 1)
                    GW_high = GW_high(1:round(segEnd*fss) + 1);
                    GW_low = GW_low(1:round(segEnd*fss) + 1);
                end
            end
        else
            if segStart
		GWc = GWc(round(segStart*fss) + 1:length(GWc(:,1)),:);
                if (fnct == 0 || fnct == 1)
                    GW_high = GW_high(round(segStart*fss) + 1:length(GW_high));
                    GW_low = GW_low(round(segStart*fss) + 1:length(GW_low));
                end
            else
		GWc = GWc(1:length(GWc(:,1)),:);
                if (fnct == 0 || fnct == 1)
                    GW_high = GW_high(1:length(GW_high));
                    GW_low = GW_low(1:length(GW_low));
                end
            end
        end
    end

    GW_1c = GWc(:,2);

    lng = length(GW_1c);

    if ~full_SNR && ~(fnct==6)% process after slicing
        if ~uSnr
		disp('Calculating SNR (post)...')
        	[~,Snr,~,~] = Snr_1(GW_1c(1:lng),fs,0,1); % calculate SNR on raw signal
        	Sr = 10^(Snr/10);
	end
    end

    if uSnr
	Snr = uSnrV;
	Sr = 10^(Snr/10);
    end

    time = lng/fs;

    %Snr = snr(GW_1c);

    if check
        if segslice
            fsCheck = fss;
        else
            fsCheck = floor(length(GWc(:,1))/(GWc(length(GWc(:,1))-1,1) - GWc(1,1)));
        end
        if fs
            disp(strcat('Expected sampling rate:...',num2str(fs)))
            disp(strcat('Calculated sampling rate:...',num2str(fsCheck)))
        else
            disp(strcat('Calculated sampling rate:...',num2str(fsCheck)))
            fs = fsCheck;
        end 

        if fsCheck ~= fs
            disp('Sampling rate mismatch! -> resetting...')
            fs = fsCheck;
        end
    end

    if check

        centSampTime = round(lng/2);
        centTime = GWc(centSampTime + 1,1);

        disp(strcat('Required central time:...',num2str(GPS_central)))
        disp(strcat('Stamped central time:...',num2str(centTime)))

        %if (centTime > (GPS_central + 1/fs)) || (centTime < (GPS_central - 1/fs))
        if centTime ~= GPS_central
            disp('Warning: Central time mismatch! -> Resetting...')
            GPS_central = centTime;
        end

    end

    clearvars GWc

    if wvlt
        GW_1c = waveletbld( GW_1c,15,0.5 );
        time_window = length(GW_1c)/fs;
    end

    if dnsmple == 1
        if fs > fsd
            GW_1c = downsample(GW_1c,(length(GW_1c)/round(time * fsd)));
            if (fnct == 0 || fnct == 1)
                GW_high = downsample(GW_high,(length(GW_high)/round(time * fsd)));
                GW_low = downsample(GW_low,(length(GW_low)/round(time * fsd)));
            end
            fs = fsd;
        end
    end

    if decm == 1
        if fs > fsd
            GW_1c = decimate(GW_1c,(round(fs/fsd)));
            if (fnct == 0 || fnct == 1)
                GW_high = decimate(GW_high,(round(fs/fsd)));
                GW_low = decimate(GW_low,(round(fs/fsd)));
            end
            fs = fsd;
        end
    end

    if slow == 1
        fs = fsd;
    end

    if autom && (time_window ~= -1)
        [offset_multiplier,timeN]=mcalc(time_window,ol,seg,fs);
        time_window = timeN;
        disp(strcat('Extended Time Window is...',num2str(time_window),'...s'))
    elseif autom && (time_window == -1)
        disp('Automation requested, but time_window set to -1.  Multiplier will be set to 1')
        offset_multiplier = 1;
    end

    if time_window ~= -1
        samp_add = round(fs * (time_window/2));

        central_time = round(length(GW_1c)/2);
        lower = round(central_time - samp_add);
        if lower < 1
            lower = 1;
        end
        upper = round(central_time + samp_add);
        if upper > length(GW_1c)
            upper = length(GW_1c);
        end
    else
        lower = 1;
        upper = length(GW_1c);
        time_window = upper/fs;
    end
    %nsamples = upper - lower
    %reTime = nsamples/fs
    %------------------------------------------

    if (fnct == 0 || fnct == 1 || fnct == 4) && filt == 1

        disp('Creating low-pass filter...')

        lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
             'PassbandFrequency',pass_band_fl,'PassbandRipple',0.2, ...
             'SampleRate',fs);

        disp('Creating high-pass filter...')     

        hpFilt = designfilt('highpassiir','FilterOrder',8, ...
             'PassbandFrequency',pass_band_fh,'PassbandRipple',0.1, ...
             'SampleRate',fs);

        disp('Performing low-pass filter...')   
        GW_low = filter(lpFilt,GW_1c(:,:));

        disp('Performing high-pass filter...')
        GW_high = filter(hpFilt,GW_1c(:,:));

    end
%-----------------------------------------------------------------

% --------------------------Calculations--------------------------
end
if fnct==0 % cross bispectrum
    bic_matrix = [GW_low.';GW_high.';GW_1c'];
    bic_matrix = bic_matrix.';

    disp('Performing cross-bispectral calculation...')
    %[bisp,freqBsp] = bispecd_single(tMatrix,0,1,16384,10);
    [bisp,freqBsp] = bispecdx_single(bic_matrix(lower:upper,1),bic_matrix(lower:upper,2),bic_matrix(lower:upper,3),0,1,16384,10,8192,0);
    
    freqBsp_F = freqBsp * fs; % save an unnormalized version of the frequency data
   
    %disp('Saving matrix...')
    %dlmwrite('bicumulant.txt',bisp,'delimiter',' ');
    
    %clearvars -except bisp freqBsp_F plot_img GW_channel fs GPS_central time_window pass_band_fl pass_band_fh
    
    if plot_img == 1 && bicoh == 1
        colormap jet
        imagesc(freqBsp_F,freqBsp_F,abs(bisp))
        GW_channel_n = strrep(GW_channel,'_','\_');
        GW_channel_n = regexprep(GW_channel_n,'\d','');
        GW_channel_n = strsplit(GW_channel_n,'.');
        title_string = strcat('Bicumulant of&',GW_channel_n{1},' (@ f_s =&',num2str(fs),'), GPS:&',num2str(GPS_central),', Dur. =&', ...
        num2str(time_window),' sec. \newline Low band <&',num2str(pass_band_fl),' Hz, High band >&',num2str(pass_band_fh),' Hz');
        title_string = strrep(title_string,'&',' ');
        title(title_string)
        xlabel('f1 (Hz)'), ylabel('f2 (Hz)')
    end
    
    %clearvars -except bisp freqBsp_F
elseif fnct==1  % cross bicoherence
    
    tiled = tile;  % 1 to enable, 0 to disable
    
    seg_length = round((length(GW_1c(lower:upper))*offset_multiplier)/seg);
    bic_matrix = [GW_low.';GW_high.';GW_1c'];
    bic_matrix = bic_matrix.';

    disp('Performing cross-bicoherence calculation...')

    if mode == 0
            disp('"Mode 0" selected.  "Periodic" Hann window will be used.')
	    [bisp,freqBsp,randVec] = bicoherdx_single(bic_matrix(lower:upper,1),bic_matrix(lower:upper,2),bic_matrix(lower:upper,3),0,1,rSequence,0,seg_length,hann(seg_length,'periodic'),seg_length,ol); %**** 200 works
    elseif mode == 1
            disp('Mode 1 (Chebyshev window) selected')
	    [bisp,freqBsp,randVec] = bicoherdx_single(bic_matrix(lower:upper,1),bic_matrix(lower:upper,2),bic_matrix(lower:upper,3),0,1,rSequence,0,seg_length,chebwin(seg_length,200),seg_length,ol); %**** 200 works
    end

    %[bicx,freqBsp] = bicoherx(bic_matrix(lower:upper,1),bic_matrix(lower:upper,2),bic_matrix(lower:upper,3),128,hann(64,'periodic'));
    freqBsp_F = freqBsp * fs; % save an unnormalized version of the frequency data
   
    %disp('Saving matrix...')
    %dlmwrite('bicoherencex.txt',bicx,'delimiter',' ');
    
    if ~tile
        time_window = seg_length/fs;
    end
    
    if plot_img == 1
        colormap jet
        imagesc(freqBsp_F,freqBsp_F,abs(bisp))
        GW_channel_n = strrep(GW_channel,'_','\_');
        GW_channel_n = regexprep(GW_channel_n,'\d','');
        GW_channel_n = strsplit(GW_channel_n,'.');
        title_string = strcat('(Cross) Bicoherence of&',GW_channel_n{1},' (@ f_s =&',num2str(fs),'), GPS:&',num2str(GPS_central),', Dur. =&', ...
        num2str(time_window),' sec.');
        title_string = strrep(title_string,'&',' ');
        title(title_string)
        xlabel('f1 (Hz)'), ylabel('f2 (Hz)')
    end
    
    T = bspcThresh(seg,seg_length,GW_1c,randVec,fs,0,Snr);
    disp(strcat('T =...',T));
    disp(strcat('SNR =...',num2str(Sr)));

    if threshMan
        Tthresh = HbicohThresh;
    else
        Tthresh = str2double(T);
    end

    final_results = bcohCrit( Tthresh , bisp, freqBsp_F, str2double(T) );

    clearvars -except bisp freqBsp_F randVec T Snr Sr final_results GPS_central path GPS_central ifoCh fnct HbicohThresh threshMan ChannelA ChannelB

    save(strcat(path,ifoCh,'_',num2str(GPS_central),'.mat'), '-v7.3') % save output
 
elseif fnct==2 % auto bicoherence
    
    tiled = tile;  % 1 to enable, 0 to disable
    
    if tiled
        % preallocate time vector matrix
        seg_length = length(GW_1c(lower:upper,1));
        
%         if rem(seg_length,2)
%             upper = upper - 1;
%             seg_length = length(GW_1c(lower:upper,1));
%         end
        %GW_1c = GW_1c(lower:upper,1);
        tMatrix = zeros(seg_length,seg);
        for i = 1:seg
        tMatrix(:,i) = GW_1c(lower:upper,1);
        % fill the matrix.  This matrix will be passed to bispecd() which will
        % phase randomize each instance during calculation
        end
    else
        seg_length = round((length(GW_1c(lower:upper))*offset_multiplier)/seg);
    end
    
    disp('Performing bicoherence calculation...')
    
    if tiled
        disp('Tiled option chosen')
        
        % mode = 1;
        
        if mode == 0
            disp('Mode 0 (Hann window) selected')
            %[bisp,freqBsp,randVec] = bispecd_single(tMatrix,0,1,tiled,seg_length,10);
            [bisp,freqBsp,randVec] = bicoherd_single_1(tMatrix,0,1,rSequence,tiled,seg_length,hann(seg_length,'periodic'),seg_length,ol);
            %[d,e,f] = bicoherd_single_follow(tMatrix,results,0,randomize,rSequence,tiled,seg_length,hann(seg_length,'periodic'),seg_length,ol);
        elseif mode == 1
            disp('Mode 1 (Chebyshev window) selected')
            [bisp,freqBsp,randVec] = bicoherd_single_1(tMatrix,0,1,rSequence,tiled,seg_length,chebwin(seg_length,200),seg_length,ol); %**** 200 works
        elseif mode == 2
            disp('Mode 2 selected')
            %[bisp,freqBsp,randVec] = bispecd_single_div(tMatrix,1,1,seg_length,10,2);
            [bisp,freqBsp,randVec] = bicoherd_single_1(tMatrix,0,1,rSequence,tiled,seg_length);
        elseif mode == 3
            disp('Mode 3 selected')
            % *FFT will be zero-padded if nfft exceeds segment length
            %[bisp,freqBsp,randVec] = bispecd_single(tMatrix,0,1,tiled,16384,10);
            [bisp,freqBsp,randVec] = bicoherd_single_1(tMatrix,0,1,rSequence,tiled,16384);
        end
    else
        disp('No-tile option selected')
        disp(strcat('Segments are...',num2str(seg_length/fs),'...seconds'))
        %fres = 0.5;
        %nft = round(fs * (1/fres));
        disp(strcat('Frequency Resolution:...',num2str(fs/seg_length),'...Hz'))
        %disp(strcat('Frequency Resolution:...',num2str(fres),'...Hz'))
        %[bisp,freqBsp,randVec] = bispecd_single(GW_1c(lower:upper),0,1,tiled,16384,10,seg_length,0);  % notile
        %ol = 50; % 0 for no overlap
        %[bisp,freqBsp,randVec] = bispecd_single(GW_1c(lower:upper),0,1,tiled,nft,10,seg_length,ol);
        %mode = 0;
        if mode == 0
            disp('"Mode 0" selected.  "Periodic" Hann window will be used.')
            [bisp,freqBsp,randVec] = bicoherd_single_1(GW_1c(lower:upper),0,1,rSequence,tiled,seg_length,hann(seg_length,'periodic'),seg_length,ol);
        elseif mode == 1
            disp('Mode 1 (Chebyshev window) selected')
            [bisp,freqBsp,randVec] = bicoherd_single_1(GW_1c(lower:upper),0,1,rSequence,tiled,seg_length,chebwin(seg_length,200),seg_length,ol); %**** 200 works
        end
    end
    
%     if (upper - lower) >= fs
%         samp_seg = round((upper - lower)/30);
%         
%         disp(strcat('Sample segment length:',num2str(upper - lower),'/30 = ',num2str(samp_seg)))
%     else
%         samp_seg = round((upper - lower)/8);
%         disp(strcat('Sample segment length: ',num2str(upper - lower),'/8 = ',num2str(samp_seg)))
%     end
    
    %[bisp,freqBsp] = bicoher(bic_vec(lower:upper,1),0,128,40,samp_seg);
    
    %[bisp,freqBsp,randVec] = bicoherd_single_1(tMatrix,0,1,tiled,16384,0);
    freqBsp_F = freqBsp * fs; % save an unnormalized version of the frequency data
    
    %disp('Saving matrix...')
    %dlmwrite('bicoherence.txt',bisp,'delimiter',' ');
    
    clearvars -except bisp freqBsp_F plot_img GW_channel fs GPS_central time_window pass_band_fl pass_band_fh randVec seg seg_length GW_1c tile Snr Sr path GPS_central ifoCh fnct threshMan HbicohThresh
    
    if ~tile
        time_window = seg_length/fs;
    end
    
    if plot_img == 1
        colormap jet
        imagesc(freqBsp_F,freqBsp_F,abs(bisp))
        GW_channel_n = strrep(GW_channel,'_','\_');
        GW_channel_n = regexprep(GW_channel_n,'\d','');
        GW_channel_n = strsplit(GW_channel_n,'.');
        title_string = strcat('(Auto) Bicoherence of&',GW_channel_n{1},' (@ f_s =&',num2str(fs),'), GPS:&',num2str(GPS_central),', Dur. =&', ...
        num2str(time_window),' sec.');
        title_string = strrep(title_string,'&',' ');
        title(title_string)
        xlabel('f1 (Hz)'), ylabel('f2 (Hz)')
    end
    
    %Bispectrum of OAF-CAL\_DARM\_DQ (@ f_s = 2048), GPS: 1102757206
    %Low band < 20 Hz, High band > 20 Hz 
    T = bspcThresh(seg,seg_length,GW_1c,randVec,fs,0,Snr);
    disp(strcat('T =...',T));
    disp(strcat('SNR =...',num2str(Sr)));
   
    if threshMan
    	Tthresh = HbicohThresh;
    else
    	Tthresh = str2double(T);
    end
 
    final_results = bcohCrit( Tthresh , bisp, freqBsp_F, str2double(T) );
    %T = bspcThresh(seg,seg_length,GW_1c(lower:upper,1),randVec);
    
    %disp(strcat('T value:...',num2str(T)))
    
    clearvars -except bisp freqBsp_F randVec T Snr Sr final_results GPS_central path GPS_central ifoCh fnct HbicohThresh threshMan

    save(strcat(path,ifoCh,'_',num2str(GPS_central),'.mat'), '-v7.3') % save output

elseif fnct==3 % auto bispectrum
     % mode = 1;
%     if (upper - lower) >= fs
%         samp_seg = round((upper - lower)/64);
%         
%         disp(strcat('Sample segment length:',num2str(upper - lower),'/64 = ',num2str(samp_seg)))
%     else
%         samp_seg = round((upper - lower)/4);
%         disp(strcat('Sample segment length: ',num2str(upper - lower),'/4 = ',num2str(samp_seg)))
%     end

    %Snr = snr(GW_1c(lower:upper));
    %[~,Snr,~,~] = Snr_1(GW_1c(lower:upper),fs,0,1);
    %Sr = 10^(Snr/10);
    
    tiled = tile;  % 1 to enable, 0 to disable
    
    if tiled
        % preallocate time vector matrix
        seg_length = length(GW_1c(lower:upper,1));
        if check && (mode==1)
            if mod(log2(seg_length),1)
                disp('WARNING: Your segment length (and hence fft length) is not a power of 2!')
            end
        end
                
%         if rem(seg_length,2)
%             upper = upper - 1;
%             seg_length = length(GW_1c(lower:upper,1));
%         end
        %GW_1c = GW_1c(lower:upper,1);
        tMatrix = zeros(seg_length,seg);
%     else
%         disp('Creating unique segment matrix')
%         seg_length = length(GW_1c)/seg;
%         % preallocate time vector matrix
%         tMatrix = zeros(seg_length,seg);
%         locseg = [1:seg_length]';
        for i = 1:seg
        %if tiled
            tMatrix(:,i) = GW_1c(lower:upper,1);
%         else
%             disp(strcat('Creating segment...',num2str(i),'...of...',num2str(seg)))
%             %locseg(end)
%             %length(GW_1c)
%             tMatrix(:,i) = GW_1c(locseg);
%             locseg = locseg + seg_length;
%         end
        % fill the matrix.  This matrix will be passed to bispecd() which will
        % phase randomize each instance during calculation
        end
    else
        seg_length = round((length(GW_1c(lower:upper))*offset_multiplier)/seg);
        if check
            if mod(log2(seg_length),1)
                disp('WARNING: Your segment length (and hence fft length) is not a power of 2!')
            end
        end
    end
    
    disp('Performing bispectral calculation...')
    %[bisp,freqBsp] = bispecd(bic_vec(lower:upper,1),0,128,hann(64,'periodic'),samp_seg);
    %[bisp,freqBsp] = bispecd(GW_1c(lower:upper,1),0,256,40);
    %[bisp,freqBsp] = bispecd(GW_1c(lower:upper,1),0,randomize,1024,3,samp_seg,75);
    
    if tiled
        disp('Tiled option chosen')
        disp(strcat('Segments are:...',num2str(seg_length),'...samples'))
        fres = fs/seg_length;
        disp(strcat('Frequency Resolution:...',num2str(fres),'...Hz'))
        
        if mode == 0 % zero-padded if time series is less than nfft
            [bisp,freqBsp,randVec,DVnfft] = bispecd_single(tMatrix,0,randomize,rSequence,tiled,16384,10);
        elseif mode == 1  % no zero padding  (normal)
            [bisp,freqBsp,randVec,DVnfft] = bispecd_single(tMatrix,0,randomize,rSequence,tiled,seg_length,10);
            %bispecd_single(tMatrix,0,randomize,rSequence,tilevec,nsamples,10,nsamples,ol);
        elseif mode == 2
            [bisp,freqBsp,randVec,DVnfft] = bispecd_single_div(tMatrix,1,1,seg_length,10,2);
        end
    else
        disp('No-tile option selected')
        disp(strcat('Segments are...',num2str(seg_length/fs),'...seconds'))
        disp(strcat('Segments are:...',num2str(seg_length),'...samples'))
        fres = fs/seg_length;
        %fres = 0.5;
        %nft = round(fs * (1/fres));
        nft = seg_length;
        %disp(strcat('Minimum frequency:...',num2str(1/(seg_length/fs)),'...Hz'))
        disp(strcat('Frequency Resolution:...',num2str(fres),'...Hz'))
        %[bisp,freqBsp,randVec] = bispecd_single(GW_1c(lower:upper),0,1,tiled,16384,10,seg_length,0);  % notile
        %ol = 50; % 0 for no overlap
        [bisp,freqBsp,randVec,DVnfft] = bispecd_single(GW_1c(lower:upper),0,1,tiled,nft,10,seg_length,ol);
    end
    %[bisp,freqBsp] = bispecd(tMatrix,0,1,16384,10); %working settings; use either 5 or 10 for window height
    % ################
    freqBsp_F = freqBsp * fs; % save an unnormalized version of the frequency data
    %##########
    %[~,freqBsp_F,bisp] = bicoherence_N1(tMatrix,64,fs,1024) %working settings; use either 5 or 10 for window height
    
    %[bisp,freqBsp] = bispecd(bic_matrix(lower:upper,:)) % coarse version
    
    %disp('Saving matrix...')
    %dlmwrite('bispectrum.txt',bisp,'delimiter',' ');
    %DVnfft = seg_length;
    
    %clearvars -except Snr DVnfft bisp freqBsp_F plot_img GW_channel fs GPS_central time_window pass_band_fl ...
    %    pass_band_fh randVec GW_1c seg tile seg_length tMatrix tiled
    
    if ~tile
        time_window = seg_length/fs;
    end
    
    if plot_img == 1
        figure;
        colormap jet
        imagesc(freqBsp_F,freqBsp_F,abs(bisp))
        GW_channel_n = strrep(GW_channel,'_','\_');
        GW_channel_n = regexprep(GW_channel_n,'\d','');
        GW_channel_n = strsplit(GW_channel_n,'.');
        title_string = strcat('Bispectrum of&',GW_channel_n{1},' (@ f_s =&',num2str(fs),'), GPS:&',num2str(GPS_central),', Dur. =&', ...
        num2str(time_window),' sec.');
        title_string = strrep(title_string,'&',' ');
        title(title_string)
        xlabel('f1 (Hz)'), ylabel('f2 (Hz)')
    end
   
    %results = peak_find(bisp,freqBsp_F,1);
    
    %T = bspcThresh(seg,seg_length,GW_1c(lower:upper,1),randVec);
    
    %disp(strcat('T value:...',num2str(T)))
    T = bspcThresh(seg,seg_length,GW_1c,randVec,fs,0,Snr);
    disp(strcat('T =...',T));
    disp(strcat('SNR =...',num2str(Sr)));
    
    resultsC = bpCrit(seg,bisp,freqBsp_F,T);
    % results = bpCrit(seg,a,b*fs,T);
    goflag = 0;
    results = [resultsC{1};resultsC{2};resultsC{3};resultsC{4}];
    if ~isempty(results)
        goflag = 1;
    end
    
    %results = bpCritTheor(seg, seg_length,Snr, bisp, freqBsp_F); %******

%     [d,e,f] = bicoherd_single_follow(tMatrix,results,0,1,tiled,nft,hann(seg_length,'periodic'),seg_length,ol); %****
    if tiled && goflag
        %ol = 0; % 0 for no overlap
        if mode == 0
            [d,e,f] = bicoherd_single_follow(tMatrix,results,0,randomize,rSequence,tiled,16384,hann(seg_length,'periodic'),seg_length,ol); %****
        elseif mode == 1
            [d,e,f] = bicoherd_single_follow(tMatrix,results,0,randomize,rSequence,tiled,seg_length,hann(seg_length,'periodic'),seg_length,ol); %****
        end
        % [Bspec,waxis,randVec] = bicoherd_single_follow (y,  res, plot, randomize, rvector, tile, nfft, wind, nsamp, overlap)
    end
    
    %T = bspcThresh(seg,seg_length,GW_1c(lower:upper),f);
    %T = bspcThresh(seg,seg_length,GW_1c,f);
    if goflag
        %disp(strcat('Threshold for (strong) bicoherence:...',num2str(T)))
        disp(strcat('Threshold for (strong) bicoherence:...',T))
        %final_results = bcohCrit( T, d, freqBsp_F );
	if threshMan
		Tthresh = HbicohBspcThresh;
	else
		Tthresh = str2double(T);
	end
        final_results = bcohCrit( Tthresh, d, freqBsp_F, str2double(T) );
        %results = 0;
    end
    clearvars -except bisp freqBsp_F fs randVec resultsC final_results d Snr Sr path GPS_central ifoCh fnct HbicohBspcThresh threshMan
    save(strcat(path,ifoCh,'_',num2str(GPS_central),'.mat'), '-v7.3') % save output 
elseif fnct==4 % magnitude squared coherence between upper and lower bands
    disp('Performing coherence analysis between high and low bands...')
    [cxy,f] = mscohere(GW_low,GW_high,[],[],16384,fs);
    cohere_matrix = [f cxy];
    save(strcat(path,ifoCh,'_',num2str(GPS_central),'.mat'), '-v7.3') % save output
elseif fnct==5 % autocoherence
    [acoh,freq] = autocoherence(GW_1c(lower:upper), length(GW_1c(lower:upper)));
    figure; plot(fs*freq,abs(acoh))
    save(strcat(path,ifoCh,'_',num2str(GPS_central),'.mat'), '-v7.3') % save output
elseif fnct==6 % stft based "bispectrogram"
    font_size = 20;
    %xlimit = [0 2000];
    %ylimit = [0 2000];
    
    pkVal = 0;

    xlimit = [0 round(fs/2)];%*
    ylimit = [0 round(fs/2)];%*
    
%   xlimit = [0 round(fs/16)];
%   ylimit = [0 round(fs/16)];
    
    xl = 'f3 (Hz)';
    yl = 'f2 (Hz)';
    
    bspc = 1; % 1 for bispectrum, 0 for biphase
    bicoh = 0; % 1 for bicoherence, 0 for bispectrum/biphase 
    
    sve = 0; % set whether to save each matrix or not
    seg_size = 2048;  %1024 works
    %seg_size = 1024;
    [spec,freq,tme] = stft(GW_1c(lower:upper,1),seg_size,50,seg_size,fs);
    
    if bspc && amax
        [upr_sp, aindx] = max(abs(spec(:)));
        %upr_sp = max(abs(spec(:)));
        %upr = (upr_sp)^(2);
        upr = (spec(aindx)*conj(spec(aindx)));
    end

    if ~amax
	disp(strcat('Manual maximum set:...',num2str(upr)))
    end

    %upr = abs(spec(aindx)*spec(aindx)*conj(spec(aindx)));
    %figure; colormap jet; imagesc(c,b,ln(abs(a)))
    [rw,cl] = size(spec);
    %flag = 1;6
    
    fctrs = factor(cl);
    if length(fctrs) == 1
        disp('Only one factor: Adjusting...')
        cl = cl - 1;
        fctrs = factor(cl);
        disp(strcat('New total sample length in time:...',num2str(cl),'...s'))
    else
        disp(strcat('Total sample length in time:...',num2str(cl),'...s'))
    end
    disp(strcat('Factors:...', num2str(fctrs)))
    %slices = 0;
    slices    = input('Select number of slices       ---> [max(Factors)] ');
    
    if (isempty(slices)) 
        slices = max(fctrs); 
    end
    
    if rem(cl,slices)
        disp(strcat('Warning: Not a factor of -> ',num2str(cl),'!  Using max(Factors)'))
        slices = max(fctrs);
    end
    
    if slices <=0
        disp('Warning: Number of slices chosen <= 0! Using max(Factors)')
        slices = max(fctrs);
    end
    
    figure; 
    %slices = min(fctrs);
    %function [Bspec,randVec] = bispecd_spectro (ft, nrecs ,plot, randomize, rvector, tile, wind, nsamp)
    %[bisp,freqBsp,randVec,DVnfft] = bispecd_single(tMatrix,0,randomize,rSequence,tiled,seg_length,10);
    for slice = 1:slices
        disp(strcat('Calculating slice:...',num2str(slice),'...of...',num2str(slices)))
        if ~bicoh
            [bisp,~] = bispecd_spectro(spec(:,(slice-1)*(cl/slices) + 1), seg ,0, randomize, rSequence, 1, 10, 1024);
        else
            [bisp,~] = bicoherd_spectro(spec(:,(slice-1)*(cl/slices) + 1), seg ,0, randomize, rSequence, 1, 10, 1024);
        end
        if bspc
            idx = find(abs(bisp)>upr);
            if idx
                bisp(idx) = upr;
            else
                bisp(1,1) = upr;
            end
        end
%         if flag == 1
%             upr = max(abs(bisp));
%             flag = 0;
%         end
        if sve
            bispectro{slice} = bisp;
        end

	if max(abs(bisp(:))) > pkVal
		pkVal = max(abs(bisp(:)));
	end

        if plot_img
            
            if sve
                colormap jet; imagesc(freq,freq,abs(bispectro{1,slice}))
                title(strcat(num2str(tme((slice-1)*(cl/slices) + 1)),'...s'))
                print(strcat(num2str(tme((slice-1)*(cl/slices) + 1)),'.png'),'-dpng')
            else
                if bspc && ~bicoh
                    colormap jet; imagesc(freq,freq,abs(bisp),[0 upr])
                elseif ~bspc && ~bicoh
                    colormap jet; imagesc(freq,freq,angle(bisp))
                else
                    colormap jet; imagesc(freq,freq,abs(bisp),[0 1])
                end
                %colormap jet; imagesc(freq,freq,abs(bisp))
                colorbar
                
                if bspc && ~bicoh
                    title(strcat('Bispectrum: ',num2str(tme((slice-1)*(cl/slices) + 1)),'...s'))
                elseif ~bspc && ~bicoh
                    title(strcat('Biphase: ',num2str(tme((slice-1)*(cl/slices) + 1)),'...s'))
                else
                    title(strcat('Bicoherence: ',num2str(tme((slice-1)*(cl/slices) + 1)),'...s'))
                end
                set(gca,'YDir','normal');
                set(gca,'FontSize',font_size);
                set(gca,'xLim',xlimit);
                set(gca,'yLim',ylimit);
                xlabel(xl), ylabel(yl)
                %print(strcat(num2str(tme((slice-1)*(cl/slices) + 1)),'.png'),'-dpng')
                print(strcat('img_',num2str(slice),'.png'),'-dpng')
            end
        end
    end
elseif fnct==7
    GW_1 = load(GW_channel);
    GW_2 = load(GW_channel_1);
    
    autoalign = 1; % attempt auto align of central times
    
    t_window = 25; % seconds
    num_shift = 20; % number of time-shifts
    lmin = -0.03; %s
    lmax = 0.03; %s
   
    %xlimit
    ylimit = [0.65 1];
    xl = 'Frequency (Hz)';
    yl = 'Coherence';
    
    results = {};
    inc = (lmax - lmin)/num_shift;
    inc_ind = round(inc*fs);
    lmin_ind = floor(lmin * fs);
    lmax_ind = ceil(lmax * fs);
    
    mid1 = ceil(length(GW_1(:,2))/2);
    mid2 = ceil(length(GW_2(:,2))/2);
    
    % sanity check
    tm1 = GW_1(mid1,1);
    tm2 = GW_2(mid2,1);
    
    disp(tm1)
    disp(tm2)
    
    if tm1 ~= tm2
        disp('Warning: times are not aligned!!')
        if autoalign
            disp('Attempting auto-align...')
            mid2_offset = tm1 - tm2; % sec
            %mid2_offset = 1; % sec
            mid2 = mid2 + round(mid2_offset * fs);
            tm1 = GW_1(mid1,1);
            tm2 = GW_2(mid2,1);
            if tm1 ~= tm2
                disp('Warning: auto-align failed!')
            else
                disp('Auto-align successful...')
            end
        end
    end
    % ---------------------
    
    samp_upper = floor((t_window * fs)/2);
    samp_lower = ceil((t_window * fs)/2);
    
    up_index1 = mid1 + samp_upper;
    lower_index1 = mid1 - samp_lower;
    
    up_index2 = mid2 + samp_upper;
    lower_index2 = mid2 - samp_lower;
    
    ui2 = up_index2;
    li2 = lower_index2;
    
    %(up_index1 - lower_index1)/fs
    
    ts1 = GW_1(lower_index1:up_index1,2);
    
    ntests = ceil((lmax_ind - lmin_ind)/inc_ind);
    
    j = 1;
    
    tmes = zeros(ntests,1);
    
    for i = lmin_ind:inc_ind:lmax_ind
        
        up_index2 = up_index2 + i;
        lower_index2 = lower_index2 + i;
    
        ts2 = GW_2(lower_index2:up_index2,2);
        
        disp(strcat('Performing test...',num2str(j),'...of...',num2str(ntests)))
        %disp (samp_upper - samp_lower)
        [Coher,F] = mscohere(ts1,ts2,[],[],samp_upper + samp_lower,fs);
        
        results{j} = [sqrt(Coher) F];
        tmes(j) = i;
        
        j = j + 1;
        up_index2 = ui2;
        lower_index2 = li2;
    end
    
    trm = 8;
    c = linspace(1,results{1,1}(round(length(results{1,1}(:,2))/trm),2),round(length(results{1,1}(:,2))/trm));
    figure;
    colormap jet;
                
    %cols = 3;
    
%     if rem(length(results),cols)
%         row_sbplt = length(results) + 1;
%     else
%         row_sbplt = length(results);
%     end
    
    for i=1:length(results)
        %subplot(row_sbplt/2,cols,i);
        disp(strcat('Plotting...',num2str(i),'...of...',num2str(ntests)))
        scatter(results{1,i}(1:round(length(results{1,i}(:,2))/trm),2),results{1,i}(1:round(length(results{1,i}(:,2))/trm),1),1.5,c)
        title(strcat(num2str(tmes(i)/fs),'...s'))
        %set(gca,'xLim',xlimit);
        set(gca,'yLim',ylimit);
        xlabel(xl), ylabel(yl)
        print(strcat(num2str(i),'.png'),'-dpng')
    end
end

disp('Done!')
