sequenceN = 'seq_552015_2.txt'; % Acceptable
%sequenceN = 'seq_552015_6.txt';
rSequence = load(sequenceN);

%fs = 1024; %higher sample rate gives greater resolution, but also increases processing time
%fs = 2048;
fs = 256;
window = 0;  % bicoherence: 0 for Hann window, 1 for Chebyshev
%fs = 4096;

%fs = 16384;
%time = 30.0;
%time = 120.0;
autom = 1; % 1 for autoset, 0 for manual set (below)
calTrim = 0; % trim time series to calculation length
%time = 30;
%----------------------------------manual set
if ~autom % use with manual set
    time = 30; % full time (for SNR)
    timeN = 4.12; % calculation time
    multiplier = 29.46;
else % use with auto set
    time = 30;
    timeC = 4; % desired segment time
end
%----------------------------------

%time = 64;
rebuild = 1;
%N = fs * time;
%multiplier = 1;
%ol = 99.9;  % normally 50(%); max 99.9
ol = 98;

tilevec = 0;
randomize = 1;

seg = 30;  % number of time segments (experiments) to use in analysis

if autom
    [multiplier,timeN]=mcalc(timeC,ol,seg,fs);
    %time = timeN;
else
end

%segT = seg * (100/ol) * (1/multiplier);
%seg = 60;
%amplititudes --------------
% For QPC, it seems that we require 
a1 = 1.0;

a2 = 0.8;

%a2 = 0.5; 
%a2 = 2;
a3 = 0;
%a4 = 0;
a5 = 0.05; %noise amplititude; 0.05 -> SNR ~ 34.4
%a5 = 0;
%a5 = 0.3;
%a5 = 0.001;
%a5 = 0;
%a5 = 0;
%a6 = 0;

% frequencies ----------------------------

% f1 = 1000.0;
% f2 = 4000.0;

% bad choice, no QPC -------------------------------
%f1 = 1000.0;  % f1 + f2 = 2600 > 0.5 fs
%f2 = 1600.0; 

%f1 = 1000.0; % better, but still weak coupling
%f2 = 1500.0;
% ------------------------------------------

%f1 = 517;
%f2 = 1023;

%f1 = 516;
%f2 = 1021;

%working----------------------
%f1 = 513;
%f2 = 4;

%f1 = 400;
%f2 = 1600;

%f1 = 100;
%f2 = 400;

%f1 = 47; % ~ 517/1023
%f2 = 93; % ~ 517/1023

% f1 = 51.2;
% f2 = 102.4;

%f1 = 51;
% f1 = 61;
% f2 = 102;

%f1 = 512.5;
%f2 = 0.75;

%f1= 60;
%f2 = 46;
%f2 = 53;
%f1 = 512.25;
%f2 = 4;

% f1 = 515.25;
% f2 = 4;
%f2 = 16;

%f1 = 513;
%f2 = 0.125;

f1 = 4;
f2 = 60;

% f1 = 515.3;
% f2 = 509.7;

% f1 = 515.25;
% f2 = 510;
%f2 = 0.25;

%f1 = 300;
%f2 = 724;

% !! f1 + f2 = f3 <= fs/2 !!   <---required condition
% -----------------------------------------

if rebuild  % option to disable if this variable (ts) already exists in memory
    
    %x=randn(1,N+1);  %gaussian vector;
    if a5 ~= 0
        xs = 1;
    else
        xs = 0;
    end
    %ts = zeros([N,1]);
    % 
    % phi1 = pi * [(-1)*fliplr(rand(1,0.5*N)) rand(1,0.5*N)]';
    % phi2 = pi * [(-1)*fliplr(rand(1,0.5*N)) rand(1,0.5*N)]';
    % phi3a = pi * [(-1)*fliplr(rand(1,0.5*N)) rand(1,0.5*N)]';
    % 
 
    %ts = sig_gen1(time,fs,f1,f2,a1,a2,a5,xs);
    ts = sig_gen(time,fs, f1, f2, a1, a2,a5,xs);
    %ts1 = sig_gen(time,fs, 4, 56, a1, a2,a5,xs);
    %ts = ts + ts1;
end

if autom || calTrim
    smcnt = ceil(fs * timeN);
    ts1 = ts;
    ts = ts(1:smcnt);
    time = smcnt/fs;
else
    ts1 = ts;
end

[~,S,~,~] = Snr_1(ts1,fs,0,1);
Sr = 10^(S/10);

if tilevec
    % preallocate time vector matrix
    nsamples = length(ts);
    if rem(nsamples,2)
       ts = ts(1:end-1);
       nsamples = length(ts);
    end
    
    tMatrix = zeros(length(ts),seg);
    segT = seg;
    %tiled = 0;

    for i = 1:seg

            tMatrix(:,i) = ts;

       % fill the matrix.  This matrix will be passed to bispecd() which will
       % phase randomize each instance during calculation
    end
else
        nsamples = floor((time/seg)* fs);
        if rem(nsamples,2)
            %ts = ts(1:end-1);
            %nsamples = length(ts);
            nsamples = nsamples - 1;
        end
        
        nsamples = round(nsamples * multiplier);
        
        %********************************************
        
        olm  = fix(nsamples * ol / 100);
        nadvance = nsamples - olm;
        segT    = fix ( ((length(ts)) - olm) / nadvance);
        disp(strcat('Calculated number of segments:...',num2str(segT)))
        
        %**********************************************
        
        tMatrix = ts;
end
%[a,b] = bispecd(ts,0,1024,64,seg,50); % Low fs (?)
%[a,b] = bispecd(ts,0,0,1024,5,32,50); % works well, higher fs; greater
%[a,b] = bispecd(ts,0,0,1024,3,128,75); % gives 36 segments at fs = 1024
%[a,b] = bispecd(ts,0,0,1024,5,72,75); % works well, higher fs; greater  %updated, still seems to work

%number of segments significantly improves performance: 30 suggested; 
%increasing nfft gives better resolution; a value of 5 for the
%window seems to be the best value--too high causes distortion; 50% overlap
%seems, as well, to be best

%[a,b] = bispecd(ts,0,1024,3,90,75);  % best settings, phase randomization, fs = 1024
%[a,b] = bispecd(ts,0,1,1024,3,72,75); %working, with phase randomization
%[a,b] = bispecd(ts,0,1,1024,3,72,0); %experimental, with phase randomization
% 
%[~,b,a] = bicoherence_N1(ts,64,fs,1024)
%[~,b,a] = bicoherence_N1(tMatrix,64,fs,1024)

% --------------------------------------------------------------
% perform bispectral calculation
disp('Performing bispectral analysis...')

%[a,b] = bispecd(tMatrix,0,1,1024,10); %working settings; use either 5 or 10 for; window height

%[a,b] = bispecd_single(tMatrix,0,1,1,16384,10); %working settings; use either 5 or 10 for window height
%[a,b] = bispecd_single(tMatrix,0,1,1,nsamples,10);

[a,b,c] = bispecd_single(tMatrix,0,randomize,rSequence,tilevec,nsamples,10,nsamples,ol); %*******

%S = snr(ts);

if randomize
    T = bspcThresh(segT,nsamples,ts1,c,fs,0,Sr);
else
    T = '0.9999900000000000';
end
%results = bpCritTheor(seg,nsamples,S,a,b*fs);

disp(strcat('T =...',T));
disp(strcat('SNR =...',num2str(Sr)));
%disp(strcat('T = ',num2str(T)));
results = bpCrit(segT,a,b*fs,T);

figure;
colormap jet;
imagesc(b*fs,b*fs,abs(a));
colorbar
set(gca,'YDir','normal');
%bicoherd_single_follow (y,  res,
%[d,e,f] = bicoherd_single_follow(tMatrix,results,0,1,tilevec,4*fs,hann(nsamples,'periodic'),nsamples,ol); %****
%[a,b,c] = bispecd_single_div(tMatrix,1,1,nsamples,10,4);

%T = bspcThresh(seg,nsamples,ts,f);

%final_results = bcohCrit( 0.8, d, b * fs,str2double(T) );
%[~,b,a] = bicoherence_N1(tMatrix,64,fs,16384);

%[a,b,c] = bispecd_single(tMatrix,0,randomize,rSequence,tilevec,nsamples,10,nsamples,ol);
if window == 0
    [d,e,~] = bicoherd_single_1(tMatrix,0,randomize,rSequence,tilevec,nsamples,hann(nsamples,'periodic'),nsamples,ol); %****
else
    [d,e,~] = bicoherd_single_1(tMatrix,0,randomize,rSequence,tilevec,nsamples,chebwin(nsamples,200),nsamples,ol); %**** 200 works
end
% y,  plot, randomize, tile, nfft, wind, nsamp, overlap

%[a,b] = bicoher(tMatrix,0,1,4096,hann(nsamp,'periodic'),nsamp,0);% 3 second chunks
 % make the plots -----------------------------
%  
figure; %bicoherence magnitude plot
colormap jet
% 
imagesc(e*fs,e*fs,abs(d));
colorbar
set(gca,'YDir','normal');
 %imagesc(b*fs,b*fs,abs(a));
%  %imagesc(b*fs,b*fs,a); %HOSA bicoherence
%  
%  imagesc(b,b,abs(a));

final_results = bcohCrit(str2double(T),d,b*fs,str2double(T));

%  figure; %biphase plot
%  colormap jet
%  imagesc(b*fs,b*fs,angle(a));
%  %imagesc(b,b,angle(a));

disp('Done!')
% ------------------------------------------------------------------

%imagesc(b,b,angle(a));
%figure;
%colormap jet
%imagesc(b,b,abs(a));

%[a,b] = bicoher(ts,0,1,1024,3,72,75);
%  [a,b,~] = bicoherence_N1(ts,64,fs,1024)
%   figure;
%   colormap jet
%   imagesc(b,b,abs(a));
%   disp('Done!')
