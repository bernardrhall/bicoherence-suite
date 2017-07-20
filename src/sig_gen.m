function [ts] = sig_gen (time,  fs, f1, f2, a1, a2, a3,a5,v)
%function [ts] = sig_gen (time,  fs, f1, f2, a1, a2,a5,v)
% fs = 1024; %higher sample rate gives greater resolution, but also increases processing time
% time = 5.0;
% f1 = 400.0; 
% f2 = 100.0;
%a3 = 0.1;
% a1 = 1;
% a2 = 1;
% a3 = 1;
% a4 = 0;
% a5 = 0.001; %noise amplititude
% a6 = 0;
if ~f1
    f1 = 1 + ((fs/2)-1)*0.2*rand;
end

rng('shuffle')  %reseed based on time

inc = 1/fs;

t=0:inc:time;

phase_mode = 4; % 0 for normal, 1 for phase-modulated, 2 for linear
% f3 = 102.4;
% f4 = 204.8;
% f6a = 450;
%f3a = f1 + f2;
%hfs = fs/2;

% f1a = hfs*rand; 
% f2a = hfs*rand;
%f4a = hfs*rand;

if phase_mode == 0
    phi1a = -pi; 
    phi2a = -pi;
    phi3a = -pi;

    while phi1a == -pi
        phi1a = -pi + (pi+pi)*rand;
    end
    while phi2a == -pi
        phi2a = -pi + (pi+pi)*rand;
    end
    while phi3a == -pi
        phi3a = -pi + (pi+pi)*rand;
    end
    
    disp('Calculating time series...')
    % QPC
%     ts = a1*cos(2*pi*f1*t + phi1a) + a2*cos(2*pi*f2*t + phi2a) + a3*cos((2*pi*(f1 + f2)*t)+(phi1a+phi2a)) + ...
%     a3*cos((2*pi*f4a*t)+phi3a) + a5*v;
elseif phase_mode == 1
    phi1a = -pi; 

    while phi1a == -pi
        phi1a = -pi + (pi+pi)*rand;
    end
elseif phase_mode == 2
    phi1a = -pi; 
    phi2a = -pi;

    while phi1a == -pi
        phi1a = -pi + (pi+pi)*rand;
    end
    while phi2a == -pi
        phi2a = -pi + (pi+pi)*rand;
    end
elseif phase_mode == 3
    
    phi1a = -pi; 
    phi2a = -pi;
    phi3a = -pi;
    phi4a = -pi; 
    phi5a = -pi;
    phi6a = -pi;

    while (phi1a == -pi) && (phi1a == phi4a)
        phi1a = -pi + (pi+pi)*rand;
        phi4a = -pi + (pi+pi)*rand;
    end
    while (phi2a == -pi) && (phi2a == phi5a)
        phi2a = -pi + (pi+pi)*rand;
        phi5a = -pi + (pi+pi)*rand;
    end
    while (phi3a == -pi) && (phi3a == phi6a)
        phi3a = -pi + (pi+pi)*rand;
        phi6a = -pi + (pi+pi)*rand;
    end
elseif phase_mode == 4
    %phi1a = 0.7; 
    %phi2a = 0.3;
    %phi3a = 0.53;

    phi1a = 0;
    %phi2a = pi/2;
    phi2a = 0;
    phi3a = 0;
end

if v == 0
    %v = zeros(length(ts1),1)';
    %N = floor(fs * time);
    N = length(t);
    v = zeros(N,1)';
elseif v == 1
    %N = floor(fs * time);
    N = length(t);
    v=randn(1,N);
end

disp('Calculating time series...')
    % QPC
%     ts = a1*cos(2*pi*f1*t + phi1a) + a2*cos(2*pi*f2*t + phi2a) + ...
%     a3*cos((2*pi*(f1 + f2)*t)+(phi1a+phi2a)) + ...
%     a1*cos(2*pi*f3*t + phi4a) + a2*cos(2*pi*f4*t + phi5a) + ...
%     a3*cos((2*pi*(f3 + f4)*t)+(phi4a+phi5a)) + ...
%     a3*cos((2*pi*f6a*t)+phi6a); % + a5*v;
% 
%     ts = awgn(ts,1);

% phi1a = 0; 
% phi2a = 0;
% phi3a = -pi/2;

% disp('Calculating time series...')
% % QPC
%a3 = 2; % temp
%ts = a1*cos(2*pi*f1*t + phi1a) + a2*cos(2*pi*f2*t + phi2a) + a3*cos((2*pi*(f1 + f2)*t)+(phi1a+phi2a)) + ...
%     a3*cos((2*pi*f4a*t)+phi3a) + a5*v;
ts = a1*cos(2*pi*f1*t + phi1a) + a2*cos(2*pi*f2*t + phi2a) + a3*cos((2*pi*(f1 + f2)*t)+(phi1a+phi2a)) + a5*v;
 
 %ts = a1*sin(2*pi*f1*t + phi1a) + a2*sin(2*pi*f2*t + phi2a) + a3*sin((2*pi*(f1 + f2)*t)+(phi1a+phi2a)) + ...
 %    a3*sin((2*pi*f4a*t)+phi3a) + a5*v;
 
% no QPC
%a3 = 1;
%a3 = 5;

%ts = a1*cos(2*pi*f1*t + phi1a) + a2*cos(2*pi*f2*t + phi2a) + a3*cos((2*pi*(f1 + f2)*t)+phi3a) + a5*v;

%ts = a1*sin(2*pi*f1*t + phi1a) + a2*sin(2*pi*f2*t + phi2a) + a3*sin((2*pi*(f1 + f2)*t)+phi3a) + a5*v;
% no QPC, no frequency coupling
%a3 = 1;
%ts = a1*cos(2*pi*f1*t + phi1a) + a2*cos(2*pi*f2*t + phi2a) + a3*cos((2*pi*f3a*t)+phi3a) + a5*v;

% linear
%ts = a1*cos(2*pi*f1*t + phi1a) + a2*cos(2*pi*f2*t + phi2a) + a5*v;
%ts = a1*sin(2*pi*f1*t + phi1a) + a2*sin(2*pi*f2*t + phi2a) + a5*v;
% bilinear
%%%%%%%%%%%%%%%
% ts = ((a1*cos(2*pi*f1*t + phi1a)) .* (a2*cos(2*pi*f2*t + phi2a))) + a5*v;
%  tsLin = a1*cos(2*pi*f1*t + phi1a) + a2*cos(2*pi*f2*t + phi2a);
%  ts = ts + tsLin;

%ts = ((a1*sin(2*pi*f1*t + phi1a)) .* (a2*sin(2*pi*f2*t + phi2a))) + a5*v;
%%%%%%%%%%%%%%%%

%ts = ((a1*sin(2*pi*f1*t + phi1a)) .* (a2*sin(2*pi*f2*t + phi2a))) + a5*v;
% phase modulated

%ts = a1*cos(2*pi*f1*t + (a2 * cos(2*pi*f2*t + phi1a))) + a5*v;

%ts = a1*cos(2*pi*(f1 + (a2 * cos(2*pi*f2*t))).*t + phi1a) + a5*v;

disp('Done!...')
%disp(strcat('f1...',num2str(f1a),'...f2...',num2str(f2a),'...f3...',num2str(f3a)))
end
