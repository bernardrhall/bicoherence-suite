function [ts] = sig_gen1 (time,  fs, f1, f2, a1, a2, a5,v)

% fs = 1024; %higher sample rate gives greater resolution, but also increases processing time

% if ~f1
%     f1 = 1 + ((fs/2)-1)*0.2*rand;
% end

rng('shuffle')  %reseed based on time

inc = 1/fs;

t=0:inc:time;

% signals = 3;

%phase_mode = 3; % 0 for normal, 1 for phase-modulated, 2 for linear
% f3 = 102.4;
% f4 = 204.8;
% f6a = 450;

phase_mode = 2;

if phase_mode == 1
  
phi1a = -pi; 
phi2a = -pi;
    while phi1a == -pi
        phi1a = -pi + (pi+pi)*rand;
    end
    while phi2a == -pi
        phi2a = -pi + (pi+pi)*rand;
    end
elseif phase_mode == 2
    phi1a = 0.7; 
    phi2a = 0.3;
    %phi3a = 0.5;
end

disp('Calculating time series...')
% QPC

% ts1 = (1-a1) * (cos(2*pi*f1*t + phi1a) + cos(2*pi*f2*t + phi2a));%********
% ts2 = a1 * (cos(2*pi*f1*t + phi1a) + cos(2*pi*f2*t + phi2a)).^2;%**********

% t1 = (a1 * cos(2*pi*f1*t + phi1a)) + (a2 * cos(2*pi*f2*t + phi2a));
% t2 = 0.1 * cos(4*pi*f1*t + 2*phi1a);
% t3 = 0.1 * cos(4*pi*f2*t + 2*phi2a);
% t4 = 0.5 * cos(2*pi*(f1+f2)*t + (phi1a+phi2a));
% t5 = 0.1 * cos(2*pi*(f1-f2)*t + (phi1a-phi2a));
phi_fm = 0;
fm = 0.125;
am = 0.75;
intercept = am + f2;
tfm = (am*cos(2*pi*fm*t + phi_fm)+intercept);
phivec = ones(length(t),1);
phivec = (phivec * phi2a)';

phi_fm1 = 0;
fm1 = 0.1;
am1 = 2;
intercept1 = f1;
tfm1 = (am1*cos(2*pi*fm1*t + phi_fm1)+intercept1);
phivec1 = ones(length(t),1);
phivec1 = (phivec1 * phi1a)';

%ts1 = a1*((cos(2*pi*(tfm1.*t) + phivec1) + cos(2*pi*(tfm.*t) + phivec)).^2);%********
%ts2 = a2*(cos(2*pi*(tfm1.*t) + phivec1) + cos(2*pi*(tfm.*t) + phivec));%**********

ts1 = a1*((cos(2*pi*f1*t + phi1a) + cos(2*pi*f2*t + phi2a)).^2);%********
ts2 = a2*(cos(2*pi*f1*t + phi1a) + cos(2*pi*f2*t + phi2a));%**********
if v == 0
    v = zeros(length(ts1),1)';
elseif v == 1
    N = fs * time;
    v=randn(1,N+1);
end
%ts2 = zeros(length(ts1),1)';

% ts = 1 + cos(2*pi*f1*t + phi1a) + cos(2*pi*f2*t + phi2a) + ...
%     0.5*cos(4*pi*f1*t + 2*phi1a) + 0.5*cos(4*pi*f2*t + 2*phi2a) + ...
%     cos(2*pi*(f1 + f2)*t + (phi1a + phi2a)) + ...
%     cos(2*pi*(f1 - f2)*t + (phi1a - phi2a)) + a5*v;

% ts1 = a1*cos(2*pi*f1*t + phi1a);
% ts2 = (a1*cos(2*pi*f1*t + phi1a)).^2;

ts = ts2 + ts1 + a5*v;
%ts = t1 + t2 + t3 + t4 + t5 + a5*v;

disp('Done!...')
%disp(strcat('f1...',num2str(f1a),'...f2...',num2str(f2a),'...f3...',num2str(f3a)))
end