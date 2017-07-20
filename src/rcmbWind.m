function [Bspec_1] = rcmbWind( winsize ,ndivnfft1, ndivnfft2 ,opwind, sections )
% rcmbWind( winsize ,nfft ,opwind )
disp('Recombining matrices...')

Bspec_1 = load('Bspec_1.mat');
Bspec_2 = load('Bspec_2.mat');
Bspec_1 = [Bspec_1.Bspec];
Bspec_2 = [Bspec_2.Bspec];

if sections == 4
    Bspec_3 = load('Bspec_3.mat');
    Bspec_4 = load('Bspec_4.mat');
    Bspec_3 = [Bspec_3.Bspec];
    Bspec_4 = [Bspec_4.Bspec];
end

%opwind = load('opwind.mat');
%opwind = d;

Bspec_1 = [Bspec_1 Bspec_2];

clearvars Bspec_2

if sections == 4
    Bspec_3 = [Bspec_3 Bspec_4];

    clearvars Bspec_4

    Bspec_1 = [Bspec_1;Bspec_3];

    clearvars Bspec_3
end
%----------------- frequency-domain smoothing ------------------------
disp('Applying window...')
% 
%[x,y] = size (Bspec_1);
%nfft = 2 * x;
%winsize = 10;
%winsize = winsize - rem(winsize,2) + 1;  % make it odd
lby2 = (winsize-1)/2;
Bspec_1 = conv2(Bspec_1,opwind);

if sections == 2
    %[x,y] = size (Bspec_1)
    lby2+ndivnfft1
    lby2+ndivnfft2
    Bspec_1 = Bspec_1(lby2+1:lby2+ndivnfft1,lby2+1:lby2+ndivnfft2);
elseif sections == 4
    Bspec_1 = Bspec_1(lby2+1:lby2+ndivnfft1,lby2+1:lby2+ndivnfft2);
end

end
