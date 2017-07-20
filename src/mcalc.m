function [ multiplier, newtime ] = mcalc( time,overlap,segments,fs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Add = (((100 - overlap)/100)*time) * (segments);

newtime = time + Add;

nsamples = ceil((newtime/segments)* fs);
sampadd = ceil(time * fs)-nsamples;

multiplier = (sampadd + nsamples)/nsamples;

if rem(nsamples,2)
    nsamples = nsamples - 1;
end
nsamples = nsamples * multiplier;

olm  = fix(nsamples * overlap / 100);
nadvance = nsamples - olm;
if nadvance < 3
    newtime = newtime + 0.1;
end
% segT    = fix ( (round(newtime * fs) - olm) / nadvance);
% disp(strcat('Expected number of segments:...',num2str(segT)))

end

