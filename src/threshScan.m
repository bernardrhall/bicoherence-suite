function threshScan( fName, mName, threshMan, HbThresh, pathN)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    threshMan = str2num(threshMan);
    HbThresh = str2double(HbThresh);
    disp('Loading variables...')
    load(fName)
    disp('Done...')
    path = strcat(pathN,'/');
    mName = eval(mName);
    if threshMan
        bicohThresh = HbThresh;
    else
        bicohThresh = T;
        %bicohThresh = bicohThresh(1:18);
        bicohThresh = str2double(bicohThresh);
    end
    disp('Analyzing...')
    final_results = bcohCrit( bicohThresh, mName, freqBsp_F, str2double(T) );
    save(strcat(path,ifoCh,'_',num2str(GPS_central),'.mat'), '-v7.3') % save output
end

