function se2cen(dur,funct,split,limit,ctimes,nwsegs,txname,outpath)

dur = str2num(dur);
funct = str2num(funct);
split = str2num(split);
limit = str2num(limit);
ctimes = str2num(ctimes);
nwsegs = str2num(nwsegs);

% defaults--------------------
if dur == 0
	dur = 12;
end
if funct == 0
	funct = 1;
end
if limit == 0
	limit = 60;
end
if strcmp(txname,'')
	if funct == 1
		txname = 'segments.txt';
	else
		txname = 'times_raw.txt';
	end
end
%split = 1;
%ctimes = 1;
%nwsegs = 0;

% --------------------------
if strcmp(outpath,'.')
	outptc = 'centralTimesS';
        outptn = 'newSegments.txt';
else
	outptc = strcat(outpath,'centralTimesS');
	outptn = strcat(outpath,'newSegments.txt'); 
end

if funct == 1
	segs = load(txname);
	num = length(segs(:,1));

	ct = zeros(num,1);
	newsegs = zeros(num,2);
	disp('Calculating new segments and central times...')

	for i = 1:num
		ct(i) = segs(i,2) + (segs(i,2) - segs(i,1))/2;
		newsegs(i,1) = ct(i) - (dur/2);
		newsegs(i,2) = ct(i) + (dur/2);
	end
elseif funct ==2
	ct = load(txname);
	ct = ct(:,1);
	num = length(ct);
	newsegs = zeros(num,2);
        disp('Calculating segments...')

        for i = 1:num
                newsegs(i,1) = ct(i) - (dur/2);
                newsegs(i,2) = ct(i) + (dur/2);
        end

end

if ctimes
	dlmwrite(strcat(outptc,'.txt'),ct,'precision','%.6f')
end
if nwsegs
	dlmwrite(outptn,newsegs,'delimiter',' ','precision','%.6f')
end

if split
	disp('Writing central time segmented files(s)')
	%endval = ct(num);
	evn = mod(num,limit);
	if (num - limit) >= 0 
		slice = zeros(limit,1);
	else
		slice = zeros(num,1);
	end
	j = 1;
	k = 1;
	for i = 1:num
		slice(j) = ct(i);
		j = j + 1;
		if ~mod(i,limit)
			disp(strcat('Writing...',outptc,num2str(k),'.txt'))
			dlmwrite(strcat(outptc,num2str(k),'.txt'),slice,'precision','%.6f')
			k = k + 1;
			j = 1;
			if (num - (i + limit)) >= 0
                		slice = zeros(limit,1);
        		else
                		slice = zeros(num - i,1);
        		end

		end
	end
	if evn
		disp(strcat('Writing...',outptc,num2str(k),'.txt'))
                dlmwrite(strcat(outptc,num2str(k),'.txt'),slice,'precision','%.6f')
	end
end
disp('done!')

end
