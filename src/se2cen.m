dur = 12;
funct = 1
split = 1;
limit = 60;

ctimes = 1;
nwsegs = 0;

if funct == 1
	segs = load('segments.txt');
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
	ct = load('times_raw.txt');
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
	dlmwrite('centralTimesS.txt',ct,'precision','%.6f')
end
if nwsegs
	dlmwrite('newSegments.txt',newsegs,'delimiter',' ','precision','%.6f')
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
			disp(strcat('Writing...','centralTimesS',num2str(k),'.txt'))
			dlmwrite(strcat('centralTimesS',num2str(k),'.txt'),slice,'precision','%.6f')
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
		disp(strcat('Writing...','centralTimesS',num2str(k),'.txt'))
                dlmwrite(strcat('centralTimesS',num2str(k),'.txt'),slice,'precision','%.6f')
	end
end
disp('done!')
