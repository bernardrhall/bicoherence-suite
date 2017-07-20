function addsignals(sig1,path,outpath,count)

sqpc=load(sig1);
%a = load('add.ini');
count = str2num(count);
listing = dir(strcat(path,'/*.txt'));

num = length(listing);
%num = 1;
for i = 1:num
	%listing(i)
	k = i + ((count - 1)*10);
	sigg = load(strcat(path,'/',listing(i).name));
	N = length(sigg(:,2));
	for j = 1:N
		sigg(j,2) = sigg(j,2) + sqpc(j);  
	end
	disp('Writing to file...')
    	outfile = fopen(strcat(outpath,'/QPC_',num2str(k),'.txt'),'w');

    	for j = 1:N
        	fprintf(outfile,'%s\r\n',char(strcat(num2str(sigg(j,1),16),{' '},num2str(sigg(j,2),16))));
    	end
    	fclose(outfile);
end

end
