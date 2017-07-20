function matArchive(ipath)

ipath = 'hveto_1/output/'
dirString = strcat(ipath,'*.mat');
directory = dir(dirString);
nfiles = length(directory);

for k = 1:nfiles
	clearvars -except dirString directory nfiles k term ipath
	%name = '1135187444+L1:GDS-CALIB_STRAIN+16384+8_1135187448.mat'
	name = directory(k).name
	%ipath
	disp(strcat('Loading...',num2str(k),'...of...',num2str(nfiles)))
	load(strcat(ipath,name));
	disp('Done...')
	%m = matfile(name,'Writable',true);

	[x,y] = size(bisp);

	term = 2;
	
	disp('Reducing...')
	for i = 1:(y - 1)
		%disp(strcat('Column...',num2str(i)))
		%disp(term)
		bisp(term:x,i) = 0;
		term = term + 1;
	end
	disp('Done...')
	disp('Making sparse...')

	bisp = sparse(bisp);

	disp('Done...')

	disp('Saving...')
	save(strcat(ipath,name), '-regexp', '^(?!(directory|dirString|nfiles|x|y|term|name|k|i)$).','-v7.3');
	%save(name,'-v7.3');
	disp('Done!')
end

end
