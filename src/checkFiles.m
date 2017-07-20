dirString = '*.mat';
directory = dir(dirString);
nfiles = length(directory);

for k = 1:nfiles
	clearvars -except dirString directory nfiles k
	name =  directory(k).name;
	ch = matfile(name);
	try
		char(ch.final_results);
		%disp(name);
	catch
		disp(name)
	end
end
