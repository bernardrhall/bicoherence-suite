function checkRed(path)

%path = 'gaussian_T_0_8';

listing = dir(strcat(path,'/clusters*.dat'));

num = length(listing);
clusters = cell(num,3);

intsc = [];

for i = 1:num
	clusters{i,1} = load(strcat(path,'/',listing(i).name));
	clusters{i,2} = listing(i).name;
	clusters{i,3} = strsplit(listing(i).name,'_');
	clusters{i,3} = clusters{i,3}(1:length(clusters{i,3})-1);
	clusters{i,1}=unique(sort(clusters{i,1},2),'rows','stable'); % unique just a double check
	%clusters{i} = clusters{i}(a,:);
	%C = intersect(clusters{6},clusters{7},'rows')
end

for i = 1:num
	for j = 1:num
		if (i ~= j) && strcmp(char(clusters{i,3}),char(clusters{j,3}))
			C = intersect(clusters{i},clusters{j},'rows');
			[a,b]=size(C);
			if ~(a==0 || b==0)
				%C
				%intsc = [intsc;[i j]];
				[c,d]=size(intsc);
				%overlap{i,j} = char(strcat(clusters{i,2},{' '},clusters{j,2})); 
				if c==0 && d==0
					intsc=[i j];
				else
					intsc=[intsc;[i j]];
				end
			end
		end
	end
end

intsc=unique(sort(intsc,2),'rows','stable')
[x,~]=size(intsc);
for i=1:x
	overlap = char(strcat(clusters{intsc(i,1),2},{' '},clusters{intsc(i,2),2},'\n'))
	fid = fopen(strcat(path,'/intersection_',num2str(i),'.dat'), 'w');
	fprintf(fid, overlap);
	fclose(fid);
	%save(strcat('intersection_',num2str(i),'.dat'),'overlap','-ascii')
end

end
