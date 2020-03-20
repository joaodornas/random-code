filenames = dir;

filenames(1:2) = [];

q = length(filenames);

xlswrite('onexls',{'Slice' 'Type 1' 'Type 2' 'Type 3' 'Type 4'},'A1:E1');

title = cell(1,q);
title{1,q} = [];
line = zeros(q,4);
newline = zeros(q,4);

for i=1:q
   
    file = filenames(i,1).name;
    
    title{i} = filenames(i,1).name(1:end-4);
    
    line(i,:) = xlsread(file,'B2:E2');
    
    xlswrite('onexls',{title{i}},strcat('A',int2str(i),':A',int2str(i)));
    
    xlswrite('onexls',line(i,:),strcat('B',int2str(i),':E',int2str(i)));

end


