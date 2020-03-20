
filename = input('Which is the file name prefix?,'s');

number = input('Which is the number of the last file?');

list = dir(strcat('./',filename,'*.mat'); % # assuming the file names are like stack1200.jpg

for idx = 1:length(list)   % # we go through every file name

    name = list(idx).name;   
   
    movefile(name,[filename num2str(number + idx) '.mat' ]);   % # we rename the file.

end