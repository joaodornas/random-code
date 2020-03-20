

files = dir('*.*');

for ifile=3:length(files)
   
    old_name = files(ifile).name;
    
    old_name(1:18) = [];
    
    new_name = strcat('LH',old_name);
    
    movefile(files(ifile).name,new_name);
    
end