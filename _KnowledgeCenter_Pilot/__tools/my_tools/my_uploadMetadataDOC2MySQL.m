
conn = database('docs_papers','joaodornas','Senh@1111','Vendor','MySQL','Server','141.44.42.20','PortNumber',7979);

full_path_pdf = 'L:\Dornas\Reference-Database\pdfs\pdfs';
full_path_metadata = 'L:\Dornas\Reference-Database\pdfs\metadata';

list = dir(strcat(full_path_pdf,'*.pdf'));

nFiles = length(list);

for iFile=1:nFiles
   
    disp(strcat('iFile:',int2str(iFile),':',list(iFile).name));
    
    
    
    
end