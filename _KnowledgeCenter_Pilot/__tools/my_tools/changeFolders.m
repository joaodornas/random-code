
source = '/Users/joaodornas/Downloads/__Papers-for-Mac/pdfs-with-PMID-from-MISSING/pdfs/';

target = '/Users/joaodornas/Downloads/__Papers-for-Mac/pdfs-with-PMID-from-MISSING/found/';

nFiles = length(pii_file);

for iFile=2:nFiles
    
   file_name = pii_file{iFile}; 
   
   system(sprintf('mv %s%s.pdf %s%s.pdf',source,file_name(1:end-5),target,file_name(1:end-5))); 
    
end
