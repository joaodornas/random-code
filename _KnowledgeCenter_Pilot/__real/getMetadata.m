
app = 'tika-app-1.14.jar';
app_folder = '/Users/joaodornas/Downloads/MISSING/';

pdf_folder = '/Users/joaodornas/Downloads/MISSING/PDF/';

pdf_files = dir(strcat(pdf_folder,'*.pdf'));

nFiles = length(pdf_files);

%for iFile=1:nFiles
for iFile=1:nFiles
    
    sprintf('java -jar %s -m %s%s > %s-txt.html',app,pdf_folder,pdf_files(iFile).name,pdf_files(iFile).name(1:end-4))
    
%    system(sprintf('java -jar %s -m %s%s > %s-txt.html',app,pdf_folder,pdf_files(iFile).name,pdf_files(iFile).name(1:end-4))); 
%    system(sprintf('java -jar %s -j %s%s > %s-json.html',app,pdf_folder,pdf_files(iFile).name,pdf_files(iFile).name(1:end-4))); 
%    system(sprintf('java -jar %s -y %s%s > %s-xmp.html',app,pdf_folder,pdf_files(iFile).name,pdf_files(iFile).name(1:end-4))); 
    
end