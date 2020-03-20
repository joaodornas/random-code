

function getPDFimages

folder_pdf = 'L:\Dornas\Reference-Database\pdfs\pdfs-with-PMID-from-Papers-for-Mac-Tika-ExtractImages\';

folder_app = 'L:\Dornas\Reference-Database\Install\';

PDF_exe = 'pdfbox-app-2.0.5.jar';
   
all_pdfs = dir(strcat(folder_pdf,'\*.pdf'));

nPDFs = length(all_pdfs);

for iPDF=1:nPDFs

    system(sprintf('java -jar %s%s ExtractImages %s%s',folder_app,PDF_exe,folder_pdf,all_pdfs(iPDF).name));

end

end