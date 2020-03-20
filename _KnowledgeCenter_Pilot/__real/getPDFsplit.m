

function getPDFsplit

% folder_pdf = 'L:\Dornas\Reference-Database\pdfs\pdfs-with-PMID-from-ReadCube\Tika-PDFSplit\';
folder_pdf = '/Users/joaodornas/Downloads/MISSING/PDF-SPLIT/';

% folder_app = 'L:\Dornas\Reference-Database\Install\';
folder_app = '/Users/joaodornas/Downloads/';

PDF_exe = 'pdfbox-app-2.0.5.jar';

all_pdfs = dir(strcat(folder_pdf,'*.pdf'));

nPDFs = length(all_pdfs);

for iPDF=1:nPDFs

    system(sprintf('java -jar %s%s PDFSplit %s%s',folder_app,PDF_exe,folder_pdf,all_pdfs(iPDF).name));

end

end