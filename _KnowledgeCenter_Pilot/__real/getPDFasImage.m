

% function getPDFasImage

% folder_pdf = '/Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/Visual-Attention-Awareness/2. Retina Model/_ArticleTitle_v4/__PDF-PapersForMac/';
folder_pdf = '/Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/__data/DOCS_PAPERS/__encrypted/PDF/';

folder_app = '/Users/joaodornas/Dropbox (joaodornas)/_Research/_CODES/reference-database/Reference-Database/__tools/PDFBox/';

PDF_exe = 'pdfbox-app-2.0.0.jar';

all_pdfs = dir(strcat(folder_pdf,'*.pdf'));

dpi = int2str(300);
    
nPDFs = length(all_pdfs);
    
for iPDF=1:nPDFs

    sprintf('java -jar "%s%s" PDFToImage -imageType jpg -dpi %s -color rgb -outputPrefix %s- "%s%s"',folder_app,PDF_exe,dpi,all_pdfs(iPDF).name(1:end-4),folder_pdf,all_pdfs(iPDF).name)

    system(sprintf('java -jar "%s%s" PDFToImage -imageType jpg -dpi %s -color rgb -outputPrefix %s- "%s%s"',folder_app,PDF_exe,dpi,all_pdfs(iPDF).name(1:end-4),folder_pdf,all_pdfs(iPDF).name));

end
    
% end
