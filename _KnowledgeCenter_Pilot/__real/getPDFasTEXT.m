

% function getPDFasImage

% folder_pdf = '/Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/__data/DOCS_PAPERS/__encrypted/Image2PDF/';
folder_pdf = '/Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/Visual-Attention-Awareness/2. Retina Model/_ArticleTitle_v5/missing_pmc/';

folder_app = '/Users/joaodornas/Dropbox (joaodornas)/_Research/_CODES/reference-database/Reference-Database/__tools/PDFBox/';

PDF_exe = 'pdfbox-app-2.0.0.jar';

all_pdfs = dir(strcat(folder_pdf,'*.pdf'));
    
nPDFs = length(all_pdfs);
    
for iPDF=1:nPDFs

    sprintf('java -jar "%s%s" ExtractText -html -encoding utf-8 %s %s.txt',folder_app,PDF_exe,all_pdfs(iPDF).name,all_pdfs(iPDF).name(1:end-4));

    system(sprintf('java -jar "%s%s" ExtractText -html -encoding utf-8 %s %s.txt',folder_app,PDF_exe,all_pdfs(iPDF).name,all_pdfs(iPDF).name(1:end-4)));

end
    
% end
