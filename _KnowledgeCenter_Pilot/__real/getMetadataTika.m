

%function getMetadataTika

% folder_pdf = '/Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/Visual-Attention-Awareness/2. Retina Model/_ArticleTitle_v4/__PDF-2-DOWNLOADED/';
% folder_pdf = '/Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/__data/DOCS_PAPERS/__encrypted/';
folder_pdf = '//Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/Visual-Attention-Awareness/2. Retina Model/_ArticleTitle_v5/missing_pmc_javascript/';

folder_app = '/Users/joaodornas/Dropbox (joaodornas)/_Research/_CODES/reference-database/Reference-Database/__tools/Tika/';

% folder_meta = '/Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/Visual-Attention-Awareness/2. Retina Model/_ArticleTitle_v4/__MetaData-2-DOWNLOADED/';
% folder_meta = '/Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/__data/DOCS_PAPERS/__encrypted/';
folder_meta = '//Users/joaodornas/Dropbox (joaodornas)/_Research/_PROJECTS/Visual-Attention-Awareness/2. Retina Model/_ArticleTitle_v5/meta_pmc_javascript/';

tika_exe = 'tika-app-1.4.jar';
    
all_pdfs = dir(strcat(folder_pdf,'*.pdf'));

nPDFs = length(all_pdfs);

for iPDF=1:nPDFs

    sprintf('%s',all_pdfs(iPDF).name);

    % system(sprintf('java -jar "%s%s" -v -m --encoding=utf-8 "%s%s" > "%s%s-metadata-txt.html"',folder_app,tika_exe,folder_pdf,all_pdfs(iPDF).name,folder_meta,all_pdfs(iPDF).name(1:end-4)));

    % system(sprintf('java -jar "%s%s" -v -m -y --encoding=utf-8 "%s%s" > "%s%s-metadata-xmp.html"',folder_app,tika_exe,folder_pdf,all_pdfs(iPDF).name,folder_meta,all_pdfs(iPDF).name(1:end-4)));

    system(sprintf('java -jar "%s%s" -v -m -j --encoding=utf-8 "%s%s" > "%s%s-metadata-json.html"',folder_app,tika_exe,folder_pdf,all_pdfs(iPDF).name,folder_meta,all_pdfs(iPDF).name(1:end-4)));

end
