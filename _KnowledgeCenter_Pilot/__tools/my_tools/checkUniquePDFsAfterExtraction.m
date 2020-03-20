
folder_pdf = 'D:\__PROJECTS\__Retina\_ArticleTitle_v4\__PDF-1-PapersForMac\';
folder_html = 'D:\__PROJECTS\__Retina\_ArticleTitle_v4\__HTML-1-PapersForMac\';
folder_jpg = 'D:\__PROJECTS\__Retina\_ArticleTitle_v4\__PDFasImage-1-PapersForMac\';

pdfs = dir(strcat(folder_pdf,'*.pdf'));
html = dir(strcat(folder_html,'*.html'));
jpg = dir(strcat(folder_jpg,'*.jpg'));

PMID_html = [];
PMID_pdf = [];

for iPMID=1:length(pdfs)
   
    PMID_pdf(iPMID) = str2num(pdfs(iPMID).name(1:end-4));
    
end

for iPMID=1:length(html)
   
    PMID_html(iPMID) = str2num(html(iPMID).name(1:end-4));
    
end

idx_member = ismember(PMID_pdf,PMID_html);
idx_non_member = idx_member==0;
idx_pdf = find(idx_non_member);

missing_pdfs = PMID_pdf(idx_pdf);

for iPMID=1:length(jpg)
   
    name = jpg(iPMID).name(1:end-4);
    idx_slash = strfind(name,'-');
    PMID_jpg(iPMID) = str2num(name(1:idx_slash-1));
    
end

PMID_jpg = unique(PMID_jpg);

