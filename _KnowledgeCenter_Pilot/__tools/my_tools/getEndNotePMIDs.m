
dd = dir('*.xml');

iiFile = 0;
for iFile=1:length(dd)
    
   [tree,~,~] = xml_read(dd(iFile).name);
   
   nRecords = length(tree.records.record);
   
   for iRecord=1:nRecords
       
       if isfield(tree.records.record(iRecord).urls,'related_DASH_urls')
           
           url = tree.records.record(iRecord).urls.pdf_DASH_urls.url;
           
           if isstruct(url) || iscell(url)
               
               for iURL=1:length(url)

                   iiFile = iiFile + 1;

                   tmp_pdf = url{iURL};
                   idx_slash = strfind(tmp_pdf,'/');

                   pdf_file{iiFile} = tmp_pdf((idx_slash(end)+1):end);

                   tmp_pmid = tree.records.record(iRecord).urls.related_DASH_urls.url.style.CONTENT;
                   idx_slash = strfind(tmp_pmid,'/');

                   pdf_pmid{iiFile} = tmp_pmid((idx_slash(end)+1):end);

               end
           
           else
              
               iiFile = iiFile + 1;

               tmp_pdf = url;
               idx_slash = strfind(tmp_pdf,'/');

               pdf_file{iiFile} = tmp_pdf((idx_slash(end)+1):end);

               tmp_pmid = tree.records.record(iRecord).urls.related_DASH_urls.url.style.CONTENT;
               idx_slash = strfind(tmp_pmid,'/');

               pdf_pmid{iiFile} = tmp_pmid((idx_slash(end)+1):end);

      
           end
           
       end
       
   end
    
end


nFiles = length(pdf_file);

pdfs_folder = '/Users/joaodornas/Downloads/__Papers-for-Mac/pdfs-with-NO-PMID-from-Papers-for-Mac-doing-EndNote/pdfs/';
pmids_folder = '/Users/joaodornas/Downloads/__Papers-for-Mac/pdfs-with-NO-PMID-from-Papers-for-Mac-doing-EndNote/pmids/';


for iFile=1:nFiles
   
    system(sprintf('mv %s%s %s%s.pdf',pdfs_folder,pdf_file{iFile},pmids_folder,pdf_pmid{iFile}));
    
end