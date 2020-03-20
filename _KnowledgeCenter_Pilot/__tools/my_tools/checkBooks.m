
function checkBooks(query_label)

% query_label = 'retina';
query_id = 1;
start_year = 1900;
end_year = 2016;

time_to_pause = 60;

full_path = 'L:\Dornas\Reference-Database\references\utf-8';

iPMID_Article = 0;
iPMID_Book = 0;
    
for iYear=start_year:end_year
    
    full_path_year = strcat(full_path,'\',query_label,'-',int2str(query_id),'\',int2str(iYear),'\');
    
    list = dir(strcat(full_path_year,'*.xml'));
    
    if ~isempty(list)
        
       nFiles = length(list);
       
       for iFile=1:nFiles
           
           disp(list(iFile).name);
           
           noLoaded = true;
           
           while noLoaded
               
               try

                    [xml_tree, ~, ~] = xml_read(strcat(full_path_year,list(iFile).name));
                
                    noLoaded = false;
                    
               catch ME
                   
                   msg = getReport(ME);
                   
                   fid = fopen(strcat(query_label,'.txt'),'a+');
                   fprintf(fid,'%s\n',msg);
                   fclose(fid);
                   
                   disp(msg);
                   disp(strcat('...pausing:',int2str(time_to_pause),'s'));
                   
                   pause(time_to_pause);
                   
               end
               
           end
           
            if ~isempty(xml_tree) && isstruct(xml_tree)
              
                apriori_fields = fieldnames(xml_tree);
              
                all_fields{iFile} = apriori_fields;
                
                if sum(strcmp('PubmedBookArticle',apriori_fields)) ~= 0
                    
                    nPublications = length(xml_tree.PubmedBookArticle);

                    if nPublications > 0
                        
                        for iPub=1:nPublications
                        
                            iPMID_Book = iPMID_Book + 1;
                        
                            Books.PMID_ID(iPMID_Book) = xml_tree.PubmedBookArticle(iPub).BookDocument.PMID;
                            Books.PMID_Year(iPMID_Book) = iYear;
                            Books.PMID_File{iPMID_Book} = list(iFile).name;
                            Books.PMID_Pub(iPMID_Book) = iPub;
                        
                        end
                        
                    end
                    
                end
                
                 if sum(strcmp('PubmedArticle',apriori_fields)) ~= 0
                    
                    nPublications = length(xml_tree.PubmedArticle);

                    if nPublications > 0
                        
                        for iPub=1:nPublications
                        
                            iPMID_Article = iPMID_Article + 1;
                        
                            Articles.PMID_ID(iPMID_Article) = xml_tree.PubmedArticle(iPub).MedlineCitation.PMID;
                            Articles.PMID_Year(iPMID_Article) = iYear;
                            Articles.PMID_File{iPMID_Article} = list(iFile).name;
                            Articles.PMID_Pub(iPMID_Article) = iPub;
                        
                        end
                        
                    end
                    
                end
                
            end
           
       end
       
    end
          
end

save(strcat('checkBooks','-',query_label,'.mat'),'Books','Articles','all_fields');

end