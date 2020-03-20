
function getLinks(label,query_id)

%%%% I NEED TO CONSULT PEOPLE IN PUBMED HOW OFTEN I CAN QUERY THE WEBSITE

full_path = 'L:\Dornas\Reference-Database\references';
% full_path = '/Volumes/DATA2/Dornas/Reference-Database/references';

start_year = 1900;
end_year = 2015;

for iYear=start_year:end_year
    
    disp(int2str(iYear)); 
    
    files = dir(strcat(full_path,'\',label,'-',int2str(query_id),'\',int2str(iYear),'\'));
    % files = dir(strcat(full_path,'/',label,'-',int2str(query_id),'/',int2str(iYear),'/'));
    
    if ~isempty(files)
        
        nFiles = length(files);
        
        for iFile=3:nFiles
            
            idx_article = 0;
            
            name = files(iFile).name;
            
            [tree, RootName, DOMnode] = xml_read(strcat(full_path,'\',label,'-',int2str(query_id),'\',int2str(iYear),'\',name));
            % [tree, RootName, DOMnode] = xml_read(strcat(full_path,'/',label,'-',int2str(query_id),'/',int2str(iYear),'/',name));
            
            if ~isempty(tree) && isstruct(tree)
              
              apriori_fields = fieldnames(tree);
              
              if ~strcmp('ERROR',apriori_fields)
          
                  nPublications = length(tree.PubmedArticle);

                  if nPublications > 0

                      for iPub=1:nPublications
                          
                          idx_article = idx_article + 1;
                          
                          PMID = tree.PubmedArticle(iPub).MedlineCitation.PMID.CONTENT;
                          
                          [medlineText, status] = urlread(strcat('https://www.ncbi.nlm.nih.gov/pubmed/',int2str(PMID)));
                          
                          if status

                              idx_FullTextSources = strfind(medlineText,'Full Text Sources');
                              idx_Https = strfind(medlineText,'https');
                              idx_Aspas = strfind(medlineText,'"');

                              if isempty(idx_FullTextSources)

                                  articles{idx_article,1} = PMID;
                                  articles{idx_article,2} = 'NoField';

                              else

                                  articles{idx_article,1} = PMID;

                                  new_idx_Https = idx_Https - idx_FullTextSources(1);
                                  new_idx_Https(new_idx_Https<=0) = [];
                                  min_Https = min(new_idx_Https);

                                  idx_Aspas = idx_Aspas - (min_Https + idx_FullTextSources(1));
                                  idx_Aspas(idx_Aspas<=0) = [];
                                  min_Aspas = min(idx_Aspas);

                                  articles{idx_article,2} = medlineText(idx_FullTextSources(1)+min_Https:idx_FullTextSources(1)+min_Https+min_Aspas-1);

                              end

                          else
                            
                              articles{idx_article,1} = PMID;
                              articles{idx_article,2} = 'NoURL';
                            
                          end
                          
                      end
                      
                  end
                  
              end
              
            end
             
            save(strcat(name(1:end-4),'-','journalLink','.mat'));
                
        end
        
    end
    
end
    
end