
function getFieldNames(label,query_id)

field_names_MedlineCitation = cell.empty;
field_names_MedlineCitation_Article = cell.empty;
field_names_MedlineCitation_Article_Abstract = cell.empty;
field_names_MedlineCitation_Article_Journal = cell.empty;

full_path = 'L:\Dornas\Reference-Database\references';
% full_path = '/Volumes/DATA2/Dornas/Reference-Database/references';

start_year = 1900;
end_year = 2015;

for iYear=start_year:end_year
    
    disp(int2str(iYear));
    
    all_Medline = cell.empty;
    all_Medline_Article = cell.empty;
    all_Medline_Article_Abstract = cell.empty;
    all_Medline_Article_Journal = cell.empty;
    
    % full_path_year = strcat(full_path,'/',label,'-',int2str(query_id),'/',int2str(iYear),'/');
    full_path_year = strcat(full_path,'\',label,'-',int2str(query_id),'\',int2str(iYear),'\');
    
    list = dir(strcat(full_path_year,'*'));
    
    if ~isempty(list)
        
       nFiles = length(list);

       for iFile=3:nFiles
           
          disp(strcat('iFile:',int2str(iFile)));
           
          [tree, RootName, DOMnode] = xml_read(strcat(full_path_year,list(iFile).name));
          
          if ~isempty(tree) && isstruct(tree)
              
              apriori_fields = fieldnames(tree);
              
              if ~strcmp('ERROR',apriori_fields)
   
                nPublications = length(tree.PubmedArticle);     

                if nPublications > 0

                    for iPub=1:nPublications
                        
                        % if mod(iPub,10) == 0; disp(strcat('iPub:',int2str(iPub))); end

                        field_names_MedlineCitation{end+1} = fieldnames(tree.PubmedArticle(iPub).MedlineCitation);
                        
                        field_names_MedlineCitation_Article{end+1} = fieldnames(tree.PubmedArticle(iPub).MedlineCitation.Article);
                        
                        tmp_fields = fieldnames(tree.PubmedArticle(iPub).MedlineCitation.Article);
                        
                        if sum(strcmp('Abstract',tmp_fields)) ~= 0
                        
                            field_names_MedlineCitation_Article_Abstract{end+1} = fieldnames(tree.PubmedArticle(iPub).MedlineCitation.Article.Abstract);
                            
                        end
                        
                        if sum(strcmp('Journal',tmp_fields)) ~= 0
                        
                            field_names_MedlineCitation_Article_Journal{end+1} = fieldnames(tree.PubmedArticle(iPub).MedlineCitation.Article.Journal);
                            
                        end

                    end

                end
                
              end
          
          end
          
       end
       
    end
    
    k=0;
    for iMed=1:length(field_names_MedlineCitation)
        ans = field_names_MedlineCitation{iMed};
        for j=1:length(ans)
            k=k+1;
            all_Medline{k} = ans{j};
        end
    end

    k=0;
    for iMed=1:length(field_names_MedlineCitation_Article)
        ans = field_names_MedlineCitation_Article{iMed};
        for j=1:length(ans)
            k=k+1;
            all_Medline_Article{k} = ans{j};
        end
    end

    k=0;
    for iMed=1:length(field_names_MedlineCitation_Article_Abstract)
        ans = field_names_MedlineCitation_Article_Abstract{iMed};
        for j=1:length(ans)
            k=k+1;
            all_Medline_Article_Abstract{k} = ans{j};
        end
    end
    
    k=0;
    for iMed=1:length(field_names_MedlineCitation_Article_Journal)
        ans = field_names_MedlineCitation_Article_Journal{iMed};
        for j=1:length(ans)
            k=k+1;
            all_Medline_Article_Journal{k} = ans{j};
        end
    end
    
    [u_Medline,c] = unique(all_Medline);
    [u_Medline_Article,c] = unique(all_Medline_Article);
    [u_Medline_Article_Abstract,c] = unique(all_Medline_Article_Abstract);
    [u_Medline_Article_Journal,c] = unique(all_Medline_Article_Journal);

    save(strcat('unique','-',label,'-',int2str(query_id),'-',int2str(iYear),'.mat'),'u_Medline','u_Medline_Article','u_Medline_Article_Abstract','u_Medline_Article_Journal');
    
end

end
