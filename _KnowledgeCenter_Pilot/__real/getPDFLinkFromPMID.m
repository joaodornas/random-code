
function getPDFLinkFromPMID(query_label,PMIDs)

time_to_pause = 60;

nPMIDs = length(PMIDs);

for iPMID=1:nPMIDs
    
    disp(strcat(int2str((iPMID/nPMIDs)*100),'%'));
                          
    % PMID = PMIDs{iPMID};
    PMID = PMIDs(iPMID);
                          
    [medlineText, status] = urlread(strcat('https://www.ncbi.nlm.nih.gov/pubmed/',PMID));
                          
     if status

     	idx_FullTextSources = strfind(medlineText,'Full Text Sources');
        idx_Https = strfind(medlineText,'https');
        idx_Aspas = strfind(medlineText,'"');

        if isempty(idx_FullTextSources)

            articles{iPMID,1} = PMID;
            articles{iPMID,2} = 'NoField';

        else

        	articles{iPMID,1} = PMID;

           new_idx_Https = idx_Https - idx_FullTextSources(1);
           new_idx_Https(new_idx_Https<=0) = [];
           min_Https = min(new_idx_Https);

           idx_Aspas = idx_Aspas - (min_Https + idx_FullTextSources(1));
           idx_Aspas(idx_Aspas<=0) = [];
           min_Aspas = min(idx_Aspas);

           articles{iPMID,2} = medlineText(idx_FullTextSources(1)+min_Https:idx_FullTextSources(1)+min_Https+min_Aspas-1);

        end

     else
                            
           articles{iPMID,1} = PMID;
           articles{iPMID,2} = 'NoURL';
                            
     end
     
     save(strcat(query_label,'-','journalLink','.mat'));
     
     pause(time_to_pause);
                          
end           
                   
end