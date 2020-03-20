
%%% CHECK PMIDs ON SEARCH PMIDs by TITLE

files = dir('*.html');

nFiles = length(files);

%%% UNIQUE PMIDS

% iiFile = 0;
% for iFile=1:nFiles
%     
%    fid = fopen(files(iFile).name,'r');
%    
%    whole = fscanf(fid,'%s');
%    
%    idx_PMID_start = strfind(whole,'<PMID'); 
%    idx_PMID_end = strfind(whole,'</PMID');
% 
%    nPMIDs = length(idx_PMID_start);
%    
%    all_files{1,1} = 'MISSING';
%    all_files{1,2} = 'PMID'; 
%      
% %    for iPMID=1:nPMIDs
% %     
% %        all_files{iFile,iPMID} = whole(idx_PMID_start(iPMID)+17:idx_PMID_end(iPMID)-1);
% %        
% %    end
% 
%     if nPMIDs == 1
%         
%         iiFile = iiFile + 1;
%         
%         all_files{iiFile+1,1} = files(iFile).name(1:end-5);
%         
%         all_files{iiFile+1,2} = whole(idx_PMID_start+17:idx_PMID_end-1);
%         
%     end
%    
%    fclose(fid);
%    
% end

%%% EMPTY RESULTS

% iiFile = 0;
% for iFile=1:nFiles
%     
%    fid = fopen(files(iFile).name,'r');
%    
%    whole = fscanf(fid,'%s');
%    
%    idx_ERROR = strfind(whole,'Emptyresult-nothingtodo'); 
% 
%    all_files{1,1} = 'ERROR';
%  
%     if ~isempty(idx_ERROR)
%         
%         iiFile = iiFile + 1;
%         
%         all_files{iiFile+1,1} = files(iFile).name(1:end-5);
%         
%     end
%    
%    fclose(fid);
%    
% end

%%% MULTIPLE PMIDS

iiFile = 0;
for iFile=1:nFiles
    
   fid = fopen(files(iFile).name,'r');
   
   whole = fscanf(fid,'%s');
   
   idx_PMID_start = strfind(whole,'<PMIDVersion="1">'); 
   idx_PMID_end = strfind(whole,'</PMID');

   nPMIDs = length(idx_PMID_start);
   
   all_files{1,1} = 'MISSING';
   all_files{1,2} = 'PMID'; 
  
   if nPMIDs > 1
     
       iiFile = iiFile + 1;
       
       all_files{iiFile+1,1} = files(iFile).name(1:end-5);
       
       for iPMID=1:nPMIDs

           all_files{iiFile+1,iPMID+1} = whole(idx_PMID_start(iPMID)+17:idx_PMID_end(iPMID)-1);

       end

   end
   
   fclose(fid);
   
end


