
nError = size(all_files_error,1);
nMultiple = size(all_files_multiple,1);
nUnique = size(all_files_unique,1);
nNewMultiple = size(new_all_files_multiple,1);

% for iFile=1:nError
%     
%    system(sprintf('move "D:\\__PDFs-Papers-for-Mac\\MISSING\\PDF\\%s.pdf" "D:\\__PDFs-Papers-for-Mac\\MISSING\\PDF-DIVIDED\\ERROR\\%s.pdf"',num2str(all_files_error{iFile+1,1}),num2str(all_files_error{iFile+1,1}))); 
%     
% end

% for iFile=1:nUnique
%     
%    system(sprintf('move "D:\\__PDFs-Papers-for-Mac\\MISSING\\PDF\\%s.pdf" "D:\\__PDFs-Papers-for-Mac\\MISSING\\PDF-DIVIDED\\UNIQUE\\%s.pdf"',num2str(all_files_unique{iFile+1,1}),num2str(all_files_unique{iFile+1,1}))); 
%     
% end

% for iFile=1:nMultiple
%     
%    system(sprintf('move "D:\\__PDFs-Papers-for-Mac\\MISSING\\PDF\\%s.pdf" "D:\\__PDFs-Papers-for-Mac\\MISSING\\PDF-DIVIDED\\MULTIPLE\\%s.pdf"',num2str(all_files_multiple{iFile+1,1}),num2str(all_files_multiple{iFile+1,1}))); 
%     
% end


%%% CHANGE UNIQUE TO PMID ID
% for iFile=1:nUnique
%     
%    system(sprintf('move "D:\\__PDFs-Papers-for-Mac\\MISSING\\PDF-DIVIDED\\UNIQUE\\%s.pdf" "D:\\__PDFs-Papers-for-Mac\\MISSING\\PDF\\%s.pdf"',num2str(all_files_unique{iFile+1,1}),num2str(all_files_unique{iFile+1,2}))); 
%     
% end

for iFile=1:nNewMultiple
    
   system(sprintf('move "D:\\__PDFs-Papers-for-Mac\\MISSING\\PDF-DIVIDED\\MULTIPLE\\%s.pdf" "D:\\__PDFs-Papers-for-Mac\\MISSING\\PDF\\%s.pdf"',num2str(new_all_files_multiple{iFile+1,1}),num2str(new_all_files_multiple{iFile+1,2}))); 
    
end