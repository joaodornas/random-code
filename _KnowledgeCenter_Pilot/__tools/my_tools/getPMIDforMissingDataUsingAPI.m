
server_root = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/';

load('metaDataFromMissing.mat');

% nPII = length(pii);
% 
% for iPII=2:nPII
%     
%     DOI_string = strcat(server_root,'?ids=10.1016/',pii(iPII));
%     
%     urldata{iPII} = urlread(DOI_string{:});
%     
% end

% nDOI = length(doi);
% 
% for iDOI=2:nDOI
%     
%     DOI_string = strcat(server_root,'?ids',pii(iDOI));
%     
%     urldata{iDOI} = urlread(DOI_string{:});
%     
% end