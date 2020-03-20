function allLatenciesForBack

% registroV = importdata('FB-v-protocols.txt');
% 
% for r=1:length(registroV)
%     
%    calcLatencyForBack(registroV{r});
%     
% end

registroG = importdata('FB-g-protocols.txt');

for r=6:6
    
   calcLatencyForBack(registroG{r});
    
end

end

