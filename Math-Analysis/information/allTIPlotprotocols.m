
function allTIPlotprotocols(protocols,start_stimulus_time,end_stimulus_time,folder)

%% PROTOCOLS

if strcmp(protocols(end-2:end),'mat')
    
    registros = protocols;
    nReg = 1;
    
elseif strcmp(protocols(end-2:end),'txt')
    
    registros = importdata(protocols);
    nReg = length(registros);

end

%% BINS RESOLUTIONS

% size_of_bins = 1:(end_stimulus_time - start_stimulus_time);
% 
% j = 1;
% 
% for i=1:length(size_of_bins)
%     
%     if mod(max(size_of_bins),size_of_bins(i)) == 0
% 
%         bins_resolutions(j) = size_of_bins(i);
%         
%         j = j + 1;
%     
%     end
%     
% end

%% EACH LOOP FOR A PROTOCOL

for r=1:nReg
    
    if nReg == 1
        
        protocol = registros;
        name = registros(2:end-4);
        
    else
        
        protocol = registros{r};
        name = registros{r}(2:end-4);
        
    end
    
    disp(strcat('Protocol:',name));

    data = load(strcat(name,'-mutual-information.mat'));
    
    nConditions = size(data.info.Condition,2);  
    
    %% EACH LOOP FOR A CONDITION
    
    for i=1:nConditions
        
        disp(strcat('Condition:',int2str(i)));
        
        latency_time = getLatency(name,i);
         
        for WO=0:1
        
            disp(strcat('WOLatency:',int2str(WO)));
            
            mutualInfoData = data.info.Condition(i).Latency(WO+1).MutualInfo;
            
            f = figure;
            
            plot(1:length(mutualInfoData),mutualInfoData);
                                 
            print(f,strcat(folder,name,'-condition-',int2str(i),'-latency-',int2str(uint8(latency_time*WO*1000)),'-mutualInfo-all'));
            
            g = figure;
            
            plot(1:(length(mutualInfoData)/2),mutualInfoData(1:(length(mutualInfoData)/2)));
            
            print(g,strcat(folder,name,'-condition-',int2str(i),'-latency-',int2str(uint8(latency_time*WO*1000)),'-mutualInfo-bins-selected'));
     
            close all;
            
        end

    end
    
    
end
    
end


