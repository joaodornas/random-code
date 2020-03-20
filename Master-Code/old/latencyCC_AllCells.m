
function latencyCC_AllCells

cells_data_file = get_all_cells_CC;

nCells = length(cells_data_file);

latency_start_time = 500;

latency_end_time = latency_start_time + 300;

p = 10000;

start_time = 0;

end_time = 10000;

iiRegistro = 0;

for iCell=1:nCells
    
    nRegistro = length(cells_data_file(iCell).registro_video);
    
    for iRegistro=1:nRegistro
        
        iiRegistro = iiRegistro + 1;
        
        Spass = load(char(strcat(cells_data_file(iCell).registro_video(iRegistro).datafile,'.mat')));
        
        spike_times = Spass.spike_times ./ 32000;
        
        nConditions = max(Spass.stimIds);
        
        for iCondition=1:nConditions
            
            trials_label(iCondition).labels = find(Spass.stimIds == iCondition); 
    
            trials_spikes = spike_times(trials_label(iCondition).labels(:),:);

%             if iCondition == 3
% 
%                 latency_time(iCondition) = 0;
%                 
%                 latency_distribution(iiRegistro,3) = latency_time(3);
%             
%             else
            
                latency_time(iCondition) = latency(trials_spikes,latency_start_time,latency_end_time,p);
                
                latency_distribution(iiRegistro,iCondition) = latency_time(iCondition);

%             end
            
        end

        save(strcat('CC-Cell','-',int2str(iCell),'-','Video','-',int2str(cells_data_file(iCell).registro_video(iRegistro).video),'-',cells_data_file(iCell).registro_video(iRegistro).datafile,'-','latency','.mat'),'latency_time');
           
    end

end

save('CC-latency-distribution.mat','latency_distribution');

