function calcLatencyForBack(registro)

Spass = load(registro);

% nConditions = 3;

start_time = 500;

latency_start_time = start_time;

latency_end_time = latency_start_time + 300;

p = 10000;

nConditions = Spass.parameters.nconditions;

for i=1:nConditions
   
    trials_label(i).labels = find(Spass.stimIds == i); 
    
end

spike_times = Spass.spike_times./ 32000;

for i=1:nConditions
    
    trials_spikes = spike_times(trials_label(i).labels(:),:); 
    
    alltrials = trials_spikes;

    alltrials = reshape(alltrials.',[],1);

    alltrials = alltrials.'; 

    if i ~= 3
        
        latency_time(i) = latency(alltrials,latency_start_time,latency_end_time,p); 
        
    else
        
        latency_time(i) = 0;
        
    end

    clear alltrials;
    
end

name = registro(2:end-4);

results = struct('name',name,'latency',latency_time);

folder = 'C:\Users\Dornas\Documents\Dropbox\_Research\_DATA\';

save(strcat(folder,'Forward-Backward\latency\',name,'-latencies-for-conditions'),'results');



end

