function calcLatencyCC(registro)

Spass = load(registro);

nConditions = Spass.parameters.nconditions;

start_stimulus_time = 500;

%end_time = 9500;

latency_start_time = start_stimulus_time;

latency_end_time = latency_start_time + 300;

p = 10000;

trials_label(1:nConditions) = struct('labels',[]);

for i=1:nConditions
   
    trials_label(i).labels = find(Spass.stimIds == i); 
    
end

spike_times = Spass.spike_times./ 32000;

latency_time = zeros(nConditions,1);

for i=1:nConditions
    
    disp(strcat('Condition:',int2str(i)));
    
    trials_spikes = spike_times(trials_label(i).labels(:),:); 
    
    alltrials = trials_spikes;

    alltrials = reshape(alltrials.',[],1);

    alltrials = alltrials.'; 
    
    alltrials = alltrials(alltrials>=0 & alltrials<=(latency_end_time/1000));

    disp('...Calculate Latency');    
    
    latency_time(i) = latency(alltrials,latency_start_time,latency_end_time,p); 

    clear alltrials;
    
end

name = registro(2:end-4);

results = struct('name',name,'latency',latency_time);

save(strcat('/Users/joaodornas/Documents/_Research/_DATA/Center-Surround/latency/',name,'-latencies-for-conditions'),'results');

end

