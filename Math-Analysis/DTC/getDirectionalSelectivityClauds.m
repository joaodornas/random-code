function [m, A, preferida, rate] = getDirectionalSelectivityClauds(DTC,start_time,end_time)


Spass = load(DTC);

nConditions = 16;

angles = 22.5;

bin_size = (end_time - start_time) / 1000;

for i=1:nConditions
   
    trials_label(i).labels = find(Spass.gas056b01_1A.datafiles.stimIds == i); 
    
end

spike_times = Spass.gas056b01_1A.datafiles.spike_times./ 32000;

for i=1:nConditions
    
    trials_spikes = spike_times(trials_label(i).labels(:),:);
    
    nTrials = size(trials_spikes,1);
    
    spike_train = reshape(trials_spikes.',[],1);
       
    spike_train = spike_train(spike_train>0);
    
    ae = spike_train(spike_train>=0/1000 & spike_train<=start_time);
       
    spike_train = spike_train(spike_train>=start_time/1000 & spike_train<=end_time);
    
    all_m(i) = numel(ae)/( (nTrials) * bin_size );
    
    rate(i) = numel(spike_train)/( (nTrials) * bin_size );
    
end

idx = find(rate==max(rate));

preferida = (idx - 1) * angles;

A = max(rate);

m = mean(all_m);


end

