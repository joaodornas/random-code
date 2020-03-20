function spike_trains = takeLatencyOff(all_trials,nTrials,start_stimulus_time,end_stimulus_time,latency_start_time,latency_time)


for t=1:nTrials

    Trial(t).spike_train = all_trials(t,:);

    a_e = Trial(t).spike_train(Trial(t).spike_train>=(0/1000) & Trial(t).spike_train<(latency_start_time/1000));

    after_ae = Trial(t).spike_train(Trial(t).spike_train>=(latency_start_time/1000) & Trial(t).spike_train<(end_stimulus_time/1000)) - latency_time;

    Trial(t).spike_train = [];

    Trial(t).spike_train = [a_e, after_ae];

    Trial(t).spike_train = Trial(t).spike_train(Trial(t).spike_train>0);

    Trial(t).spike_train = Trial(t).spike_train(Trial(t).spike_train>=(start_stimulus_time/1000) & Trial(t).spike_train<=(end_stimulus_time/1000));
    
    nSpikes(t) = numel(Trial(t).spike_train);

end

maxNSpikes = max(nSpikes);

spike_trains = zeros(t,maxNSpikes);

for t=1:nTrials
    
    spike_trains(t,:) = [Trial(t).spike_train, zeros(1,maxNSpikes - numel(Trial(t).spike_train))];
    
end

end

