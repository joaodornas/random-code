function spike_train = convertLabData( nconditions, stimIds, spike_times, taxa_amostragem )

spike_train = zeros(size(spike_times,2));

for c=17:32
    
    trials = find(stimIds == c);
    ntrials = size(trials,2);
   
    for t=1:ntrials
        
       i = trials(t);
           
       spike_train = spike_times(i,1:size(spike_times,2));
       spike_train(spike_train==-1) = [];
       spike_train = spike_train / taxa_amostragem;
       %spike_train(spike_train<1) = [];
       %spike_train(spike_train>5) = [];
       dlmwrite('spike_trains.txt',spike_train,'delimiter',' ','-append');
       
    end
    
end

end

