function WordMatrix = gpuGetWordMatrix(spike_trains,nTrials,nBins,bin_size,start_stimulus_time)

nWords = nBins;

WordMatrix = [];

for t=1:nTrials
    
   train = spike_trains(t,:);
   
   sentence = [];
   
   for w=1:nWords
       
       word = numel(train(train>=start_stimulus_time/1000 + (w-1)*bin_size/1000 & train<=start_stimulus_time/1000 + (w)*bin_size/1000));
       
       sentence = [sentence, word];
       
       word = 0;
       
   end
   
   WordMatrix = [WordMatrix; sentence];
    
end

end

