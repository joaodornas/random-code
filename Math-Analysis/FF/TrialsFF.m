function TB = TrialsFF(spike_times,start_time,end_time,time_bins)

nTrials = size(spike_times,1);

TB = zeros(length(time_bins),nTrials,length(time_bins));

for i=1:nTrials
    
   spikes = spike_times(i,:);
   
   spikes = spikes(spikes>0);
   
   spikes = spikes(spikes>=start_time/1000 & spikes<=end_time/1000);
   
   spikes = sort(spikes);
   
   for j=1:length(time_bins)
      
       nBins = j;
       
       time_step = (end_time - start_time)/nBins;

       for b=1:nBins
           
            TB(j,i,b) = numel(spikes(spikes>=(start_time/1000 + (b-1)*time_step/1000) & spikes<=(start_time/1000 + (b)*time_step/1000)));
           
       end

   end
 
end

end

