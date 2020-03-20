function [results TB] = FF(spike_times,start_time,end_time,time_bins)

nTrials = size(spike_times,1);

TB = zeros(time_bins,nTrials,(end_time - start_time)/time_bins(end));

for i=1:nTrials
    
   spikes = spike_times(i,:);
   
   spikes = spikes(spikes>0);
   
   spikes = spikes(spikes>=start_time/1000 & spikes<=end_time/1000);
   
   spikes = sort(spikes);
   
   for j=1:length(time_bins)
      
       nBins = (end_time - start_time)/time_bins(j);
       
       bin = zeros(1,nBins);
       
       for b=1:nBins
           
            TB(j,i,b) = length(spikes(spikes>=(start_time/1000 + (b-1)*time_bins(j)/1000) & spikes<=(start_time/1000 + (b)*time_bins(j)/1000)));
           
       end

   end
 
end


for j=1:length(time_bins)
    
    nBins = (end_time - start_time)/time_bins(j);
    
    for i=1:nTrials
 
        for n=1:nBins
           
            bins(n) = TB(j,i,n);
            
        end
        
        media_per_trial_across_bins(j,i) = mean(bins);
        
        variance_per_trial_across_bins(j,i) = var(bins);
        
        FF_per_trial_across_bins(j,i) = variance_per_trial_across_bins(j,i) ./ media_per_trial_across_bins(j,i) ;
        
        bins = [];
      
    end
    
    FF_per_trial_across_bins(isnan(FF_per_trial_across_bins)) = [];
    
    media_FF_per_trial_across_bins = mean(FF_per_trial_across_bins);
    
    for n=1:nBins
        
        for i=1:nTrials
            
            bins(i) = TB(j,i,n);
            
        end
        
        media_across_trials_per_time_bin(j,n) = mean(bins);
        
        variance_across_trials_per_time_bin(j,n) = var(bins);
        
        FF_across_trials_per_time_bin(j,n) = variance_across_trials_per_time_bin(j,n) ./ media_across_trials_per_time_bin(j,n) ;      
        
        bins = [];
        
    end
    
    FF_across_trials_per_time_bin(isnan(FF_across_trials_per_time_bin)) = [];
    
    media_FF_across_trials_per_time_bin = mean(FF_across_trials_per_time_bin);

end

results.media_per_trial_across_bins = media_per_trial_across_bins;
results.variance_per_trial_across_bins = variance_per_trial_across_bins;
results.FF_per_trial_across_bins = FF_per_trial_across_bins;
results.media_FF_per_trial_across_bins = media_FF_per_trial_across_bins;

results.media_across_trials_per_time_bin = media_across_trials_per_time_bin;
results.variance_across_trials_per_time_bin = variance_across_trials_per_time_bin;
results.FF_across_trials_per_time_bin = FF_across_trials_per_time_bin;
results.media_FF_across_trials_per_time_bin = media_FF_across_trials_per_time_bin;

end

