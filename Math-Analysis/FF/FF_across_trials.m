function results = FF_across_trials(TB,start_time,end_time,time_bins,nTrials)

for j=1:length(time_bins)
    
    nBins = j;
       
    %time_step = (end_time - start_time)/nBins;
    
    for n=1:nBins
        
        for i=1:nTrials
            
            bins(i) = TB(j,i,n);
            
        end
        
        bins(isnan(bins)) = 0;
        
        media_across_trials_per_time_bin(j,n) = mean(bins);
        
        variance_across_trials_per_time_bin(j,n) = var(bins);
        
        FF_across_trials_per_time_bin(j,n) = var(bins) ./ mean(bins) ;      
        
        bins = [];
        
    end
    
    FF_across_trials_per_time_bin(isnan(FF_across_trials_per_time_bin)) = 0;
    
    media_per_trial_FF_across_trials_per_time_bin = mean(FF_across_trials_per_time_bin,2);
    
    media_per_bin_FF_across_trials_per_time_bin = mean(FF_across_trials_per_time_bin,1);
    
    media_total_FF_across_trials_per_time_bin = mean(mean(FF_across_trials_per_time_bin));

end

results.media_across_trials_per_time_bin = media_across_trials_per_time_bin;
results.variance_across_trials_per_time_bin = variance_across_trials_per_time_bin;
results.FF_across_trials_per_time_bin = FF_across_trials_per_time_bin;
results.media_per_trial_FF_across_trials_per_time_bin = media_per_trial_FF_across_trials_per_time_bin;
results.media_per_bin_FF_across_trials_per_time_bin = media_per_bin_FF_across_trials_per_time_bin;
results.media_total_FF_across_trials_per_time_bin = media_total_FF_across_trials_per_time_bin;

clear media_across_trials_per_time_bin;
clear variance_across_trials_per_time_bin;
clear FF_across_trials_per_time_bin;
clear media_FF_across_trials_per_time_bin;
clear media_per_trial_FF_across_trials_per_time_bin;
clear media_per_bin_FF_across_trials_per_time_bin;
clear media_total_FF_across_trials_per_time_bin;

end

