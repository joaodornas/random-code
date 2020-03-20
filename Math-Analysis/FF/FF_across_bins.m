function results = FF_across_bins(TB,start_time,end_time,time_bins,nTrials)

for j=1:length(time_bins)
    
    nBins = j;
       
    %time_step = (end_time - start_time)/nBins;
    
    for i=1:nTrials
 
        for n=1:nBins
           
            bins(n) = TB(j,i,n);
            
        end
        
        bins(isnan(bins)) = 0;
        
        media_per_trial_across_bins(j,i) = mean(bins);
        
        variance_per_trial_across_bins(j,i) = var(bins);
        
        FF_per_trial_across_bins(j,i) = var(bins) ./ mean(bins) ;
        
        bins = [];
      
    end
    
    FF_per_trial_across_bins(isnan(FF_per_trial_across_bins)) = 0;
    
    media_per_trial_FF_per_trial_across_bins = mean(FF_per_trial_across_bins,2);
    
    media_per_bin_FF_per_trial_across_bins = mean(FF_per_trial_across_bins,1);
    
    media_total_FF_per_trial_across_bins = mean(mean(FF_per_trial_across_bins));
    
end

results.media_per_trial_across_bins = media_per_trial_across_bins;
results.variance_per_trial_across_bins = variance_per_trial_across_bins;
results.FF_per_trial_across_bins = FF_per_trial_across_bins;
results.media_per_trial_FF_per_trial_across_bins = media_per_trial_FF_per_trial_across_bins;
results.media_per_bin_FF_per_trial_across_bins = media_per_bin_FF_per_trial_across_bins;
results.media_total_FF_per_trial_across_bins = media_total_FF_per_trial_across_bins;

clear media_per_trial_across_bins;
clear variance_per_trial_across_bins;
clear FF_per_trial_across_bins;
clear media_FF_per_trial_across_bins;
clear media_per_trial_FF_per_trial_across_bins;
clear media_per_bin_FF_per_trial_across_bins;
clear media_total_FF_per_trial_across_bins;

end

