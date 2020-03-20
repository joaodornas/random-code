function plotKernel(spike_times,stimIds,condition,start_time,end_time,p,cur_color)

labels = find(stimIds == condition); 

spike_times = spike_times./ 32000;

trials = spike_times(labels(:),:);
        
trials = reshape(trials.',[],1);

trials = trials.';

trials = sort(trials);

trials = trials(trials>0);    

trials = trials(trials>=(start_time/1000) & trials<=(end_time/1000));
 
time_points = linspace(start_time/1000, end_time/1000, (end_time - start_time)/1000)*p;
    
[kernel.density, kernel.timePoints, kernel.optimalBinWidth, kernel.WBinsTested, kernel.Cost, kernel.confb95] = ssvkernel(trials,time_points);

f = figure;
plot(kernel.timePoints,kernel.density,cur_color);
print(f,'-depsc',strcat('condition-',int2str(condition)));
    
end

