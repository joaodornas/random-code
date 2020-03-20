function calcKernelCC(registro,start_time,end_time)

Spass = load(registro);

nConditions = Spass.parameters.nconditions;

p = 1000;

for i=1:nConditions
   
    trials_label(i).labels = find(Spass.stimIds == i); 
    
end

spike_times = Spass.spike_times./ 32000;

for i=1:nConditions  
        
        trials_spikes = spike_times(trials_label(i).labels(:),:); 

        protocol = registro(2:end-4);

        nTrials(i) = min(size(trials_spikes,1),size(trials_spikes,2));

        allSpikes = [];
        
        for k=1:nTrials(i)

            spikes = trials_spikes(k,:);

            spikes = spikes(spikes>0);

            spikes = sort(spikes);

            spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));

            allSpikes = [allSpikes, spikes];

        end

        allSpikes = reshape(allSpikes.',[],1);

        allSpikes = sort(allSpikes);

        allSpikes = allSpikes(allSpikes>0);
        
        time_points = linspace(start_time/1000, end_time/1000, ( end_time/1000 - start_time/1000 ) * p);
        
        [kernel(i).density, kernel(i).timePoints, kernel(i).optimalBinWidth, kernel(i).WBinsTested, kernel(i).Cost, kernel(i).confb95] = ssvkernel(allSpikes,time_points);
    
end

datakernel = struct('kernel',kernel,'start_time',start_time,'end_time',end_time,'p',p);

save(strcat('/Users/joaodornas/Documents/_Research/_DATA/Center-Surround/kernel/',protocol,'-kernel-density-function-',int2str(start_time),'-',int2str(end_time)),'datakernel');

end

