function calcKernelForBack(registro,start_time,end_time)


Spass = load(registro);

nConditions = 3;

latency_start_time = 500;

p = 1000;

for i=1:nConditions
   
    trials_label(i).labels = find(Spass.stimIds == i); 
    
end

spike_times = Spass.spike_times./ 32000;

for i=1:(nConditions + 1)
    
    if i == (nConditions + 1)
        
        j = 2;
        
    else
        
        j = i;
        
    end
        
        
        trials_spikes = spike_times(trials_label(j).labels(:),:); 

        protocol = registro(2:end-4);

        condition = j;

        latency_time(i) = getLatencyForBack(protocol,condition);

        nTrials(j) = min(size(trials_spikes,1),size(trials_spikes,2));

        allSpikes = [];

        if (j == 2), disp('Begin invert spike train...'); end
        
        for k=1:nTrials(j)

            spikes = trials_spikes(k,:);

            a_e = spikes(spikes>=(start_time/1000) & spikes<(latency_start_time/1000));

            after_ae = spikes(spikes>=(latency_start_time/1000) & spikes<(end_time/1000)) - latency_time(i);

            spikes = [];

            spikes = [a_e, after_ae];

            spikes = spikes(spikes>0);

            spikes = sort(spikes);

            spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));

            if (j == 2)

                bin_size = 1;

                nBins = (end_time - start_time)/bin_size;

                spikes_inverted = invertSpikeTrain(spikes,start_time,nBins,bin_size);
                spikes = spikes_inverted.spikes;

            end

            allSpikes = [allSpikes, spikes];

        end

        allSpikes = reshape(allSpikes.',[],1);

        allSpikes = sort(allSpikes);

        allSpikes = allSpikes(allSpikes>0);
        
        time_points = linspace(start_time/1000, end_time/1000, ( end_time/1000 - start_time/1000 ) * p);
        
        [kernel(i).density, kernel(i).timePoints, kernel(i).optimalBinWidth, kernel(i).WBinsTested, kernel(i).Cost, kernel(i).confb95] = ssvkernel(allSpikes,time_points);
    
end

datakernel = struct('kernel',kernel,'start_time',start_time,'end_time',end_time,'p',p,'latency_time',latency_time);

save(strcat('/Volumes/Data/DATA/Forward-Backward/kernel/',protocol,'-kernel-density-function-',int2str(start_time),'-',int2str(end_time)),'datakernel');

end

