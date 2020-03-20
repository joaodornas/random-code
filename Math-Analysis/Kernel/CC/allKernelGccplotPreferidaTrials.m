function allKernelGccplotPreferidaTrials

registro{1} = '_nsp021a03_1a.mat';
registro{2} = '_nsp023a03_1b.mat';

preferida(1) = 4;
preferida(2) = 13;

start_time = 0;
end_time = 3000;

nConditions = 16;

p = 1000;

for r=1:length(registro)

    disp(registro{r});
    
    protocol = registro{r}(2:end-4);
    
    Spass = load(registro{r});
    
    labels = find(Spass.stimIds == preferida(r)); 
    
    spike_times = Spass.spike_times./ 32000;
    
    nTrials = numel(labels);
    
    trials_spikes = spike_times(labels(:),:); 

    for k=1:nTrials

        spikes = trials_spikes(k,:);

        spikes = spikes(spikes>0);

        spikes = sort(spikes);

        spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));
        
        time_points = linspace(start_time/1000, end_time/1000, ( end_time/1000 - start_time/1000 ) * p);
        
        [kernel(k).density, kernel(k).timePoints, kernel(k).optimalBinWidth, kernel(k).WBinsTested, kernel(k).Cost, kernel(k).confb95] = ssvkernel(spikes,time_points);

    end       
    
    filepath = strcat('/Users/joaodornas/Documents/_Research/_DATA/Center-Surround/kernel/',protocol,'-kernel-density-function-',int2str(start_time),'-',int2str(end_time),'-per-preferida-trial');
    
    datakernel = struct('kernel',kernel,'start_time',start_time,'end_time',end_time,'p',p);

    save(filepath,'datakernel');
    
    H = figure;
    
    for i=1:nTrials
        
        subplot(nTrials,1,i);
        
        plot(time_points,kernel(i).density,'r');
        
        hold on;
        
    end
    
    print(H,'-depsc',filepath);
       
end


end

