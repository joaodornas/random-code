function  singleEventAnalysis(registro,start_time,end_time,stimulus_start_time,stimulus_end_time,method,Visible,OS)

Spass = load(registro);

nConditions = 3;

PSTH_bin_size = 1;

frame_bin_size = 30;

latency_start_time = 500;

latency_end_time = latency_start_time + 300;

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

        latency_time = getLatencyForBack(protocol,condition);

        nTrials(j) = min(size(trials_spikes,1),size(trials_spikes,2));

        allSpikes = [];

        if (j == 2), disp('Begin invert spike train...'); end
        
        for k=1:nTrials(j)

            spikes = trials_spikes(k,:);

            a_e = spikes(spikes>=(start_time/1000) & spikes<(latency_start_time/1000));

            after_ae = spikes(spikes>=(latency_start_time/1000) & spikes<(end_time/1000)) - latency_time;

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

            Spike_Trains_Conditions(i).allTrials(k).trial = spikes;
            
            Spike_Trains_Conditions(i).allTrials(k).rate = numel(spikes(spikes>=(stimulus_start_time/1000) & spikes<=(stimulus_end_time/1000))) ./ ((stimulus_end_time - stimulus_start_time)/1000) ; 

            allSpikes = [allSpikes, spikes];

        end

        allSpikes = reshape(allSpikes.',[],1);

        allSpikes = sort(allSpikes);

        allSpikes = allSpikes(allSpikes>0);
        
        totalSpikes(i) = max(size(allSpikes,1),size(allSpikes,2));

        nBins_frames = (stimulus_end_time - stimulus_start_time) / frame_bin_size;

        for k=1:nBins_frames

            PSTH_Conditions_frames(i).rate(k) = numel(allSpikes(allSpikes>=((stimulus_start_time/1000) + (k-1)*(frame_bin_size/1000)) & allSpikes<=((stimulus_start_time/1000) + (k)*(frame_bin_size/1000)) ) ) / ( nTrials(j) * frame_bin_size/1000 );

            PSTH_Conditions_frames(i).spikes(k).times = allSpikes(allSpikes>=((stimulus_start_time/1000) + (k-1)*(frame_bin_size/1000)) & allSpikes<=((stimulus_start_time/1000) + (k)*(frame_bin_size/1000)) ).';

            PSTH_Conditions_frames(i).sd(k) = std(PSTH_Conditions_frames(i).spikes(k).times);

        end
        
        meanFiringRate(i) = mean(PSTH_Conditions_frames(i).rate);
        maxFiringRate(i) = max(PSTH_Conditions_frames(i).rate);
        stdFiringRate(i) = std(PSTH_Conditions_frames(i).rate);
       
end

for n=1:nTrials(3)
    
    AE_rate(n) = Spike_Trains_Conditions(3).allTrials(1,n).rate;
    
end

[AE_h,AE_p] = lillietest(AE_rate);

kernelMat = getKernel(registro,start_time,end_time,OS);

p = kernelMat.datakernel.p;

kernel = kernelMat.datakernel.kernel;

clear kernelMat;



for i=1:nConditions

    kernelMaxRate(i) = ( totalSpikes(i).*max(kernel(i).density( ( (stimulus_start_time/1000)*p ) : (stimulus_end_time/1000)*p )) ) / nTrials(i) ;

end


if strcmp(method,'per_spike')
    
    %%%% Find Events per distance between spikes
    
    for z=0:11
    
        spike_distance = 4 + z;
        
        limiar_spikes = 10;
    
        for i=5:30

            silence_distance = i*10;

            Events_Conditions = findEventsPerSpikeTime(kernel,Spike_Trains_Conditions,spike_distance,silence_distance,stimulus_start_time,stimulus_end_time,limiar_spikes,p);

            Correlations_Events_Spikes_Conditions = getAutoCorrSpikeEvents(Events_Conditions,Spike_Trains_Conditions,stimulus_start_time,frame_bin_size,p,'per_spike');

            Correlations_Events_PSTH_Conditions = getAutoCorrPSTHEvents(Events_Conditions,PSTH_Conditions_frames,frame_bin_size,p,stimulus_start_time,'per_spike');

            Fits_Gauss_PSTH_Conditions = getGaussianFitsAutoCorrEvents(Correlations_Events_PSTH_Conditions,'PSTH'); 
            
            Fits_Gauss_Spikes_Conditions = getGaussianFitsAutoCorrEvents(Correlations_Events_Spikes_Conditions,'spike');

            if strcmp(OS,'Win')

                mainPath = strcat('Z:\DATA\Forward-Backward\Alonso\Events\per-spike\');
                bar = '\';

            elseif strcmp(OS,'Mac')

                mainPath = strcat('/Volumes/Data/DATA/Forward-Backward/Alonso/Events/per-spike/');
                bar = '/';

            end

            mkdir(strcat(mainPath,'spike-',num2str(spike_distance),'-silence-',num2str(silence_distance)));

            name = registro(2:end-4);

            if strcmp(OS,'Win')

                folderpath = strcat(mainPath,'spike-',num2str(spike_distance),'-silence-',num2str(silence_distance),'\singleEvents\');
                bar = '\';

            elseif strcmp(OS,'Mac')

                folderpath = strcat(mainPath,'spike-',num2str(spike_distance),'-silence-',num2str(silence_distance),'/singleEvents/');
                bar = '/';

            end

            mkdir(folderpath,name);

            results = struct('name',name,'Spike_Trains_Conditions',Spike_Trains_Conditions,'PSTH_Conditions_frames',PSTH_Conditions_frames,'Events_Conditions',Events_Conditions,'Correlations_Events_Spikes_Conditions',Correlations_Events_Spikes_Conditions,'Correlations_Events_PSTH_Conditions',Correlations_Events_PSTH_Conditions,'Fits_Gauss_PSTH_Conditions',Fits_Gauss_PSTH_Conditions,'Fits_Gauss_Spikes_Conditions',Fits_Gauss_Spikes_Conditions,'kernel',kernel,'kernelMaxRate',kernelMaxRate,'maxFiringRate',maxFiringRate,'meanFiringRate',meanFiringRate,'AE_h',AE_h,'AE_p',AE_p,'AE_rate',AE_rate);

            save(strcat(folderpath,name,bar,name,'-','singleEventsAutoCorr-PSTH-Spike-autocorrelation-burst-per-spike-silence-distance-',int2str(i)),'results');

            plot_name = strcat(name,'-','spike-',num2str(spike_distance),'-silence-',num2str(silence_distance));

            plotFigureFitCorrelationEvent(Events_Conditions,Fits_Gauss_PSTH_Conditions,name,plot_name,'PSTH',folderpath,bar,Visible);
            plotFigureFitCorrelationEvent(Events_Conditions,Fits_Gauss_Spikes_Conditions,name,plot_name,'spike',folderpath,bar,Visible);

            checkEventsonRasterPlots(registro,start_time,end_time,results,folderpath,bar,'per_spike',p,Visible);

            clear Events_Conditions;
            clear Correlations_Events_Spikes_Conditions;
            clear Correlations_Events_PSTH_Conditions;
            clear Fits_Gauss_PSTH_Conditions;
            clear Fits_Gauss_Spikes_Conditions;

        end
    
    end

elseif strcmp(method,'threshold')

    %%%% Find Events per limiar on PSTH
    
    for i=1:10
        
        threshold = round(maxFiringRate(3) + i*stdFiringRate(3));        

        Events_Conditions = findEventsPerLimiar(PSTH_Conditions_frames,Spike_Trains_Conditions,threshold,stimulus_start_time,frame_bin_size);

        Correlations_Events_Spikes_Conditions = getAutoCorrSpikeEvents(Events_Conditions,Spike_Trains_Conditions,stimulus_start_time,frame_bin_size,p,'per_bin');

        Correlations_Events_PSTH_Conditions = getAutoCorrPSTHEvents(Events_Conditions,PSTH_Conditions_frames,frame_bin_size,p,stimulus_start_time,'per_bin');

        Fits_Gauss_PSTH_Conditions = getGaussianFitsAutoCorrEvents(Correlations_Events_PSTH_Conditions,'PSTH');

        Fits_Gauss_Spikes_Conditions = getGaussianFitsAutoCorrEvents(Correlations_Events_Spikes_Conditions,'spike');



        if strcmp(OS,'Win')

            mainPath = strcat('Z:\DATA\Forward-Backward\Alonso\Events\threshold-sd\');
            bar = '\';

        elseif strcmp(OS,'Mac')

            mainPath = strcat('/Volumes/Data/DATA/Forward-Backward/Alonso/Events/threshold-sd/');
            bar = '/';

        end

        mkdir(strcat(mainPath,'threshold-power-',num2str(i)));

        name = registro(2:end-4);

        if strcmp(OS,'Win')

            folderpath = strcat(mainPath,'threshold-power-',num2str(i),'\singleEvents\');
            bar = '\';

        elseif strcmp(OS,'Mac')

            folderpath = strcat(mainPath,'threshold-power-',num2str(i),'/singleEvents/');
            bar = '/';

        end

        mkdir(folderpath,name);

        results = struct('name',name,'Spike_Trains_Conditions',Spike_Trains_Conditions,'PSTH_Conditions_frames',PSTH_Conditions_frames,'Events_Conditions',Events_Conditions,'Correlations_Events_Spikes_Conditions',Correlations_Events_Spikes_Conditions,'Correlations_Events_PSTH_Conditions',Correlations_Events_PSTH_Conditions,'Fits_Gauss_PSTH_Conditions',Fits_Gauss_PSTH_Conditions,'Fits_Gauss_Spikes_Conditions',Fits_Gauss_Spikes_Conditions,'kernel',kernel,'kernelMaxRate',kernelMaxRate,'maxFiringRate',maxFiringRate,'meanFiringRate',meanFiringRate,'AE_h',AE_h,'AE_p',AE_p,'AE_rate',AE_rate);

        save(strcat(folderpath,name,bar,name,'-','singleEventsAutoCorr-PSTH-Spike-autocorrelation-burst-per-threshold-',int2str(i)),'results');

        plot_name = strcat(name,'-','threshold-',int2str(i));
        
        plotFigureFitCorrelationEvent(Events_Conditions,Fits_Gauss_PSTH_Conditions,name,plot_name,'PSTH',folderpath,bar,Visible);
        plotFigureFitCorrelationEvent(Events_Conditions,Fits_Gauss_Spikes_Conditions,name,plot_name,'spike',folderpath,bar,Visible);

        checkEventsonRasterPlots(registro,start_time,end_time,results,folderpath,bar,'per_bin',p,Visible);

        clear Events_Conditions;
        clear Correlations_Events_Spikes_Conditions;
        clear Correlations_Events_PSTH_Conditions;
        clear Fits_Gauss_PSTH_Conditions;
        clear Fits_Gauss_Spikes_Conditions;
    
    end
        
end

end






