function Correlations_Spikes_Conditions = getAutoCorrSpike(Spike_Trains_Conditions,start_time,end_time)

    condicoes = max(size(Spike_Trains_Conditions,1),size(Spike_Trains_Conditions,2));
    
    for c=1:condicoes

        nTrials = max(size(Spike_Trains_Conditions(c).allTrials,1),size(Spike_Trains_Conditions(c).allTrials,2));
        
        spikes = [];
        
        for n=1:nTrials
            
            spikes = [spikes, ( (Spike_Trains_Conditions(c).allTrials(1,n).trial(:).') + ( (n-1) * (end_time - start_time)/1000 ) )];
            
        end
        
        spikes = sort(spikes);
        
        Filter_Dimension = 0;

        HSIZE = 0;

        nTrials_vector = 1;

        sigma = 0;

        start_time_segment = start_time;

        end_time_segment = (end_time - start_time)*nTrials ;

        trials_spikes(1,:) = spikes(:);

        spikes_vector = getTime_vector(trials_spikes,nTrials_vector,HSIZE,sigma,start_time_segment,end_time_segment,Filter_Dimension);

        clear trials_spikes;

        spikes_vector = spikes_vector.time_vector;

        first_spike = getFirstSpike(spikes_vector);

        last_spike = getLastSpike(spikes_vector);

        spikes_vector = spikes_vector(first_spike:last_spike);
        
        xcorrelation = xcorr(spikes_vector,spikes_vector);
        
        sd = std(spikes_vector);
        
        xcorrelation = xcorrelation ./ sd.^2;
        
        convolution = getCorrelationFFT(spikes_vector,spikes_vector);
        
        if length(xcorrelation) > (600 + 1)

            new_xcorrelation = xcorrelation((length(xcorrelation)/4):(length(xcorrelation)/2));

        else
            
            new_xcorrelation = xcorrelation;
            
        end


        if length(convolution) > (600 + 1)

            new_convolution = convolution((length(convolution)/4):(length(convolution)/2));

        else
            
            new_convolution = convolution;
            
        end

        Correlations_Spikes_Conditions(c).xcorr_Spike = new_xcorrelation;

        Correlations_Spikes_Conditions(c).correlogram_Spike = new_convolution;
        
        
    end


end