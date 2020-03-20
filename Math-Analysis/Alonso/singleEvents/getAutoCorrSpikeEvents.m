function Correlations_Events_Spikes_Conditions = getAutoCorrSpikeEvents(Events_Conditions,Spike_Trains_Conditions,start_time,frame_bin_size,p,per_bin_per_spike)

    condicoes = max(size(Spike_Trains_Conditions,1),size(Spike_Trains_Conditions,2));
    
    for c=1:condicoes

        nEvents = max(size(Events_Conditions(c).eventos,1),size(Events_Conditions(c).eventos,2));
        
        if nEvents >= 1

            for n=1:nEvents

                start_ = Events_Conditions(c).eventos(n).start;

                end_ = Events_Conditions(c).eventos(n).end;

                nTrials = max(size(Spike_Trains_Conditions(c).allTrials,1),size(Spike_Trains_Conditions(c).allTrials,2));

                spikes = [];
                for k=1:nTrials

                    spikes = [spikes, Spike_Trains_Conditions(c).allTrials(1,k).trial(:).'];

                end

                spikes = sort(spikes);

                if strcmp(per_bin_per_spike,'per_bin')
                    
                    spikes_segment = spikes(spikes>=((start_time/1000) + (start_-1)*(frame_bin_size/1000)) & spikes<=((start_time/1000) + end_*(frame_bin_size/1000)));

                elseif strcmp(per_bin_per_spike,'per_spike')
                   
                    spikes_segment = spikes(spikes>=( (start_/p) ) & spikes<=( (end_/p) ));

                end
                
                Filter_Dimension = 0;

                HSIZE = 0;

                nTrials_vector = 1;

                sigma = 0;

                if strcmp(per_bin_per_spike,'per_bin')
                    
                    start_time_segment = ((start_time) + (start_-1)*(frame_bin_size));

                    end_time_segment = ((start_time) + end_*(frame_bin_size));
                    
                elseif strcmp(per_bin_per_spike,'per_spike')
                    
                    start_time_segment = start_  ;

                    end_time_segment =  end_ ;
                    
                end

                trials_spikes(1,:) = spikes_segment(:);

                spikes_vector = getTime_vector(trials_spikes,nTrials_vector,HSIZE,sigma,start_time_segment,end_time_segment,Filter_Dimension);

                clear trials_spikes;

                spikes_vector = spikes_vector.time_vector;

                first_spike = getFirstSpike(spikes_vector);

                last_spike = getLastSpike(spikes_vector);

                spikes_vector = spikes_vector(first_spike:last_spike);
                
                %spikes_segment = spikes_vector - mean(spikes_vector)^2;
                
                spikes_segment = spikes_vector;

                sd = std(spikes_segment);

                xcorrelation = xcorr(spikes_segment,spikes_segment,'biased');

                xcorrelation = xcorrelation ./ sd.^2;

                convolution = getCorrelation(spikes_segment,spikes_segment);

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

                Correlations_Events_Spikes_Conditions(c).eventos(n).xcorr_spike = new_xcorrelation;

                Correlations_Events_Spikes_Conditions(c).eventos(n).correlogram_spike = new_convolution;

            end

        else

            Correlations_Events_Spikes_Conditions(c).eventos = [];

        end
        
        
    end


end