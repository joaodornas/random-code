function Correlations_Events_PSTH_Conditions = getAutoCorrPSTHEvents(Events_Conditions,PSTH_Conditions,frame_bin_size,p,stimulus_start_time,per_bin_per_spike)

    condicoes = max(size(PSTH_Conditions,1),size(PSTH_Conditions,2));
    
    for c=1:condicoes
        
        nEvents = max(size(Events_Conditions(c).eventos,1),size(Events_Conditions(c).eventos,2));

        if nEvents >= 1
    
            for n=1:nEvents
                
                if strcmp(per_bin_per_spike,'per_bin')
                    
                     start_ = Events_Conditions(c).eventos(n).start;

                     end_ = Events_Conditions(c).eventos(n).end;
                    
                elseif strcmp(per_bin_per_spike,'per_spike')
                    
                     start_ =  (Events_Conditions(c).eventos(n).start - stimulus_start_time) / (frame_bin_size)  ;

                     end_ =  (Events_Conditions(c).eventos(n).end - stimulus_start_time) / (frame_bin_size)  ;
                     
                     
                     if round(start_) >  start_
                         
                         start_ = round(start_) - 1;
                         
                     elseif start_ <= 1
                         
                         start_ = 1;
                         
                     else
                         
                         start_ = round(start_);
                         
                     end
                     
                     if start_ == 0, start_ = 1; end
                     
                     if round(end_) > end_
                         
                         end_ = round(end_) - 1;
                         
                     elseif end_ <= 1
                         
                         end_ = 1;
                         
                     else
                         
                         end_ = round(end_);
                         
                     end
                     
                     if end_ == 0, end_ = 1; end
                    
                end
                
                
                PSTH_segment = PSTH_Conditions(c).rate(start_:end_) - mean(PSTH_Conditions(c).rate(start_:end_))^2;

                sd = std(PSTH_segment);

                xcorrelation = xcorr(PSTH_segment,'biased');

                xcorrelation = xcorrelation ./ sd.^2;

                convolution = getCorrelation(PSTH_segment,PSTH_segment);

                if length(xcorrelation) > (600 + 1)

                    xcorrelation(:) = xcorrelation((length(xcorrelation)/4):(length(xcorrelation)/2));

                end

                if length(convolution) > (600 + 1)

                    convolution(:) = convolution((length(convolution)/4):(length(convolution)/2));

                end

                Correlations_Events_PSTH_Conditions(c).eventos(n).xcorr_PSTH = xcorrelation;

                Correlations_Events_PSTH_Conditions(c).eventos(n).correlogram_PSTH = convolution;

            end

         else

            Correlations_Events_PSTH_Conditions(c).eventos = [];

         end
            
    end


end