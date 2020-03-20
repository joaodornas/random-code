function Correlations_PSTH_Conditions = getAutoCorrPSTH(PSTH_Conditions,kind)

    condicoes = max(size(PSTH_Conditions,1),size(PSTH_Conditions,2));
    
    for c=1:condicoes

        PSTH_segment = PSTH_Conditions(c).rate - mean(PSTH_Conditions(c).rate)^2;
        
        xcorrelation = xcorr(PSTH_segment,'biased');
        
        sd = std(PSTH_segment);
        
        xcorrelation = xcorrelation ./ sd.^2;
        
        convolution = getCorrelationFFT(PSTH_segment,PSTH_segment);
        
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
        
        eval(strcat('Correlations_PSTH_Conditions(c).xcorr_',kind,' = new_xcorrelation;'));

        eval(strcat('Correlations_PSTH_Conditions(c).correlogram_',kind,' = new_convolution;'));
            
    end


end