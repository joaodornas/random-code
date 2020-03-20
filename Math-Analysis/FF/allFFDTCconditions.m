function allFFDTCconditions

tic

registro = importdata('memoryBackwardProtocols.txt');

for r=1:length(registro)
    
   [start_time, end_time, blk] = DTC_forback(char(registro{r}));
      
   if strcmp(blk,'none')
        
        Spass = NaN;
        
   else
        
        Spass = load(char(blk));
        
        getFF(Spass,start_time,end_time);
    
   end

end


function getFF(Spass,start_time,end_time)
        
    nConditions = 16;

    bins = 1000;

            
            %%% READ TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %disp('Read trials...');

            name = char(Spass.cellname);

            %start_time = char(Spass.parameters.baseline_duration);

            %end_time = str2double(char(Spass.parameters.baseline_duration))  + start_time;

            spike_times = Spass.spike_times./ 32000;
            
            time_bins = 1:bins;
            
            for i=1:bins

                size_of_bins(i) = (end_time - start_time)/i;

            end

            for n=1:nConditions
                
                labels = find(Spass.stimIds == n);
            
                allTrials = spike_times(Spass.stimIds(labels),:);
                
                nTrials = size(allTrials,1);

                TBallTrials = TrialsFF(allTrials,start_time,end_time,time_bins);

                clear allTrials;
                
                resultsDTC_across_trials = FF_across_trials(TBallTrials,start_time,end_time,time_bins,nTrials);
                
                FFDTC_across_trials = struct('name',name,'resultsDTC_across_trials',resultsDTC_across_trials,'size_of_bins',size_of_bins);
                
                save(strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF-DTC-conditions/',name,'-FF-across_trials-ForBack-DTC-conditions-',int2str(n),'.mat'),'FFDTC_across_trials');

                clear resultsDTC_across_trials;
                clear FFDTC_across_trials;
                
                resultsDTC_across_bins = FF_across_bins(TBallTrials,start_time,end_time,time_bins,nTrials);

                FFDTC_across_bins = struct('name',name,'resultsDTC_across_bins',resultsDTC_across_bins,'size_of_bins',size_of_bins);

                save(strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF-DTC-conditions/',name,'-FF-across_bins-ForBack-DTC-conditions-',int2str(n),'.mat'),'FFDTC_across_bins');

                TB_optimal_DTC = TrialsFFoptimal(TBallTrials);

                save(strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF-DTC-conditions/',name,'-FF-TB-ForBack-DTC-conditions-',int2str(n),'.mat'),'TB_optimal_DTC');

                clear resultsDTC_across_bins;
                clear FFDTC_across_bins;

                clear TBallTrials;
                clear TB_optimal_DTC;

            end
            
            clear spike_times;

            Spass.spike_times = [];
            Spass.stimIds = [];        

    
end

toc

end

