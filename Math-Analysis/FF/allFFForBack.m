function allFFForBack

tic

registro = importdata('memoryBackwardProtocols.txt');

for r=1:length(registro)    
    
    Spass(r) = load(char(registro{r}));

end

getFF(Spass);

function getFF(Spass)
        
    nConditions = 3;

    start_time = 500;

    end_time = 9500;
    
    bins = 1000;

    for h=1:length(Spass)

        %%% READ TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %disp('Read trials...');

        name = char(Spass(h).cellname);

        spike_times = Spass(h).spike_times./ 32000;
        
        forwardlabels = find(Spass(h).stimIds == 1);

        forwardtrials = spike_times(forwardlabels(:),:);

        backwardlabels = find(Spass(h).stimIds == 2);
        
        backwardtrials = spike_times(backwardlabels(:),:);
        
        AE_labels = find(Spass(h).stimIds == 3);

        AE_trials = spike_times(AE_labels(:),:);
 
        clear spike_times;
        clear forwardlabels;
        clear backwardlabels;
        clear AE_labels;
        
        Spass(h).spike_times = [];
        Spass(h).stimIds = [];

        time_bins = 1:bins;
        
        for i=1:bins
            
            size_of_bins(i) = (end_time - start_time)/i;
            
        end
        
        nTrialsForward = size(forwardtrials,1);
        
        nTrialsBackward = size(backwardtrials,1);
        
        nTrials_AE = size(AE_trials,1);
        
        TBforward = TrialsFF(forwardtrials,start_time,end_time,time_bins);
        
        resultsForward_across_trials = FF_across_trials(TBforward,start_time,end_time,time_bins,nTrialsForward);
        
        FFForward_across_trials = struct('name',name,'resultsForward_across_trials',resultsForward_across_trials,'size_of_bins',size_of_bins);
        
        save(strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF/',name,'-FF-across_trials-Forward.mat'),'FFForward_across_trials');  
        
        clear resultsForward_across_trials;
        clear FFForward_across_trials;
        
        resultsForward_across_bins = FF_across_bins(TBforward,start_time,end_time,time_bins,nTrialsForward);
        
        FFForward_across_bins = struct('name',name,'resultsForward_across_bins',resultsForward_across_bins,'size_of_bins',size_of_bins);
        
        save(strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF/',name,'-FF-across_bins-Forward.mat'),'FFForward_across_bins');
        
        TB_optimal_Forward = TrialsFFoptimal(TBforward);
        
        save(strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF/',name,'-FF-TB-Forward.mat'),'TB_optimal_Forward');
        
        clear resultsForward_across_bins;
        clear FFForward_across_bins;
        
        clear nTrialsForward;
        clear TBforward;
        clear TB_optimal_Forward;
        clear forwardtrials;
        
        TBbackward = TrialsFF(backwardtrials,start_time,end_time,time_bins);

        resultsBackward_across_trials = FF_across_trials(TBbackward,start_time,end_time,time_bins,nTrialsBackward);
        
        FFBackward_across_trials = struct('name',name,'resultsBackward_across_trials',resultsBackward_across_trials,'size_of_bins',size_of_bins);
        
        save(strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF/',name,'-FF-across_trials-Backward.mat'),'FFBackward_across_trials');
        
        clear resultsBackward_across_trials;
        clear FFBackward_across_trials;
        
        resultsBackward_across_bins = FF_across_bins(TBbackward,start_time,end_time,time_bins,nTrialsBackward);

        FFBackward_across_bins = struct('name',name,'resultsBackward_across_bins',resultsBackward_across_bins,'size_of_bins',size_of_bins);
        
        save(strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF/',name,'-FF-across_bins-Backward.mat'),'FFBackward_across_bins');
        
        TB_optimal_Backward = TrialsFFoptimal(TBbackward);
        
        save(strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF/',name,'-FF-TB-Backward.mat'),'TB_optimal_Backward');
        
        clear resultsBackward_across_bins;
        clear FFBackward_across_bins;
        
        clear nTrialsBackward;
        clear TBbackward;
        clear TB_optimal_Backward;
        clear backwardtrials;
        
        TB_AE = TrialsFF(AE_trials,start_time,end_time,time_bins);
        
        results_AE_across_trials = FF_across_trials(TB_AE,start_time,end_time,time_bins,nTrials_AE);
        
        FF_AE_across_trials = struct('name',name,'results_AE_across_trials',results_AE_across_trials,'size_of_bins',size_of_bins);
        
        save(strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF/',name,'-FF-across_trials-AE.mat'),'FF_AE_across_trials');  
        
        clear results_AE_across_trials;
        clear FF_AE_across_trials;
        
        results_AE_across_bins = FF_across_bins(TB_AE,start_time,end_time,time_bins,nTrials_AE);
        
        FF_AE_across_bins = struct('name',name,'results_AE_across_bins',results_AE_across_bins,'size_of_bins',size_of_bins);
        
        save(strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF/',name,'-FF-across_bins-AE.mat'),'FF_AE_across_bins');
        
        TB_optimal_AE = TrialsFFoptimal(TB_AE);
        
        save(strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF/',name,'-FF-TB-AE.mat'),'TB_optimal_AE');
        
        clear results_AE_across_bins;
        clear FF_AE_across_bins;
        
        clear nTrials_AE;
        clear TB_AE;
        clear TB_optimal_AE;
        clear AE_trials;
        
        clear time_bins;
        clear name;
        
        Spass(h).psth = [];
        Spass(h).parameters = [];
        
    end
    
end

toc

end

