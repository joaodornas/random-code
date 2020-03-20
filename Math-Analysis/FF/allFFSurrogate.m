function allFFSurrogate

disp('allFFSurrogateBS');

tic

registro{1} = 'allTrials-surrogate-rate-2.mat';
registro{2} = 'allTrials-surrogate-rate-3.mat';
registro{3} = 'allTrials-surrogate-rate-4.mat';
registro{4} = 'allTrials-surrogate-rate-5.mat';
registro{5} = 'allTrials-surrogate-rate-6.mat';

for s=1:length(registro)
    
    Spass(s) = load(char(registro{s}));
    
end

getFF(Spass);

function getFF(Spass)
        
    nConditions = 1;

    %bins = 1000;

    for h=1:length(Spass)

        %%% READ TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %disp('Read trials...');
        
        name = strcat('surrogate-rate-',int2str(h+6));
        
        start_time = 0;

        end_time = 10000;
    
        bins = end_time;
        
        allTrials = Spass(h).allTrials;

        time_bins = 1:bins;
        
        for i=1:bins
            
            size_of_bins(i) = (end_time - start_time)/i;
            
        end
        
        nTrials = size(allTrials,1);
        
        TBallTrials = TrialsFF(allTrials,start_time,end_time,time_bins);
        
        clear allTrials;
        
        results_across_trials = FF_across_trials(TBallTrials,start_time,end_time,time_bins,nTrials);
        
        FF_across_trials_data = struct('name',name,'results_across_trials',results_across_trials,'size_of_bins',size_of_bins);
        
        save(strcat('/Volumes/Data/DATA/Surrogate/FF/',name,'-FF-across_trials.mat'),'FF_across_trials_data','-v7.3');  
        
        clear results_across_trials;
        clear FF_across_trials;
        
        results_across_bins = FF_across_bins(TBallTrials,start_time,end_time,time_bins,nTrials);
        
        FF_across_bins_data = struct('name',name,'results_across_bins',results_across_bins,'size_of_bins',size_of_bins);
        
        save(strcat('/Volumes/Data/DATA/Surrogate/FF/',name,'-FF-across_bins.mat'),'FF_across_bins_data','-v7.3');
        
        TB_optimal = TrialsFFoptimal(TBallTrials);
        
        save(strcat('/Volumes/Data/DATA/Surrogate/FF/',name,'-FF-TB.mat'),'TB_optimal','-v7.3');
        
        clear results_across_bins;
        clear FF_across_bins;
        
        clear nTrials;
        clear TBallTrials;
        clear TB_optimal;
        
        clear time_bins;
        clear name;
        
    end
    
end

toc

end

