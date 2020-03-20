
function lifetimeSparseness(date,site_index,channel,registro,video_index,start_time,end_time,bin_size)

tic

disp('BEGIN');

%%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Spass...');

    Spass = load(strcat('_',registro,'-','v',int2str(video_index),'.mat'));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% LOAD CONDITIONS TRIALS LABELS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Conditions Trials Labels...');

    for i=1:3
   
        trials_label(i).label = find(Spass.stimIds == i); 
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%% SET RESOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Set Resolution...');

    spike_times = Spass.spike_times;
    
    spike_times = spike_times./ 32000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% BEGIN CONDITIONS LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Begin Conditions Loop...');

    for i=1:3
       
       disp(strcat('Begin Condition . ',int2str(i))); 
        
       trials_spikes = spike_times(trials_label(i).label,:);
        
       nTrials = size(trials_spikes,1);
       
       spikes_vector = reshape(trials_spikes.',[],1);

       spikes_vector = spikes_vector.';
       
       spikes_vector = sort(spikes_vector);
       
       spikes_vector = spikes_vector(spikes_vector>0);

       spikes_vector = spikes_vector(spikes_vector>(start_time/1000) & spikes_vector<(end_time/1000));
       
       nBins = (end_time - start_time)/bin_size;
       
       for k=1:nBins
          
           spikes = length(spikes_vector(spikes_vector>=((k-1)*bin_size/1000 + start_time/1000) & spikes_vector<(k*bin_size/1000 + start_time/1000)));
           
           sparseData(i).binRate(k) = spikes/(nTrials*bin_size/1000);

       end
  
       disp('...Calculate Sparseness');
       
       sparseData(i).lifetime_sparseness = kurtosis(sparseData(i).binRate) - 3;
       
       sparseData(i).trials_spikes = trials_spikes;
       sparseData(i).spikes_vector = spikes_vector;
       sparseData(i).nTrials = nTrials;
       
       disp(strcat('End Condition . ',int2str(i)));
       
    end
    
disp('End Conditions Loop');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Save data...');
    
    filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));
    filepath = strcat(filepath,'/','sparseness','/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel);
    filepath = strcat(filepath,'-bin_size-',int2str(bin_size),'-sparseness');
       
    save(filepath,'sparseData');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END');
    
toc
    
end