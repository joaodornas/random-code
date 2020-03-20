
function reliabilityCC(date,site_index,channel,registro,video_index,start_time,end_time,bin_size,bandwidth,nConditions) 

tic

disp('BEGIN');

%%%   LOAD ALL TRIALS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load All Trials...');

Spass = load(strcat('_',registro,'-','v',int2str(video_index),'.mat'));

for i=1:nConditions
   
    trials_label(i).labels = find(Spass.stimIds == i); 
    
end

Spass.spike_times = ( Spass.spike_times./ 32000 ) .* 1000 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   RELIABILITY: CONDITION 1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reliability: CONDITION 1...');

trials_spikes_CRF = Spass.spike_times(trials_label(1).labels(:),:);

reliabilityFullMovieData.CRFxCRF.trials_spikes_CRF = trials_spikes_CRF;
reliabilityFullMovieData.CRFxCRF.schreiberData = schreiberCorr(trials_spikes_CRF,bandwidth,start_time,end_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   RELIABILITY: CONDITION 2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reliability: CONDITION 2...');

trials_spikes_NonCRF = Spass.spike_times(trials_label(2).labels(:),:);

reliabilityFullMovieData.NonCRFxNonCRF.trials_spikes_reverso = trials_spikes_NonCRF;
reliabilityFullMovieData.NonCRFxNonCRF.schreiberData = schreiberCorr(trials_spikes_NonCRF,bandwidth,start_time,end_time);

if nConditions == 3
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%   RELIABILITY: CONDITION 2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Reliability: CONDITION 2...');

    trials_spikes_ExtraNonCRF = Spass.spike_times(trials_label(2).labels(:),:);

    reliabilityFullMovieData.ExtraNonCRFxExtraNonCRF.trials_spikes_reverso = trials_spikes_ExtraNonCRF;
    reliabilityFullMovieData.ExtraNonCRFxExtraNonCRF.schreiberData = schreiberCorr(trials_spikes_ExtraNonCRF,bandwidth,start_time,end_time);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Save Data...');

filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index),'/');

mkdir(filepath,'reliability');

folder = 'reliability';

filepath = strcat(filepath,'/',folder,'/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel,'-','reliability-full-movie-bandwidth-',int2str(bandwidth));

save(filepath,'reliabilityFullMovieData');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END');

toc

end

