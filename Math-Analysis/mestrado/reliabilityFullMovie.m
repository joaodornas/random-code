
function reliabilityFullMovie(date,site_index,channel,registro,video_index,start_time,end_time,bin_size,bandwidth) 

tic

disp('BEGIN');

%%%   LOAD ALL TRIALS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load All Trials...');

Spass = load(strcat('_',registro,'-','v',int2str(video_index),'.mat'));

for i=1:3
   
    trials_label(i).labels = find(Spass.stimIds == i); 
    
end

Spass.spike_times = ( Spass.spike_times./ 32000 ) .* 1000 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   RELIABILITY: DIRETO - DIRETO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reliability: Direto - Direto...');

trials_spikes_direto = Spass.spike_times(trials_label(1).labels(:),:);

reliabilityFullMovieData.DiretoDireto.trials_spikes_direto = trials_spikes_direto;
reliabilityFullMovieData.DiretoDireto.schreiberData = schreiberCorr(trials_spikes_direto,bandwidth,start_time,end_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   RELIABILITY: REVERSO - REVERSO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reliability: Reverso - Reverso...');

trials_spikes_reverso = Spass.spike_times(trials_label(2).labels(:),:);

reliabilityFullMovieData.ReversoReverso.trials_spikes_reverso = trials_spikes_reverso;
reliabilityFullMovieData.ReversoReverso.schreiberData = schreiberCorr(trials_spikes_reverso,bandwidth,start_time,end_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   RELIABILITY: INVERTIDO - INVERTIDO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reliability: Invertido - Invertido...');

trials_spikes_reverso = Spass.spike_times(trials_label(2).labels(:),:);

nTrials = size(trials_spikes_reverso,1);

nBins = (end_time - start_time)/bin_size;

    for k=1:nTrials
    
        spikes = trials_spikes_reverso(k,:);
        
        spikes = spikes(spikes>0);
        
        spikes = spikes(spikes>=(start_time) & spikes<=(end_time));

        spikes_reversed = invertSpikeTrainmili(spikes,start_time,nBins,bin_size);
        
        trials(k).spikes_reversed = sort(spikes_reversed.spikes);
        
        nSpikes(k) = length(spikes_reversed.spikes);

    end

maxSize = max(nSpikes);

    for k=1:nTrials
       
        trials_spikes_invertido(k,:) = [trials(k).spikes_reversed zeros(1,(maxSize - length(trials(k).spikes_reversed)))];
        
    end

reliabilityFullMovieData.InvertidoInvertido.trials_spikes_reverso = trials_spikes_reverso;
reliabilityFullMovieData.InvertidoInvertido.trials_spikes_invertido = trials_spikes_invertido;
reliabilityFullMovieData.InvertidoInvertido.scheiberData = schreiberCorr(trials_spikes_invertido,bandwidth,start_time,end_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   RELIABILITY: DIRETO - REVERSO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reliability: Direto - Reverso...');

trials_spikes_DiretoReverso = [trials_spikes_direto; trials_spikes_reverso];

reliabilityFullMovieData.DiretoReverso.trials_spikes_DiretoReverso = trials_spikes_DiretoReverso;
reliabilityFullMovieData.DiretoReverso.schreiberData = schreiberCorr(trials_spikes_DiretoReverso,bandwidth,start_time,end_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   RELIABILITY: DIRETO - INVERTIDO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reliability: Direto - Invertido...');

dif = size(trials_spikes_direto,2) - size(trials_spikes_invertido,2);

trials_spikes_invertido = [trials_spikes_invertido zeros(size(trials_spikes_invertido,1),dif)];

trials_spikes_DiretoInvertido = [trials_spikes_direto; trials_spikes_invertido];

reliabilityFullMovieData.DiretoInvertido.trials_spikes_DiretoInvertido = trials_spikes_DiretoInvertido;
reliabilityFullMovieData.DiretoInvertido.schreiberData = schreiberCorr(trials_spikes_DiretoInvertido,bandwidth,start_time,end_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Save Data...');

filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));

folder = 'full-movie/reliability';

filepath = strcat(filepath,'/',folder,'/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel,'-','reliability-full-movie-bandwidth-',int2str(bandwidth));

save(filepath,'reliabilityFullMovieData');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END');

toc

end

