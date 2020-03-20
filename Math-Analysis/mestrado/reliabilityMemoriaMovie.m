
function reliabilityMemoriaMovie(date,site_index,channel,registro,video_index,start_time,end_time,bandwidth) 

tic

disp('BEGIN');

%%%   LOAD ALL TRIALS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load All Trials...');

Spass = load(strcat('_',registro,'-','v',int2str(video_index),'.mat'));

Spass.spike_times = ( Spass.spike_times./ 32000 ) .* 1000 ;

for i=1:3
   
    trials_label(i).label = find(Spass.stimIds == i); 
    
end

disp('...define video pairs by history');

for i=1:3
    
    for j=1:3
        
        video_pairs(i,j).trials = zeros(1);

    end
    
end

for i=1:3
   
    for j=1:length(trials_label(i).label)
    
        if trials_label(i).label(j) ~= 1

            video_pairs(i,(Spass.stimIds(trials_label(i).label(j) - 1)) ).trials = [video_pairs(i,Spass.stimIds(trials_label(i).label(j) - 1)).trials, trials_label(i).label(j)]; 
        
        end
        
    end
    
end


for i=1:3
    
    for j=1:3
        
        video_pairs(i,j).trials(video_pairs(i,j).trials==0) = [];

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   READ THREE KINDS OF TRIALS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Read Three Kinds of Trials...');

trials_spikes_direto = Spass.spike_times(trials_label(1).label,:);
trials_spikes_ExpoAfterDireto = Spass.spike_times(video_pairs(3,1).trials,:);
trials_spikes_ExpoAfterReverso = Spass.spike_times(video_pairs(3,2).trials,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% RELIABILITY
disp('Calculate Reliability...');

trials_direto_expoAdireto = [trials_spikes_direto; trials_spikes_ExpoAfterDireto];
trials_direto_expoAreverso = [trials_spikes_direto; trials_spikes_ExpoAfterReverso];

reliabilityMemoriaMovieData.schreiberData.DiretoExpoAfterDireto = schreiberCorr(trials_direto_expoAdireto,bandwidth,start_time,end_time);
reliabilityMemoriaMovieData.schreiberData.DiretoExpoAfterReverso = schreiberCorr(trials_direto_expoAreverso,bandwidth,start_time,end_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Save Data...');

reliabilityMemoriaMovieData.trials_spikes_direto = trials_spikes_direto;
reliabilityMemoriaMovieData.trials_spikes_ExpoAfterDireto = trials_spikes_ExpoAfterDireto;
reliabilityMemoriaMovieData.trials_spikes_ExpoAfterReverso = trials_spikes_ExpoAfterReverso;

filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));

folder = 'memoria-movie';

filepath = strcat(filepath,'/',folder,'/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel,'-','reliability-memoria-movie');

save(filepath,'reliabilityMemoriaMovieData');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END');

toc

end

