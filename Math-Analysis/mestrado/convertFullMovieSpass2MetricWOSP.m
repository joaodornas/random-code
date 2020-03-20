function X = convertFullMovieSpass2MetricWOSP(date,site_index,channel,registro,video_index,start_time,end_time,invert_movie,bin_size,nRepetitions,nConditions,protocolo)

tic

disp('BEGIN');

%%%   DEFINE CATEGORY LABEL   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Define category label...');

if video_index > 100
    
    vn = round(video_index/100);
    
else
    
    vn = video_index;
    
end

if protocolo == 1
    
    switch vn

        case 1

            category{1} = 'Video 1: Swans. CENTRO';
            category{2} = 'Video 1: Swans. CONTORNO';
            category{3} = 'Video 1: Swans. EXTRA-CONTORNO';

        case 2

            category{1} = 'Video 2: Mountains. CENTRO';
            category{2} = 'Video 2: Mountains. CONTORNO';
            category{3} = 'Video 2: Mountains. EXTRA-CONTORNO';

        case 3

            category{1} = 'Video 3: Lizard. CENTRO';
            category{2} = 'Video 3: Lizard. CONTORNO';
            category{3} = 'Video 3: Lizard. EXTRA-CONTORNO';

        case 4

            category{1} = 'Video 4: Beetle. CENTRO';
            category{2} = 'Video 4: Beetle. CONTORNO';
            category{3} = 'Video 4: Beetle. EXTRA-CONTORNO';

        case 5 

            category{1} = 'Video 5: Chimp Climbing. CENTRO';
            category{2} = 'Video 5: Chimp Climbing. CONTORNO';
            category{3} = 'Video 5: Chimp Climbing. EXTRA-CONTORNO';

        case 6

            category{1} = 'Video 6: Chimp Hitting. CENTRO';
            category{2} = 'Video 6: Chimp Hitting. CONTORNO';
            category{3} = 'Video 6: Chimp Hitting. EXTRA-CONTORNO';

        case 7

            category{1} = 'Video 7: Panter. DIRECT';
            category{2} = 'Video 7: Panter. REVERSE';    

        case 8 

            category{1} = 'Video 8: Panda. CENTRO';
            category{2} = 'Video 8: Panda. CONTORNO'; 
            category{3} = 'Video 8: Panda. EXTRA-CONTORNO';

        case 9

            category{1} = 'Video 9: Crocodile 1. CENTRO';
            category{2} = 'Video 9: Crocodile 1. CONTORNO'; 
            category{3} = 'Video 9: Crocodile 1. EXTRA-CONTORNO';

        case 10

            category{1} = 'Video 10: Crocodile 2. CENTRO';
            category{2} = 'Video 10: Crocodile 2. CONTORNO'; 
            category{3} = 'Video 10: Crocodile 2. EXTRA-CONTORNO';

        case 11

            category{1} = 'Video 11: Birds. CENTRO';
            category{2} = 'Video 11: Birds. CONTORNO'; 
            category{3} = 'Video 11: Birds. EXTRA-CONTORNO';

    end
    
elseif protocolo == 2
     
    switch vn

        case 1

            category{1} = 'Video 1: Swans. DIRETO';
            category{2} = 'Video 1: Swans. REVERTO';
            category{3} = 'Atividade Espontanea';

        case 2

            category{1} = 'Video 2: Mountains. DIRETO';
            category{2} = 'Video 2: Mountains. REVERSO';
            category{3} = 'Atividade Espontanea';

        case 3

            category{1} = 'Video 3: Lizard. DIRETO';
            category{2} = 'Video 3: Lizard. REVERSO';
            category{3} = 'Atividade Espontanea';

        case 4

            category{1} = 'Video 4: Beetle. DIRETO';
            category{2} = 'Video 4: Beetle. REVERSO';
            category{3} = 'Atividade Espontanea';
            
     end
     
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   DEFINE FILE NAME AND PATH   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Define file name and path...');

filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));
%filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));


if (invert_movie == 1)

    folder = 'full-movie/full-movie-invertido';
    file_type = strcat('full-movie-invertido-bin-size-',int2str(bin_size));
    category{2} = strcat(category{2},'-invertido');

else
    
    folder = 'full-movie/full-movie-original';
    file_type = 'full-movie-original';
    
end

filepath = strcat(filepath,'/',folder,'/','v',int2str(video_index),'-','WOSP-',date,'-','sitio',int2str(site_index),'-',channel,'-',file_type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   SET METRIC SPACE PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Set Metric Space parameters...');

%Escala de tempo em segundos
time_scale = 1;

%Escala da resolu??o temporal em segundos
time_resolution = 0.0001;

%Sistema Internacional de Medidas
si_prefix = 1;

%N?mero de classes de est?mulos
M = nConditions;

%N?mero de s?tios (canais)
N = 1;

%N?mero de trials (repeti??es por condi??o)
P = nRepetitions;

%Tempos inicial e final do trial
%start_time = 0;
%end_time = 0;

%N?mero de pontos no vetor list
%Q = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   LOAD TRIALS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load trials...');
    
Spass = load(strcat('_',registro,'-','v',int2str(video_index),'.mat'));

for i=1:M
   
    trials_label(i).labels = find(Spass.stimIds == i); 
    
end

Spass.spike_times = Spass.spike_times./ 32000;

categories(1:M) = struct('label','','P',P,'trials',zeros(P));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   BEGIN CATEGORIES LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Begin categories loop...');

for i=1:M
    
    disp(strcat('Category .',int2str(i)));
    
    disp('...get trials');
    
    trials_spikes = Spass.spike_times(trials_label(i).labels(:),:);
    
    nTrials = size(trials_spikes,1);
    
    trials(N,1:nTrials) = struct('start_time',start_time/1000,'end_time',end_time/1000,'Q',0,'list',zeros(size(trials_spikes,2)));
    
    for k=1:nTrials
    
        spikes = trials_spikes(k,:);
        spikes = spikes(spikes>0);
        
        spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));
        spikes = sort(spikes);
           
        if (i == 2) & (invert_movie == 1)
           
            nBins = (end_time - start_time)/bin_size;
            
            disp(strcat('...invert trial .',int2str(k)));
            
            spikes_inverted = invertSpikeTrain(spikes,start_time,nBins,bin_size);
            spikes = spikes_inverted.spikes;
            
        end
        
        disp(strcat('...set trial .',int2str(k)));
        
        trials(N,k).Q = size(spikes,2);
        trials(N,k).list = spikes;
        
    end
    
    trials = trials';
    
    categories(i).label = category{i};
    categories(i).P = nTrials;
    categories(i).trials = trials;
    
    clear trials;
        
end

disp('End categories loop...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CREATE METRIC SPACE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Create Metric Space data...');

categories = categories';

sites = struct('label','unit_001','recording_tag','episodic','time_scale',time_scale,'time_resolution',time_resolution,'si_unit','none','si_prefix',si_prefix);

X = struct('M',M,'N',N,'sites',sites,'categories',categories);

%X = struct('M',M,'N',N,'sites',struct('label','unit_001','recording_tag','episodic','time_scale',time_scale,'time_resolution',time_resolution,'si_unit','none','si_prefix',si_prefix),'categories',struct('label','','P',P,'trials',struct('start_time',start_time,'end_time',end_time,'Q',Q,'list',list)));

%stawrite(X,strcat(filepath,'-1'));
stawrite(X,filepath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END');

toc

end