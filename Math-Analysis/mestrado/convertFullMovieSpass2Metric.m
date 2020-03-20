function X = convertFullMovieSpass2Metric(date,site_index,channel,registro,video_index,start_time,end_time,invert_movie,bin_size,nRepetitions,ae)

tic

disp('BEGIN');

%%%   DEFINE CATEGORY LABEL   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Define category label...');

if video_index > 100
    
    vn = round(video_index/100);
    
else
    
    vn = video_index;
    
end

switch vn
    
    case 1
        
        category{1} = 'Video 1: Swans. DIRECT';
        category{2} = 'Video 1: Swans. REVERSE';
        category{3} = 'Video 1: Spontaneous Activity';
                
    case 2
        
        category{1} = 'Video 2: Mountains. DIRECT';
        category{2} = 'Video 2: Mountains. REVERSE';
        category{3} = 'Video 2: Spontaneous Activity';
        
    case 3
        
        category{1} = 'Video 3: Lizard. DIRECT';
        category{2} = 'Video 3: Lizard. REVERSE';
        category{3} = 'Video 3: Spontaneous Activity';
        
    case 4
        
        category{1} = 'Video 4: Beetle. DIRECT';
        category{2} = 'Video 4: Beetle. REVERSE';
        category{3} = 'Video 4: Spontaneous Activity';
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   DEFINE FILE NAME AND PATH   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Define file name and path...');

filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));

if ae == 1
        
        mkdir(strcat(filepath,'/','full-movie/full-movie-ae'));
        folder = 'full-movie/full-movie-ae';
        file_type = 'full-movie-ae';
        
else
    
    if (invert_movie == 1)

        folder = 'full-movie/full-movie-invertido';
        file_type = strcat('bin-size-',int2str(bin_size),'-','full-movie-invertido');
        category{2} = strcat(category{2},'-invertido');

    else

        folder = 'full-movie/full-movie-original';
        file_type = 'full-movie-original';

    end
    
end

filepath = strcat(filepath,'/',folder,'/','v',int2str(video_index),'-','SP','-',date,'-','sitio',int2str(site_index),'-',channel,'-',int2str(start_time),'-',int2str(end_time),'-',file_type);

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
M = 3;

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

for i=1:3
   
    trials_label(i).labels = find(Spass.stimIds == i); 
    
end

Spass.spike_times = Spass.spike_times./ 32000;

categories(1:3) = struct('label','','P',P,'trials',zeros(P));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   BEGIN CATEGORIES LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Begin categories loop...');

for i=1:3
    
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

if ae == 1
    
    categories(2) = [] ;
    M = 2;
    
end

categories = categories';

sites = struct('label','unit_001','recording_tag','episodic','time_scale',time_scale,'time_resolution',time_resolution,'si_unit','none','si_prefix',si_prefix);

X = struct('M',M,'N',N,'sites',sites,'categories',categories);

%X = struct('M',M,'N',N,'sites',struct('label','unit_001','recording_tag','episodic','time_scale',time_scale,'time_resolution',time_resolution,'si_unit','none','si_prefix',si_prefix),'categories',struct('label','','P',P,'trials',struct('start_time',start_time,'end_time',end_time,'Q',Q,'list',list)));

stawrite(X,filepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END');

toc

end