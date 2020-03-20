function convertMemoriaSpass2Metric(date,site_index,channel,registro,video_index,start_time,end_time,nRepetitions)

tic

disp('BEGIN');

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

Spass.spike_times = Spass.spike_times./ 32000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   SET CATEGORIES LABELS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Set categories labels...');

categories(1:3) = struct('label','','P',P,'trials',zeros(P));

for i=1:3

    switch i

        case 1

            category{1} = 'Video Direct (anterior: direct)';
            category{2} = 'Video Direct (anterior: reverse)';
            category{3} = 'Video Direct (anterior: spontaneous)';
            category{4} = '-direct';

        case 2

            category{1} = 'Video Reverse (anterior: direct)';
            category{2} = 'Video Reverse (anterior: reverse)';
            category{3} = 'Video Reverse (anterior: spontaneous)';
            category{4} = '-reverse';

        case 3

            category{1} = 'Video Spontaneous (anterior: direct)';
            category{2} = 'Video Spontaneous (anterior: reverse)';
            category{3} = 'Video Spontaneous (anterior: spontaneous)';        
            category{4} = '-spontaneous';

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));

    filepath = strcat(filepath,'/','memoria-movie','/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel);
    
    filepath = strcat(filepath,'-memoria-movie',category{4});

%%%   SET X STRUCTURE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(strcat('Set X structure (category:',int2str(i),')...'));

    for j=1:3
    
        trials_spikes = Spass.spike_times(video_pairs(i,j).trials(:),:);    
        
        nTrials = size(trials_spikes,1);
   
        trials(N,1:nTrials) = struct('start_time',start_time/1000,'end_time',end_time/1000,'Q',0,'list',zeros(size(trials_spikes,2)));
    
        for k=1:nTrials
    
            spikes = trials_spikes(k,:);
            spikes = spikes(spikes>0);

            spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));
            spikes = sort(spikes);

            trials(N,k).Q = size(spikes,2);
            trials(N,k).list = spikes;
                   
        end
 
        trials = trials';
        
        categories(j).label = category{j};
        categories(j).P = nTrials;
        categories(j).trials = trials;
    
        clear trials;

    end

    categories = categories';

    sites = struct('label','unit_001','recording_tag','episodic','time_scale',time_scale,'time_resolution',time_resolution,'si_unit','none','si_prefix',si_prefix);

    X = struct('M',M,'N',N,'sites',sites,'categories',categories);

    %X = struct('M',M,'N',N,'sites',struct('label','unit_001','recording_tag','episodic','time_scale',time_scale,'time_resolution',time_resolution,'si_unit','none','si_prefix',si_prefix),'categories',struct('label','','P',P,'trials',struct('start_time',start_time,'end_time',end_time,'Q',Q,'list',list)));

    stawrite(X,filepath);

    clear X category filepath;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

disp('END');

toc

end