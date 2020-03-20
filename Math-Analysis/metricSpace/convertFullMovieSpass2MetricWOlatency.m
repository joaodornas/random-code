function X = convertFullMovieSpass2MetricWOlatency(date,site_index,channel,registro,video_index,start_time,end_time,inverter_movie,bin_size,nRepetitions,ae)

tic

disp('BEGIN');

%%%   DEFINE CATEGORY LABEL   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('Define category label...');

if ae == 0
 
    nConditions = 2;
 
elseif ae == 1

    nConditions = 3;

end

if video_index > 100
    
    vn = round(video_index/100);
    
else
    
    vn = video_index;
    
end

switch vn
    
    case 1
        
        category{1} = 'Video 1: Swans. DIRECT WO Latency';
        category{2} = 'Video 1: Swans. REVERSE WO Latency';
        category{3} = 'Video 1: Spontaneous Activity WO Latency';
                
    case 2
        
        category{1} = 'Video 2: Mountains. DIRECT WO Latency';
        category{2} = 'Video 2: Mountains. REVERSE WO Latency';
        category{3} = 'Video 2: Spontaneous Activity WO Latency';
        
    case 3
        
        category{1} = 'Video 3: Lizard. DIRECT WO Latency';
        category{2} = 'Video 3: Lizard. REVERSE WO Latency';
        category{3} = 'Video 3: Spontaneous Activity WO Latency';
        
    case 4
        
        category{1} = 'Video 4: Beetle. DIRECT WO Latency';
        category{2} = 'Video 4: Beetle. REVERSE WO Latency';
        category{3} = 'Video 4: Spontaneous Activity WO Latency';
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   DEFINE FILE NAME AND PATH   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('Define file name and path...');

filepath = strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/metricMovieWOlatency/',int2str(start_time),'-',int2str(end_time),'/',registro,'-','v',int2str(video_index));

if ae == 1
        
        mkdir(strcat(filepath,'/','full-movie/full-movie-ae'));
        folder = 'full-movie/full-movie-ae';
        file_type = 'full-movie-ae-WO-latency';
        
else
    
    if (inverter_movie == 1)

        mkdir(strcat(filepath,'/','full-movie/full-movie-invertido'));
        folder = 'full-movie/full-movie-invertido';
        file_type = strcat('bin-size-',int2str(bin_size),'-','full-movie-invertido-WO-latency');
        category{2} = strcat(category{2},'-invertido');

    else

        mkdir(strcat(filepath,'/','full-movie/full-movie-original'));
        folder = 'full-movie/full-movie-original';
        file_type = 'full-movie-original-WO-latency';

    end
    
end

filepath = strcat(filepath,'/',folder,'/',registro,'-','v',int2str(video_index),'-',file_type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   SET METRIC SPACE PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('Set Metric Space parameters...');

%Escala de tempo em segundos
time_scale = 1;

%Escala da resolu??o temporal em segundos
time_resolution = 0.0001;

%Sistema Internacional de Medidas
si_prefix = 1;

%N?mero de classes de est?mulos
M = 2;

%N?mero de s?tios (canais)
N = 1;

%N?mero de trials (repeti??es por condi??o)
P = nRepetitions;

%Tempos inicial e final do trial
%start_time = 0;
%end_time = 0;

%N?mero de pontos no vetor list
%Q = 0;

latency_start_time = 500;

latency_end_time = latency_start_time + 300;

p = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   LOAD TRIALS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('Load trials...');
    
Spass = load(strcat('_',registro,'-','v',int2str(video_index),'.mat'));

for i=1:nConditions
   
    trials_label(i).labels = find(Spass.stimIds == i); 
    
end

Spass.spike_times = Spass.spike_times./ 32000;

categories(1:nConditions) = struct('label','','P',P,'trials',zeros(P));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   BEGIN CATEGORIES LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('Begin categories loop...');

w = 0;
for i=1:nConditions
  
    %disp(strcat('Category .',int2str(i)));

    %disp('...get trials');

    trials_spikes = Spass.spike_times(trials_label(i).labels(:),:);
        
    if ae == 1
        
        if i ~= 2
            
            w = w + 1;

            the_category = getCategory(trials_spikes,start_time,end_time,latency_start_time,latency_end_time,p,i,N,bin_size,category,inverter_movie);

            categories(w) = the_category;

        end
        
    elseif ae == 0
        
        the_category = getCategory(trials_spikes,start_time,end_time,latency_start_time,latency_end_time,p,i,N,bin_size,category,inverter_movie);
        
        categories(i) = the_category;
        
    end
            
        clear the_category;
  
end


%disp('End categories loop...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CREATE METRIC SPACE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('Create Metric Space data...');

categories = categories';

sites = struct('label','unit_001','recording_tag','episodic','time_scale',time_scale,'time_resolution',time_resolution,'si_unit','none','si_prefix',si_prefix);

X = struct('M',M,'N',N,'sites',sites,'categories',categories);

%X = struct('M',M,'N',N,'sites',struct('label','unit_001','recording_tag','episodic','time_scale',time_scale,'time_resolution',time_resolution,'si_unit','none','si_prefix',si_prefix),'categories',struct('label','','P',P,'trials',struct('start_time',start_time,'end_time',end_time,'Q',Q,'list',list)));

stawrite(X,filepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END');

function the_category = getCategory(trials_spikes,start_time,end_time,latency_start_time,latency_end_time,p,i,N,bin_size,category_name,inverter_movie)

    
            alltrials = trials_spikes;

            alltrials = reshape(alltrials.',[],1);

            alltrials = alltrials.';

            latency_time = latency(alltrials,latency_start_time,latency_end_time,p); 

            if i == 3

                latency_time = 0;

            end

            clear alltrials;

            disp(strcat('Latencia:',num2str(latency_time)));

            nTrials = size(trials_spikes,1);

            trials(N,1:nTrials) = struct('start_time',start_time/1000,'end_time',end_time/1000,'Q',0,'list',zeros(size(trials_spikes,2)));

            for k=1:nTrials

                spikes = trials_spikes(k,:);

                a_e = spikes(spikes>=(0/1000) & spikes<(latency_start_time/1000));

                after_ae = spikes(spikes>=(latency_start_time/1000) & spikes<(end_time/1000)) - latency_time;

                spikes = [];

                spikes = [a_e, after_ae];

                spikes = spikes(spikes>0);

                spikes = sort(spikes);

                spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));

                if (i == 2) & (inverter_movie == 1)

                    nBins = (end_time - start_time)/bin_size;

                    %disp(strcat('...invert trial .',int2str(k)));

                    spikes_inverted = invertSpikeTrain(spikes,start_time,nBins,bin_size);
                    spikes = spikes_inverted.spikes;

                end

                %disp(strcat('...set trial .',int2str(k)));

                trials(N,k).Q = size(spikes,2);
                trials(N,k).list = spikes;

            end

            trials = trials';

            the_category.label = category_name{i};
            the_category.P = nTrials;
            the_category.trials = trials;

            clear trials;
        

end

toc

end