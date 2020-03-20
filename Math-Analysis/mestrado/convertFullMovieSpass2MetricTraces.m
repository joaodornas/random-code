function X = convertFullMovieSpass2MetricTraces(date,site_index,channel,registro,video_index,start_time,end_time,invert_movie,bin_size,nRepetitions,intervals,stimulus_start_time,remove_latency)

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
        
        switch invert_movie
            
            case 1
                
                category{2} = 'Video 1: Swans. REVERSE';
                
            case 0
                
                category{2} = 'Video 1: Swans. INVERTIDO';
                
        end
                
        %category{3} = 'Video 1: Spontaneous Activity';
                
    case 2
        
        category{1} = 'Video 2: Mountains. DIRECT';
        
        switch invert_movie
            
            case 1
                
                category{2} = 'Video 2: Mountains. REVERSE';
                
            case 0
                
                category{2} = 'Video 2: Mountains. INVERTIDO';
                
        end
                
        %category{3} = 'Video 2: Spontaneous Activity';
        
    case 3
        
        category{1} = 'Video 3: Lizard. DIRECT';
        
        switch invert_movie
            
            case 1
                
                category{2} = 'Video 3: Lizard. REVERSE';
                
            case 0
                
                category{2} = 'Video 3: Lizard. INVERTIDO';
    
        end

        %category{3} = 'Video 3: Spontaneous Activity';
        
    case 4
        
        category{1} = 'Video 4: Beetle. DIRECT';
        
        switch invert_movie
            
            case 1
                
                category{2} = 'Video 4: Beetle. REVERSE';
                
            case 0
                
                category{2} = 'Video 4: Beetle. INVERTIDO';
    
        end
        
        %category{3} = 'Video 4: Spontaneous Activity';
        
end

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
M = 2 ;

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

for i=1:2
   
    trials_label(i).labels = find(Spass.stimIds == i); 
    
end

Spass.spike_times = Spass.spike_times./ 32000;

categories(1:2) = struct('label','','P',P,'trials',zeros(P));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% round( (end_time/1000 - start_time/1000) / (interval/1000) )
%%%   BEGIN CATEGORIES LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Begin categories loop...');

cell_label = strcat('v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel);
kernelString = strcat(cell_label,'-psth-kernel.mat');
kernel = load(kernelString);
latencyPoints = kernel.histData(1).latencyPoints
latencia = latencyPoints - (stimulus_start_time/1000)
    
    for i=1:2

        disp(strcat('Category .',int2str(i)));

        disp('...get trials');

        trials_spikes = Spass.spike_times(trials_label(i).labels(:),:);

        nTrials(i) = size(trials_spikes,1);

        for k=1:nTrials(i)

            spikes = trials_spikes(k,:);
            spikes = sort(spikes);
            spikes = spikes(spikes>0);
            
            if i == 2 && remove_latency == 1, spikes = spikes - latencia; end
            
            spikes = spikes(spikes>0);

            spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));
                 
            if i == 2 && invert_movie == 1

                nBins = (end_time - start_time)/bin_size;

                disp('...invert spike train ');
                
                spikes_inverted = invertSpikeTrain(spikes,start_time,nBins,bin_size);
                spikes = [];
                spikes = spikes_inverted.spikes; 

            end        
            
            s = round( (end_time/1000 - start_time/1000) / (intervals/1000) ); 
            
            cond(i).trials(k).tracesSpikes(N,1:s) = struct('start_time',start_time/1000,'end_time',end_time/1000,'Q',0,'list',zeros(1));
            
            s = 0;
            for t=start_time/1000:intervals/1000:end_time/1000

                s = s + 1;

                cond(i).trials(k).tracesSpikes(N,s).list = spikes(spikes>=t & spikes<=(t+intervals/1000));
                cond(i).trials(k).tracesSpikes(N,s).Q = size(spikes(spikes>=t & spikes<=(t+intervals/1000)),2);
                cond(i).trials(k).tracesSpikes(N,s).start_time = t;
                cond(i).trials(k).tracesSpikes(N,s).end_time = (t+intervals/1000);

            end

        end
        
           
    end
    
    for i=1:2
        
        for k=1:s

            cond(i).traces(k).trialsSpikes(N,1:nTrials(i)) = struct('start_time',start_time/1000,'end_time',end_time/1000,'Q',0,'list',zeros(1));
            
            for j=1:nTrials(i)
                               
                cond(i).traces(k).trialsSpikes(N,j).list = cond(i).trials(j).tracesSpikes(N,k).list;

                cond(i).traces(k).trialsSpikes(N,j).Q = size(cond(i).trials(j).tracesSpikes(N,k),2);

                cond(i).traces(k).trialsSpikes(N,j).start_time = cond(i).trials(j).tracesSpikes(N,k).start_time;

                cond(i).traces(k).trialsSpikes(N,j).end_time = cond(i).trials(j).tracesSpikes(N,k).end_time;
                
                cond(i).traces(k).start_time = cond(i).trials(j).tracesSpikes(N,k).start_time;
                
                cond(i).traces(k).end_time = cond(i).trials(j).tracesSpikes(N,k).end_time;

            end

        end
        
    end
    
    
    for k=1:s
         
        for i=1:2
        
            clear trials;
       
            for j=1:nTrials(i)
                      
                trials(N,j).list = cond(i).traces(k).trialsSpikes(N,j).list;
                trials(N,j).Q = size(trials(N,j).list,2);
                trials(N,j).start_time = cond(i).traces(k).trialsSpikes(N,j).start_time;
                trials(N,j).end_time = cond(i).traces(k).trialsSpikes(N,j).end_time;
            
            end
            
            trials = trials';
            
            categories(i).label = category{i};
            categories(i).P = nTrials(i);
            categories(i).trials = trials;            
        
            start_time = cond(i).traces(k).start_time;
        
            end_time = cond(i).traces(k).end_time;
            
        end

            filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));

            if (invert_movie == 1)

                folder = 'full-movie/trecho-invertido';
                file_type = strcat('bin-size-',int2str(bin_size),'-','trecho-invertido');

            else

                folder = 'full-movie/trecho-original';
                file_type = 'trecho-original';

            end

            filepath = strcat(filepath,'/',folder,'/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel,'-',int2str(start_time*1000),'-',int2str(end_time*1000),'-',int2str(intervals),'-',file_type);

            categories = categories';

            sites = struct('label','unit_001','recording_tag','episodic','time_scale',time_scale,'time_resolution',time_resolution,'si_unit','none','si_prefix',si_prefix);

            X = struct('M',M,'N',N,'sites',sites,'categories',categories);

            stawrite(X,filepath);

    end
         


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END');

toc

end