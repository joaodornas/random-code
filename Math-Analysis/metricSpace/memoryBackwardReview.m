function memoryBackwardReview

tic

disp('Load variables...');

date = { '12-08-20' '12-08-27' '12-08-29' '12-08-29' '12-08-29' '12-08-30' '12-08-30' '12-08-30' '12-08-30' '12-08-30' '12-08-31' '12-08-31' '12-08-31' '12-08-31' '12-08-31' '12-08-31' '12-09-03' '12-09-03' '12-09-04' '12-09-04' '12-09-05' '12-09-05' '12-09-06' '12-09-06' '12-10-16' '12-12-17' '12-12-17' '12-12-17'};
site_index = [ 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1];
channel = { 'E1' 'E2' 'E1' 'E1' 'E1' 'E1' 'E3' 'E1' 'E3' 'E3' 'E1' 'E1' 'E1' 'E3' 'E3' 'E3' 'E3' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E3'};
registro = { 'nsp003a01_1b' 'nsp005a01_2b' 'nsp006a01_1b' 'nsp006a02_1b' 'nsp006b01_1b' 'nsp007a01_1b' 'nsp007a01_2b' 'nsp007a02_1b' 'nsp007a02_2b' 'nsp007a03_2a' 'nsp008a01_1b' 'nsp008a02_1a' 'nsp008a03_1b' 'nsp008a01_2b' 'nsp008a02_2b' 'nsp008a03_2a' 'nsp009a01_2b' 'nsp009b01_1a' 'nsp010a1_1b' 'nsp010a02_1b' 'nsp011a01_1a' 'nsp011a02_1b' 'nsp012a01_1b' 'nsp012a02_1a' 'nps013a01_1b' 'nsp033a09_1b' 'nsp033a09_1c' 'nsp033a09_3b' };
video_index = [ 3 3 3 4 1 4 4 1 1 3 1 3 4 1 3 4 2 4 2 4 3 4 3 4 1 311 321 331];

%start_time = 0;
start_time = 500;
end_time = [3500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500];
%end_time = [4000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000];
nRepetitions = 60; 
bin_size = 1;

Qmax = [0 0.0625 0.125 0.25 0.5 1 2 4 8 16 32 64 128 256 512];
Rmax = 1./Qmax;

caminho = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/memoryReview/500-9500/');

disp('Convert SPASS to Metric...');

for r=1:length(registro)

    datasetO{r} = filename(char(date{r}),site_index(r),char(channel{r}),char(registro{r}),video_index(r),start_time,end_time(r),0,bin_size,nRepetitions,0);
    datasetI{r} = filename(char(date{r}),site_index(r),char(channel{r}),char(registro{r}),video_index(r),start_time,end_time(r),1,bin_size,nRepetitions,0);

    datasetOae{r} = filename(char(date{r}),site_index(r),char(channel{r}),char(registro{r}),video_index(r),start_time,end_time(r),0,bin_size,nRepetitions,1);
    datasetIae{r} = filename(char(date{r}),site_index(r),char(channel{r}),char(registro{r}),video_index(r),start_time,end_time(r),1,bin_size,nRepetitions,1);

end

% for r=1:length(registro)
% 
%     convertFullMovieSpass2Metric(char(date{r}),site_index(r),char(channel{r}),char(registro{r}),video_index(r),start_time,end_time(r),0,bin_size,nRepetitions,0);
%     convertFullMovieSpass2Metric(char(date{r}),site_index(r),char(channel{r}),char(registro{r}),video_index(r),start_time,end_time(r),1,bin_size,nRepetitions,0);
% 
%     convertFullMovieSpass2Metric(char(date{r}),site_index(r),char(channel{r}),char(registro{r}),video_index(r),start_time,end_time(r),0,bin_size,nRepetitions,1);
%     convertFullMovieSpass2Metric(char(date{r}),site_index(r),char(channel{r}),char(registro{r}),video_index(r),start_time,end_time(r),1,bin_size,nRepetitions,1);
% 
% end

% disp('Calculate Metric Space...');
% 
% for r=1:length(registro)
% 
%     metricMovie(strcat(char(datasetO{r}),'.stam'), start_time, end_time(r), nRepetitions);
% 
%     close all;
% 
%     metricMovie(strcat(char(datasetI{r}),'.stam'), start_time, end_time(r), nRepetitions);
% 
%     close all;
% 
%     metricMovie(strcat(char(datasetOae{r}),'.stam'), start_time, end_time(r), nRepetitions);
% 
%     close all;
% 
%     metricMovie(strcat(char(datasetIae{r}),'.stam'), start_time, end_time(r), nRepetitions);
% 
%     close all;
% 
% end

disp('Consolidate Results...');

diminuiu = 0;
aumentou = 0;
manteve = 0;
diminuiuAE = 0;
aumentouAE = 0;
manteveAE = 0;
diminuiuAll = 0;
aumentouAll = 0;
manteveAll = 0;
diminuiuAllae = 0;
aumentouAllae = 0;
manteveAllae = 0;
for r=1:length(registro)
   
    metricO = load(strcat(char(datasetO{r}),'.mat'));
    metricI = load(strcat(char(datasetI{r}),'.mat'));
    
    entropyO(r) = metricO.metric_analysis.max_info.max_info_tpmc_T;
    entropyI(r) = metricI.metric_analysis.max_info.max_info_tpmc_T;

    hcountO(r) = metricO.metric_analysis.max_info.Hcount_tpmc_info_T;
    hcountI(r) = metricI.metric_analysis.max_info.Hcount_tpmc_info_T;

    resolutionO(r) = Rmax(metricO.metric_analysis.max_info.max_info_tpmc_idx_T);
    resolutionI(r) = Rmax(metricI.metric_analysis.max_info.max_info_tpmc_idx_T);
    
    if metricO.metric_analysis.max_info.max_info_tpmc_T > metricI.metric_analysis.max_info.max_info_tpmc_T 
        
        diminuiu = diminuiu + 1;
        
    elseif metricO.metric_analysis.max_info.max_info_tpmc_T < metricI.metric_analysis.max_info.max_info_tpmc_T 
        
        aumentou = aumentou + 1;

    else

        manteve = manteve + 1;
        
    end
    
    for e=1:length(Rmax)
    
        if metricO.metric_analysis.out_T(1,e).table.information(2,1).value > metricI.metric_analysis.out_T(1,e).table.information(2,1).value
            
            diminuiuAll = diminuiuAll + 1;
            
        elseif metricO.metric_analysis.out_T(1,e).table.information(2,1).value < metricI.metric_analysis.out_T(1,e).table.information(2,1).value
            
            aumentouAll = aumentouAll + 1;

        else
        
            manteveAll = manteveAll + 1;
            
        end
        
    end
        
    metricOae = load(strcat(char(datasetOae{r}),'.mat'));
    metricIae = load(strcat(char(datasetIae{r}),'.mat'));

    entropyOae(r) = metricOae.metric_analysis.max_info.max_info_tpmc_T;
    entropyIae(r) = metricIae.metric_analysis.max_info.max_info_tpmc_T;

    hcountOae(r) = metricOae.metric_analysis.max_info.Hcount_tpmc_info_T;
    hcountIae(r) = metricIae.metric_analysis.max_info.Hcount_tpmc_info_T;
    
    resolutionOae(r) = Rmax(metricOae.metric_analysis.max_info.max_info_tpmc_idx_T);
    resolutionIae(r) = Rmax(metricIae.metric_analysis.max_info.max_info_tpmc_idx_T);

    if metricOae.metric_analysis.max_info.max_info_tpmc_T > metricIae.metric_analysis.max_info.max_info_tpmc_T 
        
        diminuiuAE = diminuiuAE + 1;
        
    elseif metricOae.metric_analysis.max_info.max_info_tpmc_T < metricIae.metric_analysis.max_info.max_info_tpmc_T
        
        aumentouAE = aumentouAE + 1;
        
    else

        manteveAE = manteveAE + 1;

    end
    
    for e=1:length(Rmax)
    
        if metricOae.metric_analysis.out_T(1,e).table.information(2,1).value > metricIae.metric_analysis.out_T(1,e).table.information(2,1).value
            
            diminuiuAllae = diminuiuAllae + 1;
            
        elseif metricOae.metric_analysis.out_T(1,e).table.information(2,1).value < metricIae.metric_analysis.out_T(1,e).table.information(2,1).value
            
            aumentouAllae = aumentouAllae + 1;

        else

            manteveAllae = manteveAllae + 1;
            
        end
        
    end

end

results.entropyO = entropyO;
results.entropyI = entropyI;

results.lillieEntropyO = lillietest(entropyO);
results.lillieEntropyI = lillietest(entropyI);

[h,PV] = ttest(entropyO,entropyI);
results.EntropiesTTestH = h;
results.EntropiesTTestP = PV;

[PV,h] = ranksum(entropyO,entropyI);
results.EntropiesWKxH = h;
results.EntropiesWKxP = PV;

results.resolutionO = resolutionO;
results.resolutionI = resolutionI;

results.lillieResolutonO = lillietest(resolutionO);
results.lillieResolutionI = lillietest(resolutionI);

[h,PV] = ttest(resolutionO,resolutionI);
results.ResolutionTTestH = h;
results.ResolutionTTestP = PV;

[PV,h] = ranksum(resolutionO,resolutionI);
results.ResolutionWKxH = h;
results.ResolutionWKxP = PV;

results.diminuiu = diminuiu;
results.aumentou = aumentou;
results.manteve = manteve;

results.diminuiuAll = diminuiuAll;
results.aumentouAll = aumentouAll;
results.manteveAll = manteveAll;

results.entropyOae = entropyOae;
results.entropyIae = entropyIae;

results.lillieEntropyOae = lillietest(entropyOae);
results.lillieEntropyIae = lillietest(entropyIae);

[h,PV] = ttest(entropyOae,entropyIae);
results.EntropiesTTestHae = h;
results.EntropiesTTestPae = PV;

[PV,h] = ranksum(entropyOae,entropyIae);
results.EntropiesWKxHae = h;
results.EntropiesWKxPae = PV;

results.resolutionOae = resolutionOae;
results.resolutionIae = resolutionIae;

results.lillieResolutonOae = lillietest(resolutionOae);
results.lillieResolutionIae = lillietest(resolutionIae);

[h,PV] = ttest(resolutionOae,resolutionIae);
results.ResolutionTTestHae = h;
results.ResolutionTTestPae = PV;

[PV,h] = ranksum(resolutionOae,resolutionIae);
results.ResolutionWKxHae = h;
results.ResolutionWKxPae = PV;

results.diminuiuAE = diminuiuAE;
results.aumentouAE = aumentouAE;
results.manteveAE = manteveAE;

results.diminuiuAllae = diminuiuAllae;
results.aumentouAllae = aumentouAllae;
results.manteveAllae = manteveAllae;

maxVal = max(max(entropyO),max(hcountO));

f = figure;
plot(hcountO,entropyO,'ro');
xlim([0 maxVal]);
ylim([0 maxVal]);
xlabel('H-count-O');
ylabel('H-timing-O');
hold on;
plot([0 maxVal],[0 maxVal],'k');

print(f,'-depsc',strcat(caminho,'/hcountOXhtimingO-CC'));

maxVal = max(max(entropyI),max(hcountI));

g = figure;
plot(hcountI,entropyI,'ro');
xlim([0 maxVal]);
ylim([0 maxVal]);
xlabel('H-count-I');
ylabel('H-timing-I');
hold on;
plot([0 maxVal],[0 maxVal],'k');

print(g,'-depsc',strcat(caminho,'/hcountIXhtimingI-CC'));

maxVal = max(max(entropyOae),max(hcountOae));

h = figure;
plot(hcountOae,entropyOae,'ro');
xlim([0 maxVal]);
ylim([0 maxVal]);
xlabel('H-count-Oae');
ylabel('H-timing-Oae');
hold on;
plot([0 maxVal],[0 maxVal],'k');

print(h,'-depsc',strcat(caminho,'/hcountOaeXhtimingOae-CC'));

maxVal = max(max(entropyIae),max(hcountIae));

l = figure;
plot(hcountIae,entropyIae,'ro');
xlim([0 maxVal]);
ylim([0 maxVal]);
xlabel('H-count-Iae');
ylabel('H-timing-Iae');
hold on;
plot([0 maxVal],[0 maxVal],'k');

print(l,'-depsc',strcat(caminho,'/hcountIaeXhtimingIae-CC'));


save(strcat(caminho,'results-metric-movies'),'results');

disp('FINISH');

toc

function dataset = filename(date,site_index,channel,registro,video_index,start_time,end_time,invert_movie,bin_size,nRepetitions,ae)

path = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/memoryReview/500-9500/');

if ae == 1
        
    if (invert_movie == 1)

        mkdir(strcat(path,'/','full-movie-invertido-ae'));
        folder = 'full-movie-invertido-ae';
        file_type = strcat('bin-size-',int2str(bin_size),'-','full-movie-invertido-ae');

    else

        mkdir(strcat(path,'/','full-movie-original-ae'));
        folder = 'full-movie-original-ae';
        file_type = 'full-movie-original-ae';

    end
        
else
    
    if (invert_movie == 1)

        mkdir(strcat(path,'/','full-movie-invertido'));
        folder = 'full-movie-invertido';
        file_type = strcat('bin-size-',int2str(bin_size),'-','full-movie-invertido');

    else

        mkdir(strcat(path,'/','full-movie-original'));
        folder = 'full-movie-original';
        file_type = 'full-movie-original';

    end
    
end

dataset = strcat(path,'/',folder,'/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel,'-',int2str(start_time),'-',int2str(end_time),'-',file_type);

end

function convertFullMovieSpass2Metric(date,site_index,channel,registro,video_index,start_time,end_time,invert_movie,bin_size,nRepetitions,ae)

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

filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/memoryReview/500-9500/');

if ae == 1
        
    if (invert_movie == 1)

        path = filename(date,site_index,channel,registro,video_index,start_time,end_time,1,bin_size,nRepetitions,1);
        category{2} = strcat(category{2},'-invertido');
        
        M = 3;

    else

        path = filename(date,site_index,channel,registro,video_index,start_time,end_time,0,bin_size,nRepetitions,1);
        
        M = 3;

    end
        
else
    
    if (invert_movie == 1)

        path = filename(date,site_index,channel,registro,video_index,start_time,end_time,1,bin_size,nRepetitions,0);
        category{2} = strcat(category{2},'-invertido');
        
        M = 2;

    else

        path = filename(date,site_index,channel,registro,video_index,start_time,end_time,0,bin_size,nRepetitions,0);
        
        M = 2;

    end
    
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
%M = 3;

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

stawrite(X,path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END');

end


function metricMovie(dataset, start_time, end_time, nRepetitions)

disp('READING DATA');
% READ DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%
X = staread(strrep(dataset,'/',filesep));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('SETTING OPTIONS');
% OPTIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.entropy_estimation_method = {'plugin','tpmc','jack'};
%opts.variance_estimation_method = {'jack'};

%opts.unoccupied_bins_strategy = -1; % Ignore unoccupied bins
opts.unoccupied_bins_strategy = 0; % Use an unoccupied bin only if its row and column are occupied
%opts.unoccupied_bins_strategy = 1; % Use all bins

opts.parallel = 1;
opts.possible_words = 'unique';

opts.start_time = start_time / 1000;
opts.end_time = end_time / 1000;
opts.shift_cost = [0 2.^(-4:9)];
%opts.label_cost = [0 1 2];
opts.clustering_exponent = -2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = nRepetitions;

disp('METRIC TIMING');
% METRIC TIMING + SHUFFLE    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.metric_family = 0;

disp('...computing metric space timing');
[out_T,opts_used] = metric(X,opts);

for q_idx=1:length(opts.shift_cost)

    info_plugin_T(q_idx) = out_T(q_idx).table.information(1).value;
    info_tpmc_T(q_idx) = out_T(q_idx).table.information(2).value;
    info_jack_T(q_idx) = out_T(q_idx).table.information(3).value;

end

[max_info_plugin_T,max_info_plugin_idx_T] = max(info_plugin_T);
[max_info_tpmc_T,max_info_tpmc_idx_T] = max(info_tpmc_T);
[max_info_jack_T,max_info_jack_idx_T] = max(info_jack_T);

Hcount_plugin_T = info_plugin_T(1);
Hcount_tpmc_T = info_tpmc_T(1);
Hcount_jack_T = info_jack_T(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
disp('METRIC INTERVAL');
% METRIC INTERVAL + SHUFFLE    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.metric_family = 1;

disp('...computing metric space interval');
[out_I,opts_used] = metric(X,opts);

for q_idx=1:length(opts.shift_cost)

    info_plugin_I(q_idx) = out_I(q_idx).table.information(1).value;
    info_tpmc_I(q_idx) = out_I(q_idx).table.information(2).value;
    info_jack_I(q_idx) = out_I(q_idx).table.information(3).value;

end

[max_info_plugin_I,max_info_plugin_idx_I] = max(info_plugin_I);
[max_info_tpmc_I,max_info_tpmc_idx_I] = max(info_tpmc_I);
[max_info_jack_I,max_info_jack_idx_I] = max(info_jack_I);

Hcount_plugin_I = info_plugin_I(1);
Hcount_tpmc_I = info_tpmc_I(1);
Hcount_jack_I = info_jack_I(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('CALCULATE H-BIAS VIA SHUFFLE');
% CALCULATE H-BIAS VIA SHUFFLE %%%%%%%%%%%%%%%

disp('...computing metric space timing shuffle');

opts.entropy_estimation_method = {'plugin'};

rand('state',0);
opts.metric_family = 0;
[Y_T,SHUF_T,opts_used] = metric_shuf(X,opts,S);

SHUF_T = SHUF_T';

disp('...computing metric space interval shuffle');

rand('state',0);
opts.metric_family = 1;
[Y_I,SHUF_I,opts_used] = metric_shuf(X,opts,S);

SHUF_I = SHUF_I';

for q_idx=1:length(opts.shift_cost)
    
    for i=1:S
        
        HShuffle_T(i,q_idx) = SHUF_T(i,q_idx).table.information.value;
        
    end

    for i=1:S
        
        HShuffle_I(i,q_idx) = SHUF_I(i,q_idx).table.information.value;
        
    end

end

HBias_T = mean(HShuffle_T,1);
HBias_std_T = std(HShuffle_T,[],1);
HBias_I = mean(HShuffle_I,1);
HBias_std_I = std(HShuffle_I,[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% leave-one-out Jackknife 
disp('leave-one-out Jackknife');

opts.entropy_estimation_method = {'plugin'};
opts.metric_family = 0;
[out_unjk,jk,opts_used] = metric_jack(X,opts);

P_total = size(jk,1);

temp_info_jk = zeros(P_total,length(opts.shift_cost));

for q_idx=1:length(opts.shift_cost)
    
  info_unjk(q_idx)= out_unjk(q_idx).table.information.value;
  
  for p=1:P_total
      
    temp_info_jk(p,q_idx) = jk(p,q_idx).table.information.value;
    
  end
  
end

info_jk_sem = sqrt((P_total-1)*var(temp_info_jk,1,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('SAVING DATA');
% SAVE DATA    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_info = struct('max_info_plugin_idx_T',max_info_plugin_idx_T,'max_info_plugin_T',max_info_plugin_T,'max_info_plugin_idx_I',max_info_plugin_idx_I,'max_info_plugin_I',max_info_plugin_I,'Hcount_plugin_info_T',Hcount_plugin_T,'Hcount_plugin_info_I',Hcount_plugin_I,'max_info_tpmc_idx_T',max_info_tpmc_idx_T,'max_info_tpmc_T',max_info_tpmc_T,'max_info_tpmc_idx_I',max_info_tpmc_idx_I,'max_info_tpmc_I',max_info_tpmc_I,'Hcount_tpmc_info_T',Hcount_tpmc_T,'Hcount_tpmc_info_I',Hcount_tpmc_I,'max_info_jack_idx_T',max_info_jack_idx_T,'max_info_jack_T',max_info_jack_T,'max_info_jack_idx_I',max_info_jack_idx_I,'max_info_jack_I',max_info_jack_I,'Hcount_jack_info_T',Hcount_jack_T,'Hcount_jack_info_I',Hcount_jack_I,'HBias_T',HBias_T,'HBias_std_T',HBias_std_T,'HBias_I',HBias_I,'HBias_std_I',HBias_std_I);

metric_analysis = struct('X', X, 'out_T', out_T, 'out_I', out_I, 'max_info', max_info, 'opts', opts, 'Y_I', Y_I, 'Y_T', Y_T, 'SHUF_I', SHUF_I, 'SHUF_T', SHUF_T);

dataset = dataset(1:end-5);
save(dataset,'metric_analysis');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING RASTER');
% RASTER PLOT    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure;
set(gcf,'name',['Metric ' dataset]); 

staraster(X,[opts.start_time opts.end_time]);
title('Raster plot');

print(f1,'-depsc',strcat(dataset,'-raster-plot'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING DISTANCE-MATRIX TIMING');
% DISTANCE-MATRIX TIMING PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2 = figure;
imagesc(out_T(max_info_plugin_idx_T).d);
xlabel('Spike train index');
ylabel('Spike train index');
title('Distance matrix at maximum information - timing');
print(f2,'-depsc',strcat(dataset,'-distance-matrix-timing'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING CONFUSION-MATRIX TIMING');
% CONFUSION-MATRIX TIMING PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3 = figure;
imagesc(out_T(max_info_plugin_idx_T).cm);
xlabel('Category index');
ylabel('Category index');
title('Confusion matrix from clustering of distances - timing');
print(f3,'-depsc',strcat(dataset,'-confusion-matrix-timing'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING DISTANCE-MATRIX INTERVAL');
% DISTANCE-MATRIX INTERVAL PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%
f4 = figure;
imagesc(out_I(max_info_plugin_idx_I).d);
xlabel('Spike train index');
ylabel('Spike train index');
title('Distance matrix at maximum information - interval');
print(f4,'-depsc',strcat(dataset,'-distance-matrix-interval'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING CONFUSION-MATRIX INTERVAL');
% CONFUSION-MATRIX INTERVAL PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f5 = figure;
imagesc(out_I(max_info_plugin_idx_I).cm);
xlabel('Category index');
ylabel('Category index');
title('Confusion matrix from clustering of distances - interval');
print(f5,'-depsc',strcat(dataset,'-confusion-matrix-interval'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING INFO TIMING-INTERVAL');
% INFO TIMING-INTERVAL PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f6 = figure;
plot(1:length(opts.shift_cost),info_plugin_T,'b');
hold on;
plot(1:length(opts.shift_cost),info_tpmc_T,'b--');
plot(1:length(opts.shift_cost),info_jack_T,'b.-');
plot(1:length(opts.shift_cost),info_plugin_I,'r');
plot(1:length(opts.shift_cost),info_tpmc_I,'r--');
plot(1:length(opts.shift_cost),info_jack_I,'r.-');
hold off;

set(gca,'xtick',1:length(opts.shift_cost));
set(gca,'xticklabel',opts.shift_cost);
set(gca,'xlim',[1 length(opts.shift_cost)]);
set(gca,'ylim',[-0.5 2.5]);

xlabel('Temporal precision (1/sec)');
ylabel('Information (bits)');

legend('No correction - Timing','Trevez-Panzeri - Timing','Jackknife - Timing','No correction - Interval','Trevez-Panzeri - Interval','Jackknife - Interval');
%legend('No correction - Timing','Jackknife - Timing','No correction - Interval','Jackknife - Interval'); 
print(f6,'-depsc',strcat(dataset,'-plot-info-timing-interval'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING INFO TIMING-JACK-SHUFFLE');
% INFO + SHUFFLE INFO PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f7 = figure;
errorbar(1:length(opts.shift_cost),info_unjk,2*info_jk_sem);
hold on;
errorbar(1:length(opts.shift_cost),HBias_T,2*HBias_std_T,'r');
hold off;
set(gca,'xtick',1:length(opts.shift_cost));
set(gca,'xticklabel',opts.shift_cost);
set(gca,'xlim',[1 length(opts.shift_cost)]);
set(gca,'ylim',[-0.5 2.5]);
xlabel('Temporal precision (1/sec)');
ylabel('Information (bits)');
legend('Original data (\pm 2 SE via Jackknife)','Shuffled data (\pm 2 SD)',...
       'location','best');

scalefig(gcf,1.5);
print(f7,'-depsc',strcat(dataset,'-plot-info-shuffle'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


end