function metricReviewCC

tic

disp('Load variables...');

date = { '12-11-01' '12-11-01' '12-11-06' '12-11-06' '12-11-07' '12-11-07' '12-11-08' '12-11-08' '12-11-09' '12-11-09' '12-11-13' '12-11-13' '12-12-03' '12-12-03' '12-12-04' '12-12-04' '12-12-05' '12-12-05' '12-12-07' '12-12-07' '12-12-10' '12-12-10' '12-12-12' '12-12-12' '12-12-13' '12-12-13' '12-12-17' '12-12-17' '12-12-17' '12-12-17' '12-12-17' '12-12-17' '12-12-17' '12-12-17' '12-12-17' '12-12-17' '12-12-18' '12-12-18' '12-12-19' '12-12-19'};
site_index = [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
channel = { 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E3' 'E3' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E2' 'E2' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E3' 'E1' 'E1' 'E3' 'E1' 'E1' 'E1' 'E3' 'E2' 'E2' 'E2' 'E2' };
registro = { 'nsp020a01_1b' 'nsp020a02_1b' 'nsp021a01_1a' 'nsp021a02_1a' 'nsp022a01_1b' 'nsp022a02_1b' 'nsp023a01_1a' 'nsp023a02_1a' 'nsp024a01_1a' 'nsp024a02_1a' 'nsp025a01_1b' 'nsp025a02_1a' 'nsp026a01_3b' 'nsp026a02_3a' 'nsp027a03_1a' 'nsp027a04_1b' 'nsp028a03_1a' 'nsp028a04_1b' 'nsp029a03_1a' 'nsp029a04_1b' 'nsp030a03_2a' 'nsp030a04_2b' 'nsp031b03_1b' 'nsp031b04_1b' 'nsp032a03_1a' 'nsp032a04_1a' 'nsp033a04_1b' 'nsp033a04_1c' 'nsp033a04_3a' 'nsp033a05_1b' 'nsp033a05_1c' 'nsp033a05_3a' 'nsp033a06_1b' 'nsp033a06_1c' 'nsp033a07_1b' 'nsp033a07_3b' 'nsp034a04_2b' 'nsp034a05_2b' 'nsp035a03_2b' 'nsp035a04_2b'};
video_index = [ 5 6 5 6 5 6 1 4 6 3 2 1 8 11 10 9 6 8 8 3 5 8 3 4 8 2 821 811 831 511 521 531 812 822 512 532 8 11 8 1];

start_time = 500;
end_time = [9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500 9500]; 
bin_size = 30;

Qmax = [0 0.0625 0.125 0.25 0.5 1 2 4 8 16 32 64 128 256 512];
Rmax = 1./Qmax;

caminho = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/metricReviewCC/500-9500');

disp('Convert SPASS to Metric...');

for r=1:length(registro)

    datasetO{r} = filename(char(date{r}),site_index(r),char(channel{r}),video_index(r),start_time,end_time(r));

end

% for r=1:length(registro)
% 
%     convertFullMovieSpass2Metric(char(date{r}),site_index(r),char(channel{r}),char(registro{r}),video_index(r),start_time,end_time(r));
% 
% end

% disp('Calculate Metric Space...');
% 
% for r=1:length(registro)
% 
%     if r < 13, nRepetitions = 60; end
% 
%     if r > 12, nRepetitions = 40; end
% 
%     metricMovie(strcat(char(datasetO{r}),'.stam'), start_time, end_time(r), nRepetitions);
% 
%     close all;
% 
% end

disp('Consolidate Results...');


for r=1:length(registro)
   
    metricO = load(strcat(char(datasetO{r}),'.mat'));
    
    entropyO(r) = metricO.metric_analysis.max_info.max_info_tpmc_T;
    
    hcountO(r) = metricO.metric_analysis.max_info.Hcount_tpmc_info_T;
    
end

maxVal = max(max(entropyO),max(hcountO));

f = figure;
plot(hcountO,entropyO,'ro');
xlim([0 maxVal]);
ylim([0 maxVal]);
xlabel('H-count');
ylabel('H-timing');
hold on;
plot([0 maxVal],[0 maxVal],'k');

print(f,'-depsc',strcat(caminho,'/hcountXhtiming-CC'));

results.entropyO = entropyO;
results.hcountO = hcountO;

save(strcat(caminho,'/results-metric-CC'),'results');

disp('FINISH');

toc

function dataset = filename(date,site_index,channel,video_index,start_time,end_time)

path = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/metricReviewCC/500-9500');

mkdir(strcat(path,'/','full-movie-original'));
folder = 'full-movie-original';
file_type = 'full-movie-original';

dataset = strcat(path,'/',folder,'/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel,'-',int2str(start_time),'-',int2str(end_time),'-',file_type);

end

function convertFullMovieSpass2Metric(date,site_index,channel,registro,video_index,start_time,end_time)

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
        
        category{1} = 'Video 1: Swans. CRF';
        category{2} = 'Video 1: Swans. nCRF';
        category{3} = 'Video 1: Swans. ExtraNCRF';
                
    case 2
        
        category{1} = 'Video 2: Mountains. CRF';
        category{2} = 'Video 2: Mountains. nCRF';
        category{3} = 'Video 2: Mountains. ExtraNCRF';
        
    case 3
        
        category{1} = 'Video 3: Lizard. CRF';
        category{2} = 'Video 3: Lizard. nCRF';
        category{3} = 'Video 3: Lizard. ExtraNCRF';
        
    case 4
        
        category{1} = 'Video 4: Beetle. CRF';
        category{2} = 'Video 4: Beetle. nCRF';
        category{3} = 'Video 4: Beetle. ExtraNCRF';

    case 5
        
        category{1} = 'Video 5: Chimp. CRF';
        category{2} = 'Video 5: Chimp. nCRF';
        category{3} = 'Video 5: Chimp. ExtraNCRF';

    case 6
        
        category{1} = 'Video 6: Chimp. CRF';
        category{2} = 'Video 6: Chimp. nCRF';
        category{3} = 'Video 6: Chimp. ExtraNCRF';

    case 7
        
        category{1} = 'Video 7: Leopard. CRF';
        category{2} = 'Video 7: Leopard. nCRF';
        category{3} = 'Video 7: Leopard. ExtraNCRF';

    case 8
        
        category{1} = 'Video 8: Panda. CRF';
        category{2} = 'Video 8: Panda. nCRF';
        category{3} = 'Video 8: Panda. ExtraNCRF';

    case 9
        
        category{1} = 'Video 9: Javali. CRF';
        category{2} = 'Video 9: Javali. nCRF';
        category{3} = 'Video 9: Javali. ExtraNCRF';

    case 10
        
        category{1} = 'Video 10: Crocodile. CRF';
        category{2} = 'Video 10: Crocodile. nCRF';
        category{3} = 'Video 10: Crocodile. ExtraNCRF';

    case 11
        
        category{1} = 'Video 11: Birds. CRF';
        category{2} = 'Video 11: Birds. nCRF';
        category{3} = 'Video 11: Birds. ExtraNCRF';
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   DEFINE FILE NAME AND PATH   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Define file name and path...');

filepath = filename(date,site_index,channel,video_index,start_time,end_time)

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
%P = nRepetitions;

%Tempos inicial e final do trial
%start_time = 0;
%end_time = 0;

%N?mero de pontos no vetor list
%Q = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   LOAD TRIALS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load trials...');
    
Spass = load(strcat('_',registro,'-','v',int2str(video_index),'.mat'));

M = max(Spass.stimIds);

for i=1:M
   
    trials_label(i).labels = find(Spass.stimIds == i); 

    nTrials(i) = size(Spass.spike_times(trials_label(i).labels(:),:),1);
    
end

Spass.spike_times = Spass.spike_times./ 32000;

P = min(nTrials);
clear nTrials;

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

stawrite(X,filepath);

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