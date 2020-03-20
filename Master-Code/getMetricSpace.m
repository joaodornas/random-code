
function getMetricSpace(spassregistroName,startTime,endTime,nRepetitions)

ct = cputime;

idx = strfind(spassregistroName,'\');

registro = spassregistroName(idx(length(idx))+1:length(spassregistroName));

disp('READING DATA');
% READ DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%
X = staread(strrep(strcat(registro(1:length(registro)-4),'.stam'),'/',filesep));
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

opts.startTime = startTime / 1000;
opts.endTime = endTime / 1000;
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

all_info = struct('info_plugin_T',info_plugin_T,'info_plugin_I',info_plugin_I,'info_tpmc_T',info_tpmc_T,'info_tpmc_I',info_tpmc_I,'info_jack_T',info_jack_T,'info_jack_I',info_jack_I,'info_unjk',info_unjk,'HShuffle_T',HShuffle_T,'HShuffle_I',HShuffle_I,'info_jk_sem',info_jk_sem,'HBias_T',HBias_T,'HBias_I',HBias_I,'HBias_std_T',HBias_std_T,'HBias_std_I',HBias_std_I);

max_info = struct('max_info_plugin_idx_T',max_info_plugin_idx_T,'max_info_plugin_T',max_info_plugin_T,'max_info_plugin_idx_I',max_info_plugin_idx_I,'max_info_plugin_I',max_info_plugin_I,'Hcount_plugin_info_T',Hcount_plugin_T,'Hcount_plugin_info_I',Hcount_plugin_I,'max_info_tpmc_idx_T',max_info_tpmc_idx_T,'max_info_tpmc_T',max_info_tpmc_T,'max_info_tpmc_idx_I',max_info_tpmc_idx_I,'max_info_tpmc_I',max_info_tpmc_I,'Hcount_tpmc_info_T',Hcount_tpmc_T,'Hcount_tpmc_info_I',Hcount_tpmc_I,'max_info_jack_idx_T',max_info_jack_idx_T,'max_info_jack_T',max_info_jack_T,'max_info_jack_idx_I',max_info_jack_idx_I,'max_info_jack_I',max_info_jack_I,'Hcount_jack_info_T',Hcount_jack_T,'Hcount_jack_info_I',Hcount_jack_I,'HBias_T',HBias_T,'HBias_std_T',HBias_std_T,'HBias_I',HBias_I,'HBias_std_I',HBias_std_I);

metric_analysis = struct('opts',opts,'X', X, 'out_T', out_T, 'out_I', out_I, 'all_info', all_info, 'max_info', max_info, 'Y_I', Y_I, 'Y_T', Y_T, 'SHUF_I', SHUF_I, 'SHUF_T', SHUF_T);

registro = registro(1:length(registro)-5);
save(registro,'metric_analysis');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING RASTER');
% RASTER PLOT    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure;
set(gcf,'name',['Metric ' registro]); 

staraster(X,[opts.startTime opts.endTime]);
title('Raster plot');

print(f1,'-depsc',strcat(registro,'-raster-plot'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING DISTANCE-MATRIX TIMING');
% DISTANCE-MATRIX TIMING PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2 = figure;
imagesc(out_T(max_info_plugin_idx_T).d);
xlabel('Spike train index');
ylabel('Spike train index');
title('Distance matrix at maximum information - timing');
print(f2,'-depsc',strcat(registro,'-distance-matrix-timing'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING CONFUSION-MATRIX TIMING');
% CONFUSION-MATRIX TIMING PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3 = figure;
imagesc(out_T(max_info_plugin_idx_T).cm);
xlabel('Category index');
ylabel('Category index');
title('Confusion matrix from clustering of distances - timing');
print(f3,'-depsc',strcat(registro,'-confusion-matrix-timing'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING DISTANCE-MATRIX INTERVAL');
% DISTANCE-MATRIX INTERVAL PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%
f4 = figure;
imagesc(out_I(max_info_plugin_idx_I).d);
xlabel('Spike train index');
ylabel('Spike train index');
title('Distance matrix at maximum information - interval');
print(f4,'-depsc',strcat(registro,'-distance-matrix-interval'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING CONFUSION-MATRIX INTERVAL');
% CONFUSION-MATRIX INTERVAL PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f5 = figure;
imagesc(out_I(max_info_plugin_idx_I).cm);
xlabel('Category index');
ylabel('Category index');
title('Confusion matrix from clustering of distances - interval');
print(f5,'-depsc',strcat(registro,'-confusion-matrix-interval'));
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
print(f6,'-depsc',strcat(registro,'-plot-info-timing-interval'));
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
set(gca,'ylim',[0 1]);
xlabel('Temporal precision (1/sec)');
ylabel('Information (bits)');
legend('Original data (\pm 2 SE via Jackknife)','Shuffled data (\pm 2 SD)',...
       'location','best');

scalefig(gcf,1.5);
print(f7,'-depsc',strcat(registro,'-plot-info-shuffle'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING INFO TIMING-JACK-SHUFFLE-RESOLUTION');
% INFO + SHUFFLE INFO PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f8 = figure;
resolution = 1 ./ opts.shift_cost ;

errorbar(1:length(opts.shift_cost),info_unjk,2*info_jk_sem);
hold on;
errorbar(1:length(opts.shift_cost),HBias_T,2*HBias_std_T,'r');
hold off;
set(gca,'xtick',1:length(opts.shift_cost));
set(gca,'xticklabel',resolution);
set(gca,'xlim',[min(resolution) max(resolution)]);
set(gca,'ylim',[0 1]);
xlabel('Resolution (s)');
ylabel('Information (bits)');
legend('Original data (\pm 2 SE via Jackknife)','Shuffled data (\pm 2 SD)',...
       'location','best');

scalefig(gcf,1.5);
print(f8,'-depsc',strcat(registro,'-plot-info-shuffle-resolution'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cte = cputime - ct;

disp(strcat('CPUTIME: ',num2str(cte)));

end

