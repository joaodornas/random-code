function plotMetricResults(nome_do_registro)

disp('BEGIN - plotMetricResults');

%%%   LE DADOS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataset = nome_do_registro;

data = load(strcat(dataset,'.mat'));

X = data.metric_analysis.X;
opts = data.metric_analysis.opts;
out_T = data.metric_analysis.out_T;
out_I = data.metric_analysis.out_I;
max_info_plugin_idx_T = data.metric_analysis.max_info.max_info_plugin_idx_T;
max_info_plugin_idx_I = data.metric_analysis.max_info.max_info_plugin_idx_I;
info_plugin_T = data.metric_analysis.info.info_plugin_T;
info_tpmc_T = data.metric_analysis.info.info_tpmc_T;
info_jack_T = data.metric_analysis.info.info_jack_T;
info_plugin_I = data.metric_analysis.info.info_plugin_I;
info_tpmc_I = data.metric_analysis.info.info_tpmc_I;
info_jack_I = data.metric_analysis.info.info_jack_I;
info_unjk = data.metric_analysis.info.info_unjk;
info_jk_sem = data.metric_analysis.info.info_jk_sem;
HBias_T = data.metric_analysis.max_info.HBias_T;
HBias_std_T = data.metric_analysis.max_info.HBias_std_T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% RASTER PLOT    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure;
set(gcf,'name',['Metric ' dataset]); 

staraster(X,[opts.start_time opts.end_time]);
title('Raster plot');

print(f1,'-depsc',strcat(dataset,'-raster-plot'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DISTANCE-MATRIX TIMING PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2 = figure;
imagesc(out_T(max_info_plugin_idx_T).d);
xlabel('Spike train index');
ylabel('Spike train index');
title('Distance matrix at maximum information - timing');
print(f2,'-depsc',strcat(dataset,'-distance-matrix-timing'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CONFUSION-MATRIX TIMING PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f3 = figure;
imagesc(out_T(max_info_plugin_idx_T).cm);
xlabel('Category index');
ylabel('Category index');
title('Confusion matrix from clustering of distances - timing');
print(f3,'-depsc',strcat(dataset,'-confusion-matrix-timing'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DISTANCE-MATRIX INTERVAL PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f4 = figure;
imagesc(out_I(max_info_plugin_idx_I).d);
xlabel('Spike train index');
ylabel('Spike train index');
title('Distance matrix at maximum information - interval');
print(f4,'-depsc',strcat(dataset,'-distance-matrix-interval'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CONFUSION-MATRIX INTERVAL PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f5 = figure;
imagesc(out_I(max_info_plugin_idx_I).cm);
xlabel('Category index');
ylabel('Category index');
title('Confusion matrix from clustering of distances - interval');
print(f5,'-depsc',strcat(dataset,'-confusion-matrix-interval'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INFO TIMING-INTERVAL PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
print(f6,'-depsc',strcat(dataset,'-plot-info-timing-interval'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INFO + SHUFFLE INFO PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
print(f7,'-depsc',strcat(dataset,'-plot-info-shuffle'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INFO + SHUFFLE INFO PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
print(f8,'-depsc',strcat(dataset,'-plot-info-shuffle-resolution'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END - plotMetricResults');

end

