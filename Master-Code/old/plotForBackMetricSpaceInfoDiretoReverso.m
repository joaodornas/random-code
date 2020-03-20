
function plotForBackMetricSpaceInfoDiretoReverso

cells_data_file = get_all_cells_FB;

prefix = '/Users/joaodornas/Dropbox/_Research/_ANALYSIS/_PAPERS/Master of Science/Forward-Backward/metricspace/infofile';

for iCell=1:length(cells_data_file)
    
    nVideos = length(cells_data_file(iCell).registro_video);
    
    for iVideo=1:nVideos
        
        prefix_nome_metric_file =  strcat('ForBack-Cell','-',int2str(iCell),'-','Video','-',int2str(cells_data_file(iCell).registro_video(iVideo).video),'-',cells_data_file(iCell).registro_video(iVideo).datafile,'-');
        
        name = strcat('ForBack-Cell','-',int2str(iCell),'-','Video','-',int2str(cells_data_file(iCell).registro_video(iVideo).video));

        nome_metric_file = strcat('direto','-','WITHlatency','-','metricSpace');
        
        info_file = strcat(prefix,'/',prefix_nome_metric_file,nome_metric_file,'.mat');
        
        dataset_direto_withlatency = load(info_file);
        
        nome_metric_file = strcat('direto','-','WOlatency','-','metricSpace');
        
        info_file = strcat(prefix,'/',prefix_nome_metric_file,nome_metric_file,'.mat');
        
        dataset_direto_wolatency = load(info_file);

        nome_metric_file = strcat('reverso','-','WITHlatency','-','metricSpace');
        
        info_file = strcat(prefix,'/',prefix_nome_metric_file,nome_metric_file,'.mat');
        
        dataset_reverso_withlatency = load(info_file);
        
        nome_metric_file = strcat('reverso','-','WOlatency','-','metricSpace');
        
        info_file = strcat(prefix,'/',prefix_nome_metric_file,nome_metric_file,'.mat');
        
        dataset_reverso_wolatency = load(info_file);
        
        plotForBack(dataset_direto_withlatency,dataset_reverso_withlatency,'WITH-latency',prefix_nome_metric_file,name);
        
        plotForBack(dataset_direto_wolatency,dataset_reverso_wolatency,'WO-latency',prefix_nome_metric_file,name);
        
        close all

    end

end

end

function plotForBack(direto,reverso,kind,prefix,name)

info_T_direto = direto.metric_analysis.all_info.info_tpmc_T;
opts_direto = direto.metric_analysis.opts;

info_T_reverso = reverso.metric_analysis.all_info.info_tpmc_T;
opts_reverso = reverso.metric_analysis.opts;

info_shuffle_direto = direto.metric_analysis.all_info.HBias_T;
info_shuffle_reverso = reverso.metric_analysis.all_info.HBias_T;

disp('PLOTTING INFO TIMING-TPMC-SHUFFLE-RESOLUTION');
% INFO + SHUFFLE INFO PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure;
%resolution = 1 ./ opts_direto.shift_cost ;

plot(1:length(opts_direto.shift_cost),info_T_direto,'b-o');
hold on;
plot(1:length(opts_reverso.shift_cost),info_T_reverso,'b--o');
plot(1:length(opts_reverso.shift_cost),info_shuffle_direto,'r-');
plot(1:length(opts_reverso.shift_cost),info_shuffle_reverso,'r--');
set(gca,'xtick',1:length(opts_direto.shift_cost));
set(gca,'xticklabel',opts_direto.shift_cost);
set(gca,'xlim',[1 length(opts_direto.shift_cost)]);
set(gca,'ylim',[0 1]);
xlabel('Shift-Cost (1/s)');
ylabel('Information (bits)');
legend('For-Back (without inversion)','For-Back (with inversion)','Shuffle (without inversion)','Shuffle (with inversion)','location','best');
title(strcat(name,'-',kind));

scalefig(gcf,1.5);
print(f,'-depsc',strcat(prefix,kind,'-plot-info-diretoXreverso-resolution.eps'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

