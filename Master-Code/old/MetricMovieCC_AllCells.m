function MetricMovieCC_AllCells

cells_data_file = get_all_cells_CC;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prefix = '/Users/joaodornas/Dropbox/_Research/_ANALYSIS/_PAPERS/Master of Science/Center-Surround/metricspace/infofile';

for iCell=1:length(cells_data_file)
    
    nVideos = length(cells_data_file(iCell).registro_video);
    
    for iVideo=1:nVideos
        
        prefix_nome_metric_file =  strcat('CC-Cell','-',int2str(iCell),'-','Video','-',int2str(cells_data_file(iCell).registro_video(iVideo).video),'-',cells_data_file(iCell).registro_video(iVideo).datafile,'-');
        
        nTrials = 60;
        start_time = 500;
        end_time = 9500;
        
        nome_metric_file = strcat('direto','-','WITHlatency','-','metricSpace');
        
        info_file = strcat(prefix,'/',prefix_nome_metric_file,nome_metric_file,'.stam');
        
        metricMovie(info_file,start_time,end_time,nTrials);
        
        close 
        
        nome_metric_file = strcat('direto','-','WOlatency','-','metricSpace');
        
        info_file = strcat(prefix,'/',prefix_nome_metric_file,nome_metric_file,'.stam');
        
        metricMovie(info_file,start_time,end_time,nTrials);
        
        close all    
        
    end

end

