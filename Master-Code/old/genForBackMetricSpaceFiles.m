function genForBackMetricSpaceFiles

cells_data_file = get_all_cells_FB;

for iCell=1:length(cells_data_file)
    
    nVideos = length(cells_data_file(iCell).registro_video);
    
    for iVideo=1:nVideos
    
        nome_do_registro = char(strcat(cells_data_file(iCell).registro_video(iVideo).datafile));
        
        latency_file = strcat('ForBack-Cell','-',int2str(iCell),'-','Video','-',int2str(cells_data_file(iCell).registro_video(iVideo).video),'-',cells_data_file(iCell).registro_video(iVideo).datafile,'-','latency','.mat');
        
        Conditions = [1 2];
        
        prefix_nome_metric_file =  strcat('ForBack-Cell','-',int2str(iCell),'-','Video','-',int2str(cells_data_file(iCell).registro_video(iVideo).video),'-',cells_data_file(iCell).registro_video(iVideo).datafile,'-');
        
        nome_metric_file = strcat(prefix_nome_metric_file,'direto','-','WITHlatency','-','metricSpace');

        inverter_movie = 0;
        WOlatency = false;
        spass2metric(nome_do_registro,inverter_movie,latency_file,WOlatency,Conditions,nome_metric_file);
        
        nome_metric_file = strcat(prefix_nome_metric_file,'direto','-','WOlatency','-','metricSpace');
        
        inverter_movie = 0;
        WOlatency = true;
        spass2metric(nome_do_registro,inverter_movie,latency_file,WOlatency,Conditions,nome_metric_file);
        
        nome_metric_file = strcat(prefix_nome_metric_file,'reverso','-','WITHlatency','-','metricSpace');
        
        inverter_movie = 1;
        WOlatency = false;
        spass2metric(nome_do_registro,inverter_movie,latency_file,WOlatency,Conditions,nome_metric_file);
        
        nome_metric_file = strcat(prefix_nome_metric_file,'reverso','-','WOlatency','-','metricSpace');
        
        inverter_movie = 1;
        WOlatency = true;
        spass2metric(nome_do_registro,inverter_movie,latency_file,WOlatency,Conditions,nome_metric_file);
        
    end
       
end

