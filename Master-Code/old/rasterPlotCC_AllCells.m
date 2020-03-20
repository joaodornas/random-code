function rasterPlotCC_AllCells

cells_data_file = get_all_cells_CC;

nConditions = 3;

for c=1:length(cells_data_file)
    
    nVideos = length(cells_data_file(c).registro_video);
    
    for v=1:nVideos
    
        registro = char(strcat(cells_data_file(c).registro_video(v).datafile,'.mat'));
        
        name = char(strcat('CC-Cell','-',int2str(c),'-','Video','-',int2str(cells_data_file(c).registro_video(v).video),'-',cells_data_file(c).registro_video(v).datafile));
        
        latency_file = strcat('CC-Cell','-',int2str(c),'-','Video','-',int2str(cells_data_file(c).registro_video(v).video),'-',cells_data_file(c).registro_video(v).datafile,'-','latency','.mat');
    
        rasterAll(registro,nConditions,name,latency_file,true);
        rasterAll(registro,nConditions,name,latency_file,false);
        
        close all
        
    end
    
end

end

