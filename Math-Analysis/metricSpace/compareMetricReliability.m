function compareMetricReliability


registro = importdata('memoryBackwardProtocols.txt');

Qmax = [0 0.0625 0.125 0.25 0.5 1 2 4 8 16 32 64 128 256 512];

for i=1:length(registro)

    if i == 1
        
        end_time = 4000;
        
    else
        
        end_time = 10000;
        
    end
    
    if i < 27
        
        protocol = registro{i}(2:end-7);
        
    else
        
        protocol = registro{i}(2:end-9);
        
    end
    
    name = registro{i}(2:end-4);
    
    start_time = 0;
    
    filepath = strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/metricMovieWOlatency/',int2str(start_time),'-',int2str(end_time),'/',name);

    file_type_original = 'full-movie-original-WO-latency.mat';

    metricFile = strcat(filepath,'/full-movie/full-movie-original/',name,'-',file_type_original);
    
    metricData = load(metricFile);
    
    getDataFor = load(strcat(protocol,'-reliabilityOptimal-Forward.mat'));
    
    getDataBack = load(strcat(protocol,'-reliabilityOptimal-Backward.mat'));
    
    results(i).name = getDataFor.data.name;
    
    results(i).max_for_real = max(getDataFor.data.max_for_real);
    
    results(i).max_back_real = max(getDataBack.data.max_back_real);
    
    results(i).max_resolution_for = max(getDataFor.data.max_resolution_for);
    
    results(i).max_resolution_back = max(getDataBack.data.max_resolution_back);
    
    results(i).max_info_timing = metricData.metric_analysis.max_info.max_info_tpmc_T;
    
    results(i).max_info_resolution = 1./Qmax(metricData.metric_analysis.max_info.max_info_tpmc_idx_T); 
    
end

metricReliability = results;

save(strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/metricMovie/metricReliability'),'metricReliability');


end

