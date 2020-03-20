function compareMetric

bin_size = 1;

registro = importdata('memoryBackwardProtocols.txt');

for i=1:length(registro)
    
    registro{i} = registro{i}(2:end-4);
    
%     if i < 27
%         
%         registro{i} = registro{i}(2:end-7);
%         
%     else
%         
%         registro{i} = registro{i}(2:end-9);
%         
%     end

    if i == 1
        
        end_time = 4000;
        
    else
        
        end_time = 10000;
        
    end
    
    start_time = 0;
    
    filepathWOlatency = strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/metricMovieWOlatency/',int2str(start_time),'-',int2str(end_time),'/');
    
    file_type_ae_WO = 'full-movie-ae-WO-latency';
    file_type_original_WO = 'full-movie-original-WO-latency';
    file_type_invertido_WO = strcat('bin-size-',int2str(bin_size),'-','full-movie-invertido-WO-latency');

    data_ae_WO_latency = load(strcat(filepathWOlatency,char(registro{i}),'/full-movie/full-movie-ae/',char(registro{i}),'-',file_type_ae_WO,'.mat'));
    data_original_WO_latency = load(strcat(filepathWOlatency,char(registro{i}),'/full-movie/full-movie-original/',char(registro{i}),'-',file_type_original_WO,'.mat'));
    data_invertido_WO_latency = load(strcat(filepathWOlatency,char(registro{i}),'/full-movie/full-movie-invertido/',char(registro{i}),'-',file_type_invertido_WO,'.mat'));
    
    filepath = strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/metricMovieWithlatency/',int2str(start_time),'-',int2str(end_time),'/');
     
    file_type_ae = 'full-movie-ae-With-latency';
    file_type_original = 'full-movie-original-With-latency';
    file_type_invertido = strcat('bin-size-',int2str(bin_size),'-','full-movie-invertido-With-latency');
    
    data_ae = load(strcat(filepath,char(registro{i}),'/full-movie/full-movie-ae/',char(registro{i}),'-',file_type_ae,'.mat'));
    data_original = load(strcat(filepath,char(registro{i}),'/full-movie/full-movie-original/',char(registro{i}),'-',file_type_original,'.mat'));
    data_invertido = load(strcat(filepath,char(registro{i}),'/full-movie/full-movie-invertido/',char(registro{i}),'-',file_type_invertido,'.mat'));
    
    ae_WO_info(i) = data_ae_WO_latency.metric_analysis.max_info.max_info_tpmc_T;
    ae_WO_idx(i) = data_ae_WO_latency.metric_analysis.max_info.max_info_tpmc_idx_T;
    
    original_WO_info(i) = data_original_WO_latency.metric_analysis.max_info.max_info_tpmc_T;
    original_WO_idx(i) = data_original_WO_latency.metric_analysis.max_info.max_info_tpmc_idx_T;
    
    invertido_WO_info(i) = data_invertido_WO_latency.metric_analysis.max_info.max_info_tpmc_T;
    invertido_WO_idx(i) = data_invertido_WO_latency.metric_analysis.max_info.max_info_tpmc_idx_T;
  
    ae_info(i) = data_ae.metric_analysis.max_info.max_info_tpmc_T;
    ae_idx(i) = data_ae.metric_analysis.max_info.max_info_tpmc_idx_T;
    
    original_info(i) = data_original.metric_analysis.max_info.max_info_tpmc_T;
    original_idx(i) = data_original.metric_analysis.max_info.max_info_tpmc_idx_T;
    
    invertido_info(i) = data_invertido.metric_analysis.max_info.max_info_tpmc_T;
    invertido_idx(i) = data_invertido.metric_analysis.max_info.max_info_tpmc_idx_T;
    
end

ae_info_aumentou = length(ae_WO_info(ae_WO_info > ae_info));
ae_idx_mudou = length(ae_WO_idx(ae_WO_idx > ae_idx));

original_info_aumentou = length(original_WO_info(original_WO_info > original_info));
original_idx_mudou = length(original_WO_idx(original_WO_idx > original_idx));

invertido_info_aumentou = length(invertido_WO_info(invertido_WO_info > invertido_info));
invertido_idx_mudou = length(invertido_WO_idx(invertido_WO_idx > invertido_idx));

infos_WO = struct('ae_WO_info',ae_WO_info,'ae_WO_idx',ae_WO_idx,'original_WO_info',original_WO_info,'original_WO_idx',original_WO_idx,'invertido_WO_info',invertido_WO_info,'invertido_WO_idx',invertido_WO_idx);

infos = struct('ae_info',ae_info,'ae_idx',ae_idx,'original_info',original_info,'original_idx',original_idx,'invertido_info',invertido_info,'invertido_idx',invertido_idx);

resultsMetric = struct('ae_info_aumentou',ae_info_aumentou,'ae_idx_mudou',ae_idx_mudou,'original_info_aumentou',original_info_aumentou,'original_idx_mudou',original_idx_mudou,'invertido_info_aumentou',invertido_info_aumentou,'invertido_idx_mudou',invertido_idx_mudou,'infos_WO',infos_WO,'infos',infos);

save('/Volumes/Data/DATA/Forward-Backward/Jerome/metricMovie','resultsMetric');

ae = ae_WO_info ./ ae_info;
original = original_WO_info ./ original_info;
invertido = invertido_WO_info ./ invertido_info;

ae_q = ae_WO_idx ./ ae_idx;
original_q = original_WO_idx ./ original_idx;
invertido_q = invertido_WO_idx ./ invertido_idx;

f = figure;
plot(1:length(ae),ae,'r');
hold on;
plot(1:length(ae),original,'b');
hold on;
plot(1:length(ae),invertido,'g');
title('Entropy rate');
legend('Spontaneous Activity', 'Forward-Backward', 'Forward-Backward(reversed)');
print(f,'-depsc',strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/metricMovie','resultsMetric-','entropy-rate','.eps'));

g = figure;
plot(1:length(ae),ae_q,'r');
hold on;
plot(1:length(ae),original_q,'b');
hold on;
plot(1:length(ae),invertido_q,'g');
title('Resolution rate');
legend('Spontaneous Activity', 'Forward-Backward', 'Forward-Backward(reversed)');
print(g,'-depsc',strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/metricMovie','resultsMetric-','resolution-rate','.eps'));



end

