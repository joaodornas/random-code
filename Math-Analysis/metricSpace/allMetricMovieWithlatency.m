function allMetricMovieWithlatency(start_time,end_time)

trials = 60;

bin_size = 1;

registro = importdata('memoryBackwardProtocols.txt');

file_type_ae = 'full-movie-ae-With-latency';
file_type_invertido = strcat('bin-size-',int2str(bin_size),'-','full-movie-invertido-With-latency');
file_type_original = 'full-movie-original-With-latency';
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

registro{1} = registro{1}(2:end-4);

metricMovie(strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/metricMovieWithlatency/',int2str(0),'-',int2str(4000),'/',char(registro{1}),'/full-movie/full-movie-ae/',char(registro{1}),'-',file_type_ae,'.stam'),start_time,4000,40);
close all;

metricMovie(strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/metricMovieWithlatency/',int2str(0),'-',int2str(4000),'/',char(registro{1}),'/full-movie/full-movie-invertido/',char(registro{1}),'-',file_type_invertido,'.stam'),start_time,4000,40);
close all;

metricMovie(strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/metricMovieWithlatency/',int2str(0),'-',int2str(4000),'/',char(registro{1}),'/full-movie/full-movie-original/',char(registro{1}),'-',file_type_original,'.stam'),start_time,4000,40);
close all;

for i=2:length(registro)
    
    registro{i} = registro{i}(2:end-4);
    
    metricMovie(strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/metricMovieWithlatency/',int2str(start_time),'-',int2str(end_time),'/',char(registro{i}),'/full-movie/full-movie-ae/',char(registro{i}),'-',file_type_ae,'.stam'),start_time,end_time,trials);
    close all;
    
    metricMovie(strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/metricMovieWithlatency/',int2str(start_time),'-',int2str(end_time),'/',char(registro{i}),'/full-movie/full-movie-invertido/',char(registro{i}),'-',file_type_invertido,'.stam'),start_time,end_time,trials);
    close all;
    
    metricMovie(strcat('/Volumes/Data/DATA/Forward-Backward/Jerome/metricMovieWithlatency/',int2str(start_time),'-',int2str(end_time),'/',char(registro{i}),'/full-movie/full-movie-original/',char(registro{i}),'-',file_type_original,'.stam'),start_time,end_time,trials);
    close all;

end


end

