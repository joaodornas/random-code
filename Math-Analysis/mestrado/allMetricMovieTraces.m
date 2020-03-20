function allMetricMovieTraces(intervals)


date = {'12-08-20' '12-08-27' '12-08-27' '12-08-29' '12-08-29' '12-08-29' '12-08-30' '12-08-30' '12-08-30' '12-08-30' '12-08-30' '12-08-31' '12-08-31' '12-08-31' '12-08-31' '12-08-31' '12-08-31' '12-09-03' '12-09-03' '12-09-04' '12-09-04' '12-09-05' '12-09-05' '12-09-06' '12-09-06' '12-10-16' '12-12-17' '12-12-17' '12-12-17'};

sitio = {int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(2) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(2) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1) int2str(1)};

channel = {'E1' 'E1' 'E2' 'E1' 'E1' 'E1' 'E1' 'E3' 'E1' 'E3' 'E3' 'E1' 'E3' 'E1' 'E3' 'E1' 'E3' 'E3' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E1' 'E3'};

video = {int2str(3) int2str(3) int2str(3) int2str(3) int2str(4) int2str(1) int2str(4) int2str(4) int2str(1) int2str(1) int2str(3) int2str(1) int2str(1) int2str(3) int2str(3) int2str(4) int2str(4) int2str(2) int2str(4) int2str(4) int2str(2) int2str(3) int2str(4) int2str(3) int2str(4) int2str(1) int2str(311) int2str(321) int2str(331)};

label = {'nsp003a01_1b' 'nsp005a01_1b' 'nsp005a01_2b' 'nsp006a01_1b' 'nsp006a02_1b' 'nsp006b01_1b' 'nsp007a01_1b' 'nsp007a01_2b' 'nsp007a02_1b' 'nsp007a02_2b' 'nsp007a03_2a' 'nsp008a01_1b' 'nsp008a01_2b' 'nsp008a02_1a' 'nsp008a02_2b' 'nsp008a03_1b' 'nsp008a03_2a' 'nsp009a01_2b' 'nsp009b01_1a' 'nsp010a02_1b' 'nsp010a1_1b' 'nsp011a01_1a' 'nsp011a02_1b' 'nsp012a01_1b' 'nsp012a02_1a' 'nps013a01_1b' 'nsp033a09_1b' 'nsp033a09_1c' 'nsp033a09_3b'};

end_time = [4000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000];
 

for c=1:length(date)
    
    s = round( (end_time(c)/1000 ) / (intervals/1000) );
    
    for i=0:(s-1)

        disp(label(c));
        
        movie = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/',date(c),'/sitio',sitio(c),'/',channel(c),'/v',video(c),'/full-movie/trecho-invertido/v',video(c),'-',date(c),'-sitio',sitio(c),'-',channel(c),'-',int2str(intervals*i),'-',int2str(intervals*i+intervals),'-90-bin-size-1-trecho-invertido.stam');
        
        metricMovie(char(movie),intervals*i,intervals*i+intervals,60);

        close all;

    end

end

end


