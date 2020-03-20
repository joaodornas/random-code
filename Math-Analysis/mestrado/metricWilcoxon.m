function metricWillcoxon

Protocol(1).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-20/sitio1/E1/v3/full-movie/full-movie-original/v3-WOSP-20-08-12-sitio1-E1-full-movie-original.mat');
Protocol(2).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-27/sitio1/E2/v3/full-movie/full-movie-original/v3-WOSP-27-08-12-sitio1-E2-full-movie-original.mat');
Protocol(3).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio1/E1/v3/full-movie/full-movie-original/v3-WOSP-29-08-12-sitio1-E1-full-movie-original.mat');
Protocol(4).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio1/E1/v4/full-movie/full-movie-original/v4-WOSP-29-08-12-sitio1-E1-full-movie-original.mat');
Protocol(5).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio2/E1/v1/full-movie/full-movie-original/v1-WOSP-29-08-12-sitio2-E1-full-movie-original.mat');
Protocol(6).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E1/v1/full-movie/full-movie-original/v1-WOSP-30-08-12-sitio1-E1-full-movie-original.mat');
Protocol(7).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E1/v4/full-movie/full-movie-original/v4-WOSP-30-08-12-sitio1-E1-full-movie-original.mat');
Protocol(8).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v1/full-movie/full-movie-original/v1-WOSP-30-08-12-sitio1-E3-full-movie-original.mat');
Protocol(9).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v3/full-movie/full-movie-original/v3-WOSP-30-08-12-sitio1-E3-full-movie-original.mat');
Protocol(10).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v4/full-movie/full-movie-original/v4-WOSP-30-08-12-sitio1-E3-full-movie-original.mat');
Protocol(11).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v1/full-movie/full-movie-original/v1-WOSP-31-08-12-sitio1-E1-full-movie-original.mat');
Protocol(12).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v3/full-movie/full-movie-original/v3-WOSP-31-08-12-sitio1-E1-full-movie-original.mat');
Protocol(13).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v4/full-movie/full-movie-original/v4-WOSP-31-08-12-sitio1-E1-full-movie-original.mat');
Protocol(14).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v1/full-movie/full-movie-original/v1-WOSP-31-08-12-sitio1-E3-full-movie-original.mat');
Protocol(15).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v3/full-movie/full-movie-original/v3-WOSP-31-08-12-sitio1-E3-full-movie-original.mat');
Protocol(16).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v4/full-movie/full-movie-original/v4-WOSP-31-08-12-sitio1-E3-full-movie-original.mat');
Protocol(17).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-03/sitio1/E3/v2/full-movie/full-movie-original/v2-WOSP-12-09-03-sitio1-E3-full-movie-original.mat');
Protocol(18).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-03/sitio2/E1/v4/full-movie/full-movie-original/v4-WOSP-12-09-03-sitio2-E1-full-movie-original.mat');
Protocol(19).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-04/sitio1/E1/v2/full-movie/full-movie-original/v2-WOSP-12-09-04-sitio1-E1-full-movie-original.mat');
Protocol(20).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-04/sitio1/E1/v4/full-movie/full-movie-original/v4-WOSP-12-09-04-sitio1-E1-full-movie-original.mat');
Protocol(21).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-05/sitio1/E1/v3/full-movie/full-movie-original/v3-WOSP-12-09-05-sitio1-E1-full-movie-original.mat');
Protocol(22).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-05/sitio1/E1/v4/full-movie/full-movie-original/v4-WOSP-12-09-05-sitio1-E1-full-movie-original.mat');
Protocol(23).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-06/sitio1/E1/v3/full-movie/full-movie-original/v3-WOSP-12-09-06-sitio1-E1-full-movie-original.mat');
Protocol(24).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-06/sitio1/E1/v4/full-movie/full-movie-original/v4-WOSP-12-09-06-sitio1-E1-full-movie-original.mat');
Protocol(25).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-10-16/sitio1/E1/v1/full-movie/full-movie-original/v1-WOSP-12-10-16-sitio1-E1-full-movie-original.mat');
Protocol(26).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E1/v311/full-movie/full-movie-original/v311-WOSP-12-12-17-sitio1-E1-full-movie-original.mat');
Protocol(27).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E1/v321/full-movie/full-movie-original/v321-WOSP-12-12-17-sitio1-E1-full-movie-original.mat');
Protocol(28).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E3/v331/full-movie/full-movie-original/v331-WOSP-12-12-17-sitio1-E3-full-movie-original.mat');


Protocol(1).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-20/sitio1/E1/v3/full-movie/full-movie-invertido/v3-WOSP-12-08-20-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(2).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-27/sitio1/E2/v3/full-movie/full-movie-invertido/v3-WOSP-12-08-27-sitio1-E2-full-movie-invertido-bin-size-1.mat');
Protocol(3).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio1/E1/v3/full-movie/full-movie-invertido/v3-WOSP-12-08-29-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(4).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio1/E1/v4/full-movie/full-movie-invertido/v4-WOSP-12-08-29-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(5).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio2/E1/v1/full-movie/full-movie-invertido/v1-WOSP-12-08-29-sitio2-E1-full-movie-invertido-bin-size-1.mat');
Protocol(6).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E1/v1/full-movie/full-movie-invertido/v1-WOSP-12-08-30-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(7).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E1/v4/full-movie/full-movie-invertido/v4-WOSP-12-08-30-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(8).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v1/full-movie/full-movie-invertido/v1-WOSP-12-08-30-sitio1-E3-full-movie-invertido-bin-size-1.mat');
Protocol(9).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v3/full-movie/full-movie-invertido/v3-WOSP-12-08-30-sitio1-E3-full-movie-invertido-bin-size-1.mat');
Protocol(10).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v4/full-movie/full-movie-invertido/v4-WOSP-12-08-30-sitio1-E3-full-movie-invertido-bin-size-1.mat');
Protocol(11).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v1/full-movie/full-movie-invertido/v1-WOSP-12-08-31-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(12).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v3/full-movie/full-movie-invertido/v3-WOSP-12-08-31-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(13).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v4/full-movie/full-movie-invertido/v4-WOSP-12-08-31-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(14).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v1/full-movie/full-movie-invertido/v1-WOSP-12-08-31-sitio1-E3-full-movie-invertido-bin-size-1.mat');
Protocol(15).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v3/full-movie/full-movie-invertido/v3-WOSP-12-08-31-sitio1-E3-full-movie-invertido-bin-size-1.mat');
Protocol(16).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v4/full-movie/full-movie-invertido/v4-WOSP-12-08-31-sitio1-E3-full-movie-invertido-bin-size-1.mat');
Protocol(17).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-03/sitio1/E3/v2/full-movie/full-movie-invertido/v2-WOSP-12-09-03-sitio1-E3-full-movie-invertido-bin-size-1.mat');
Protocol(18).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-03/sitio2/E1/v4/full-movie/full-movie-invertido/v4-WOSP-12-09-03-sitio2-E1-full-movie-invertido-bin-size-1.mat');
Protocol(19).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-04/sitio1/E1/v2/full-movie/full-movie-invertido/v2-WOSP-12-09-04-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(20).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-04/sitio1/E1/v4/full-movie/full-movie-invertido/v4-WOSP-12-09-04-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(21).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-05/sitio1/E1/v3/full-movie/full-movie-invertido/v3-WOSP-12-09-05-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(22).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-05/sitio1/E1/v4/full-movie/full-movie-invertido/v4-WOSP-12-09-05-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(23).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-06/sitio1/E1/v3/full-movie/full-movie-invertido/v3-WOSP-12-09-06-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(24).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-06/sitio1/E1/v4/full-movie/full-movie-invertido/v4-WOSP-12-09-06-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(25).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-10-16/sitio1/E1/v1/full-movie/full-movie-invertido/v1-WOSP-12-10-16-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(26).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E1/v311/full-movie/full-movie-invertido/v311-WOSP-12-12-17-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(27).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E1/v321/full-movie/full-movie-invertido/v321-WOSP-12-12-17-sitio1-E1-full-movie-invertido-bin-size-1.mat');
Protocol(28).metricIn = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E3/v331/full-movie/full-movie-invertido/v331-WOSP-12-12-17-sitio1-E3-full-movie-invertido-bin-size-1.mat');


for i=1:size(Protocol,2)
   
    
    HTmaxIn(i) = Protocol(i).metricIn.metric_analysis.max_info.max_info_plugin_T;    
    
    HTmax(i) = Protocol(i).metric.metric_analysis.max_info.max_info_plugin_T; 
    
    
end


hin = lillietest(HTmaxIn)
h = lillietest(HTmax)


[p,h] = signrank(HTmax,HTmaxIn)

median(HTmaxIn)

median(HTmax)

end

