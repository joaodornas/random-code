function reliabilityWilcoxon


Protocol(1).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-20/sitio1/E1/v3/full-movie/reliability/v3-12-08-20-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(2).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-27/sitio1/E2/v3/full-movie/reliability/v3-12-08-27-sitio1-E2-reliability-full-movie-bandwidth-15.mat');
Protocol(3).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio1/E1/v3/full-movie/reliability/v3-12-08-29-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(4).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio1/E1/v4/full-movie/reliability/v4-12-08-29-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(5).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio2/E1/v1/full-movie/reliability/v1-12-08-29-sitio2-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(6).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E1/v1/full-movie/reliability/v1-12-08-30-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(7).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E1/v4/full-movie/reliability/v4-12-08-30-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(8).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v1/full-movie/reliability/v1-12-08-30-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
Protocol(9).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v3/full-movie/reliability/v3-12-08-30-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
Protocol(10).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v4/full-movie/reliability/v4-12-08-30-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
Protocol(11).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v1/full-movie/reliability/v1-12-08-31-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(12).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v3/full-movie/reliability/v3-12-08-31-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(13).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v4/full-movie/reliability/v4-12-08-31-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(14).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v1/full-movie/reliability/v1-12-08-31-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
Protocol(15).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v3/full-movie/reliability/v3-12-08-31-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
Protocol(16).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v4/full-movie/reliability/v4-12-08-31-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
Protocol(17).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-03/sitio1/E3/v2/full-movie/reliability/v2-12-09-03-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
Protocol(18).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-03/sitio2/E1/v4/full-movie/reliability/v4-12-09-03-sitio2-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(19).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-04/sitio1/E1/v2/full-movie/reliability/v2-12-09-04-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(20).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-04/sitio1/E1/v4/full-movie/reliability/v4-12-09-04-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(21).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-05/sitio1/E1/v3/full-movie/reliability/v3-12-09-05-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(22).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-05/sitio1/E1/v4/full-movie/reliability/v4-12-09-05-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(23).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-06/sitio1/E1/v3/full-movie/reliability/v3-12-09-06-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(24).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-06/sitio1/E1/v4/full-movie/reliability/v4-12-09-06-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(25).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-10-16/sitio1/E1/v1/full-movie/reliability/v1-12-10-16-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(26).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E1/v311/full-movie/reliability/v311-12-12-17-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(27).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E1/v321/full-movie/reliability/v321-12-12-17-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
Protocol(28).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E3/v331/full-movie/reliability/v331-12-12-17-sitio1-E3-reliability-full-movie-bandwidth-15.mat');


for i=1:size(Protocol,2)
    
    
   diretoDireto(i) = Protocol(i).reliability.reliabilityFullMovieData.DiretoDireto.schreiberData.reliability; 
   
   reversoReverso(i) = Protocol(i).reliability.reliabilityFullMovieData.ReversoReverso.schreiberData.reliability; 
   
   diretoReverso(i) = Protocol(i).reliability.reliabilityFullMovieData.DiretoReverso.schreiberData.reliability;
   
   diretoInvertido(i) = Protocol(i).reliability.reliabilityFullMovieData.DiretoInvertido.schreiberData.reliability;    
    
    
end

hDDireto = lillietest(diretoDireto)

hRReverso = lillietest(reversoReverso)

hDReverso = lillietest(diretoReverso)

hDInvertido = lillietest(diretoInvertido)

[p,h] = signrank(diretoDireto,reversoReverso)

[p,h] = signrank(diretoDireto,diretoReverso)

[p,h] = signrank(diretoDireto,diretoInvertido)

[p,h] = signrank(diretoReverso,0.55)
[p,h] = signrank(diretoReverso,0.49)
[p,h] = signrank(diretoReverso,0.42)
[p,h] = signrank(diretoReverso,0.36)

diretoDireto(isnan(diretoDireto)) = 0;
reversoReverso(isnan(reversoReverso)) = 0;
diretoReverso(isnan(diretoReverso)) = 0;
diretoInvertido(isnan(diretoInvertido)) = 0;

median(diretoDireto)
median(reversoReverso)
median(diretoReverso)
median(diretoInvertido)


end

