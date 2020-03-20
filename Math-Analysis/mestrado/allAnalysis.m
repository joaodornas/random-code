function allAnalysis

% Protocol(1).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-20/sitio1/E1/v3/full-movie/reliability/v3-12-08-20-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(2).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-27/sitio1/E2/v3/full-movie/reliability/v3-12-08-27-sitio1-E2-reliability-full-movie-bandwidth-15.mat');
% Protocol(3).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio1/E1/v3/full-movie/reliability/v3-12-08-29-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(4).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio1/E1/v4/full-movie/reliability/v4-12-08-29-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(5).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio2/E1/v1/full-movie/reliability/v1-12-08-29-sitio2-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(6).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E1/v1/full-movie/reliability/v1-12-08-30-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(7).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E1/v4/full-movie/reliability/v4-12-08-30-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(8).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v1/full-movie/reliability/v1-12-08-30-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
% Protocol(9).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v3/full-movie/reliability/v3-12-08-30-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
% Protocol(10).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v4/full-movie/reliability/v4-12-08-30-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
% Protocol(11).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v1/full-movie/reliability/v1-12-08-31-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(12).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v3/full-movie/reliability/v3-12-08-31-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(13).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v4/full-movie/reliability/v4-12-08-31-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(14).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v1/full-movie/reliability/v1-12-08-31-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
% Protocol(15).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v3/full-movie/reliability/v3-12-08-31-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
% Protocol(16).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v4/full-movie/reliability/v4-12-08-31-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
% Protocol(17).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-03/sitio1/E3/v2/full-movie/reliability/v2-12-09-03-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
% Protocol(18).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-03/sitio2/E1/v4/full-movie/reliability/v4-12-09-03-sitio2-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(19).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-04/sitio1/E1/v2/full-movie/reliability/v2-12-09-04-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(20).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-04/sitio1/E1/v4/full-movie/reliability/v4-12-09-04-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(21).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-05/sitio1/E1/v3/full-movie/reliability/v3-12-09-05-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(22).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-05/sitio1/E1/v4/full-movie/reliability/v4-12-09-05-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(23).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-06/sitio1/E1/v3/full-movie/reliability/v3-12-09-06-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(24).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-06/sitio1/E1/v4/full-movie/reliability/v4-12-09-06-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(25).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-10-16/sitio1/E1/v1/full-movie/reliability/v1-12-10-16-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(26).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E1/v311/full-movie/reliability/v311-12-12-17-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(27).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E1/v321/full-movie/reliability/v321-12-12-17-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocol(28).reliability = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E3/v331/full-movie/reliability/v331-12-12-17-sitio1-E3-reliability-full-movie-bandwidth-15.mat');


% Protocolo(1).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v1/reliability/v1-12-10-22-sitio1-E2-reliability-full-movie-bandwidth-15.mat');
% Protocolo(2).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v3/reliability/v3-12-10-22-sitio1-E2-reliability-full-movie-bandwidth-15.mat');
% Protocolo(3).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v2/reliability/v2-12-10-23-sitio1-E2-reliability-full-movie-bandwidth-15.mat');
% Protocolo(4).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v4/reliability/v4-12-10-23-sitio1-E2-reliability-full-movie-bandwidth-15.mat');
% Protocolo(5).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E1/v3/reliability/v3-12-10-24-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(6).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E2/v2/reliability/v2-12-10-24-sitio1-E2-reliability-full-movie-bandwidth-15.mat');
% Protocolo(7).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-26/sitio2/E1/v3/reliability/v3-12-10-26-sitio2-E1-reliability-full-movie-bandwidth-15.mat');
% % load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v5/reliability/v5-12-10-31-sitio1-E3-reliability-full-movie-bandwidth-15.mat')
% % load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v6/reliability/v6-12-10-31-sitio1-E3-reliability-full-movie-bandwidth-15.mat')
% Protocolo(8).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v5/reliability/v5-12-11-01-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(9).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v6/reliability/v6-12-11-01-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(10).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v5/reliability/v5-12-11-06-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(11).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v6/reliability/v6-12-11-06-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(12).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v5/reliability/v5-12-11-07-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(13).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v6/reliability/v6-12-11-07-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(14).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v1/reliability/v1-12-11-08-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(15).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v4/reliability/v4-12-11-08-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(16).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v3/reliability/v3-12-11-09-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(17).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v6/reliability/v6-12-11-09-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(18).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v1/reliability/v1-12-11-13-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(19).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v2/reliability/v2-12-11-13-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(20).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v8/reliability/v8-12-12-03-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
% Protocolo(21).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v11/reliability/v11-12-12-03-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
% Protocolo(22).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v9/reliability/v9-12-12-04-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(23).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v10/reliability/v10-12-12-04-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(24).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-05/sitio1/E1/v6/reliability/v6-12-12-05-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(25).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-05/sitio1/E1/v8/reliability/v8-12-12-05-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(26).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-07/sitio1/E1/v3/reliability/v3-12-12-07-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(27).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-07/sitio1/E1/v8/reliability/v8-12-12-07-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(28).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-10/sitio1/E2/v5/reliability/v5-12-12-10-sitio1-E2-reliability-full-movie-bandwidth-15.mat');
% Protocolo(29).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-10/sitio1/E2/v8/reliability/v8-12-12-10-sitio1-E2-reliability-full-movie-bandwidth-15.mat');
% Protocolo(30).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-12/sitio2/E1/v3/reliability/v3-12-12-12-sitio2-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(31).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-12/sitio2/E1/v4/reliability/v4-12-12-12-sitio2-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(32).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-13/sitio1/E1/v2/reliability/v2-12-12-13-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(33).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-13/sitio1/E1/v8/reliability/v8-12-12-13-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(34).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v511/reliability/v511-12-12-17-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(35).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v512/reliability/v512-12-12-17-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(36).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v521/reliability/v521-12-12-17-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(37).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v811/reliability/v811-12-12-17-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(38).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v812/reliability/v812-12-12-17-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(39).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v821/reliability/v821-12-12-17-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(40).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v822/reliability/v822-12-12-17-sitio1-E1-reliability-full-movie-bandwidth-15.mat');
% Protocolo(41).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v531/reliability/v531-12-12-17-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
% Protocolo(42).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v532/reliability/v532-12-12-17-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
% Protocolo(43).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v831/reliability/v831-12-12-17-sitio1-E3-reliability-full-movie-bandwidth-15.mat');
% Protocolo(44).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-18/sitio1/E2/v8/reliability/v8-12-12-18-sitio1-E2-reliability-full-movie-bandwidth-15.mat');
% Protocolo(45).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-18/sitio1/E2/v11/reliability/v11-12-12-18-sitio1-E2-reliability-full-movie-bandwidth-15.mat');
% Protocolo(46).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-19/sitio1/E2/v1/reliability/v1-12-12-19-sitio1-E2-reliability-full-movie-bandwidth-15.mat');
% Protocolo(47).real = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-19/sitio1/E2/v8/reliability/v8-12-12-19-sitio1-E2-reliability-full-movie-bandwidth-15.mat');

% Protocolo(1).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v1/sparseness/v1-12-10-22-sitio1-E2-bin_size-30-sparseness.mat');
% Protocolo(2).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v3/sparseness/v3-12-10-22-sitio1-E2-bin_size-30-sparseness.mat');
% Protocolo(3).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v2/sparseness/v2-12-10-23-sitio1-E2-bin_size-30-sparseness.mat');
% Protocolo(4).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v4/sparseness/v4-12-10-23-sitio1-E2-bin_size-30-sparseness.mat');
% Protocolo(5).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E1/v3/sparseness/v3-12-10-24-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(6).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E2/v2/sparseness/v2-12-10-24-sitio1-E2-bin_size-30-sparseness.mat');
% Protocolo(7).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-26/sitio2/E1/v3/sparseness/v3-12-10-26-sitio2-E1-bin_size-30-sparseness.mat');
% % Protocolo(8).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v5/sparseness/v5-12-10-31-sitio1-E3-bin_size-30-sparseness.mat');
% % Protocolo(9).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v6/sparseness/v6-12-10-31-sitio1-E3-bin_size-30-sparseness.mat');
% Protocolo(8).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v5/sparseness/v5-12-11-01-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(9).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v6/sparseness/v6-12-11-01-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(10).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v5/sparseness/v5-12-11-06-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(11).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v6/sparseness/v6-12-11-06-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(12).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v5/sparseness/v5-12-11-07-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(13).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v6/sparseness/v6-12-11-07-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(14).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v1/sparseness/v1-12-11-08-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(15).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v4/sparseness/v4-12-11-08-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(16).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v3/sparseness/v3-12-11-09-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(17).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v6/sparseness/v6-12-11-09-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(18).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v1/sparseness/v1-12-11-13-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(19).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v2/sparseness/v2-12-11-13-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(20).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v8/sparseness/v8-12-12-03-sitio1-E3-bin_size-30-sparseness.mat');
% Protocolo(21).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v11/sparseness/v11-12-12-03-sitio1-E3-bin_size-30-sparseness.mat');
% Protocolo(22).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v9/sparseness/v9-12-12-04-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(23).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v10/sparseness/v10-12-12-04-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(24).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-05/sitio1/E1/v6/sparseness/v6-12-12-05-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(25).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-05/sitio1/E1/v8/sparseness/v8-12-12-05-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(26).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-07/sitio1/E1/v3/sparseness/v3-12-12-07-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(27).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-07/sitio1/E1/v8/sparseness/v8-12-12-07-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(28).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-10/sitio1/E2/v5/sparseness/v5-12-12-10-sitio1-E2-bin_size-30-sparseness.mat');
% Protocolo(29).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-10/sitio1/E2/v8/sparseness/v8-12-12-10-sitio1-E2-bin_size-30-sparseness.mat');
% Protocolo(30).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-12/sitio2/E1/v3/sparseness/v3-12-12-12-sitio2-E1-bin_size-30-sparseness.mat');
% Protocolo(31).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-12/sitio2/E1/v4/sparseness/v4-12-12-12-sitio2-E1-bin_size-30-sparseness.mat');
% Protocolo(32).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-13/sitio1/E1/v8/sparseness/v8-12-12-13-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(33).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-13/sitio1/E1/v2/sparseness/v2-12-12-13-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(34).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v511/sparseness/v511-12-12-17-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(35).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v512/sparseness/v512-12-12-17-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(36).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v521/sparseness/v521-12-12-17-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(37).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v811/sparseness/v811-12-12-17-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(38).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v812/sparseness/v812-12-12-17-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(39).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v821/sparseness/v821-12-12-17-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(40).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v822/sparseness/v822-12-12-17-sitio1-E1-bin_size-30-sparseness.mat');
% Protocolo(41).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v531/sparseness/v531-12-12-17-sitio1-E3-bin_size-30-sparseness.mat');
% Protocolo(42).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v532/sparseness/v532-12-12-17-sitio1-E3-bin_size-30-sparseness.mat');
% Protocolo(43).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v831/sparseness/v831-12-12-17-sitio1-E3-bin_size-30-sparseness.mat');
% Protocolo(44).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-18/sitio1/E2/v8/sparseness/v8-12-12-18-sitio1-E2-bin_size-30-sparseness.mat');
% Protocolo(45).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-18/sitio1/E2/v11/sparseness/v11-12-12-18-sitio1-E2-bin_size-30-sparseness.mat');
% Protocolo(46).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-19/sitio1/E2/v1/sparseness/v1-12-12-19-sitio1-E2-bin_size-30-sparseness.mat');
% Protocolo(47).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-19/sitio1/E2/v8/sparseness/v8-12-12-19-sitio1-E2-bin_size-30-sparseness.mat');

% Protocolo(1).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v1/TI/v1-12-10-22-sitio1-E2-bin_size-30-TI.mat');
% Protocolo(2).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v3/TI/v3-12-10-22-sitio1-E2-bin_size-30-TI.mat');
% Protocolo(3).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v2/TI/v2-12-10-23-sitio1-E2-bin_size-30-TI.mat');
% Protocolo(4).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v4/TI/v4-12-10-23-sitio1-E2-bin_size-30-TI.mat');
% Protocolo(5).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E1/v3/TI/v3-12-10-24-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(6).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E2/v2/TI/v2-12-10-24-sitio1-E2-bin_size-30-TI.mat');
% Protocolo(7).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-26/sitio2/E1/v3/TI/v3-12-10-26-sitio2-E1-bin_size-30-TI.mat');
% % Protocolo(8).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v5/TI/v5-12-10-31-sitio1-E3-bin_size-30-TI.mat');
% % Protocolo(9).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v6/TI/v6-12-10-31-sitio1-E3-bin_size-30-TI.mat');
% Protocolo(8).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v5/TI/v5-12-11-01-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(9).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v6/TI/v6-12-11-01-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(10).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v5/TI/v5-12-11-06-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(11).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v6/TI/v6-12-11-06-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(12).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v5/TI/v5-12-11-07-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(13).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v6/TI/v6-12-11-07-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(14).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v1/TI/v1-12-11-08-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(15).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v4/TI/v4-12-11-08-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(16).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v3/TI/v3-12-11-09-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(17).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v6/TI/v6-12-11-09-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(18).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v1/TI/v1-12-11-13-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(19).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v2/TI/v2-12-11-13-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(20).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v8/TI/v8-12-12-03-sitio1-E3-bin_size-30-TI.mat');
% Protocolo(21).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v11/TI/v11-12-12-03-sitio1-E3-bin_size-30-TI.mat');
% Protocolo(22).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v9/TI/v9-12-12-04-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(23).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v10/TI/v10-12-12-04-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(24).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-05/sitio1/E1/v6/TI/v6-12-12-05-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(25).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-05/sitio1/E1/v8/TI/v8-12-12-05-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(26).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-07/sitio1/E1/v3/TI/v3-12-12-07-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(27).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-07/sitio1/E1/v8/TI/v8-12-12-07-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(28).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-10/sitio1/E2/v5/TI/v5-12-12-10-sitio1-E2-bin_size-30-TI.mat');
% Protocolo(29).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-10/sitio1/E2/v8/TI/v8-12-12-10-sitio1-E2-bin_size-30-TI.mat');
% Protocolo(30).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-12/sitio2/E1/v3/TI/v3-12-12-12-sitio2-E1-bin_size-30-TI.mat');
% Protocolo(31).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-12/sitio2/E1/v4/TI/v4-12-12-12-sitio2-E1-bin_size-30-TI.mat');
% Protocolo(32).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-13/sitio1/E1/v8/TI/v8-12-12-13-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(33).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-13/sitio1/E1/v2/TI/v2-12-12-13-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(34).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v511/TI/v511-12-12-17-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(35).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v512/TI/v512-12-12-17-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(36).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v521/TI/v521-12-12-17-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(37).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v811/TI/v811-12-12-17-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(38).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v812/TI/v812-12-12-17-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(39).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v821/TI/v821-12-12-17-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(40).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v822/TI/v822-12-12-17-sitio1-E1-bin_size-30-TI.mat');
% Protocolo(41).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v531/TI/v531-12-12-17-sitio1-E3-bin_size-30-TI.mat');
% Protocolo(42).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v532/TI/v532-12-12-17-sitio1-E3-bin_size-30-TI.mat');
% Protocolo(43).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v831/TI/v831-12-12-17-sitio1-E3-bin_size-30-TI.mat');
% Protocolo(44).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-18/sitio1/E2/v8/TI/v8-12-12-18-sitio1-E2-bin_size-30-TI.mat');
% Protocolo(45).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-18/sitio1/E2/v11/TI/v11-12-12-18-sitio1-E2-bin_size-30-TI.mat');
% Protocolo(46).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-19/sitio1/E2/v1/TI/v1-12-12-19-sitio1-E2-bin_size-30-TI.mat');
% Protocolo(47).TI = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-19/sitio1/E2/v8/TI/v8-12-12-19-sitio1-E2-bin_size-30-TI.mat');

Protocolo(1).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v1/suppression/v1-12-10-22-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(2).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v3/suppression/v3-12-10-22-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(3).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v2/suppression/v2-12-10-23-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(4).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v4/suppression/v4-12-10-23-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(5).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E1/v3/suppression/v3-12-10-24-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(6).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E2/v2/suppression/v2-12-10-24-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(7).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-26/sitio2/E1/v3/suppression/v3-12-10-26-sitio2-E1-bin_size-30-suppression.mat');
% Protocolo(8).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v5/suppression/v5-12-10-31-sitio1-E3-bin_size-30-suppression.mat');
% Protocolo(9).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v6/suppression/v6-12-10-31-sitio1-E3-bin_size-30-suppression.mat');
Protocolo(8).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v5/suppression/v5-12-11-01-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(9).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v6/suppression/v6-12-11-01-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(10).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v5/suppression/v5-12-11-06-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(11).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v6/suppression/v6-12-11-06-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(12).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v5/suppression/v5-12-11-07-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(13).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v6/suppression/v6-12-11-07-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(14).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v1/suppression/v1-12-11-08-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(15).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v4/suppression/v4-12-11-08-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(16).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v3/suppression/v3-12-11-09-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(17).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v6/suppression/v6-12-11-09-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(18).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v1/suppression/v1-12-11-13-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(19).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v2/suppression/v2-12-11-13-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(20).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v8/suppression/v8-12-12-03-sitio1-E3-bin_size-30-suppression.mat');
Protocolo(21).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v11/suppression/v11-12-12-03-sitio1-E3-bin_size-30-suppression.mat');
Protocolo(22).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v9/suppression/v9-12-12-04-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(23).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v10/suppression/v10-12-12-04-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(24).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-05/sitio1/E1/v6/suppression/v6-12-12-05-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(25).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-05/sitio1/E1/v8/suppression/v8-12-12-05-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(26).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-07/sitio1/E1/v3/suppression/v3-12-12-07-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(27).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-07/sitio1/E1/v8/suppression/v8-12-12-07-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(28).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-10/sitio1/E2/v5/suppression/v5-12-12-10-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(29).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-10/sitio1/E2/v8/suppression/v8-12-12-10-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(30).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-12/sitio2/E1/v3/suppression/v3-12-12-12-sitio2-E1-bin_size-30-suppression.mat');
Protocolo(31).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-12/sitio2/E1/v4/suppression/v4-12-12-12-sitio2-E1-bin_size-30-suppression.mat');
Protocolo(32).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-13/sitio1/E1/v8/suppression/v8-12-12-13-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(33).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-13/sitio1/E1/v2/suppression/v2-12-12-13-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(34).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v511/suppression/v511-12-12-17-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(35).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v512/suppression/v512-12-12-17-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(36).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v521/suppression/v521-12-12-17-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(37).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v811/suppression/v811-12-12-17-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(38).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v812/suppression/v812-12-12-17-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(39).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v821/suppression/v821-12-12-17-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(40).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v822/suppression/v822-12-12-17-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(41).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v531/suppression/v531-12-12-17-sitio1-E3-bin_size-30-suppression.mat');
Protocolo(42).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v532/suppression/v532-12-12-17-sitio1-E3-bin_size-30-suppression.mat');
Protocolo(43).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v831/suppression/v831-12-12-17-sitio1-E3-bin_size-30-suppression.mat');
Protocolo(44).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-18/sitio1/E2/v8/suppression/v8-12-12-18-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(45).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-18/sitio1/E2/v11/suppression/v11-12-12-18-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(46).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-19/sitio1/E2/v1/suppression/v1-12-12-19-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(47).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-19/sitio1/E2/v8/suppression/v8-12-12-19-sitio1-E2-bin_size-30-suppression.mat');

% Protocolo(1).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v1/full-movie/full-movie-original/v1-WOSP-12-10-22-sitio1-E2-full-movie-original.mat');
% Protocolo(2).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v3/full-movie/full-movie-original/v3-WOSP-12-10-22-sitio1-E2-full-movie-original.mat');
% Protocolo(3).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v2/full-movie/full-movie-original/v2-WOSP-12-10-23-sitio1-E2-full-movie-original.mat');
% Protocolo(4).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v4/full-movie/full-movie-original/v4-WOSP-12-10-23-sitio1-E2-full-movie-original.mat');
% Protocolo(5).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E1/v3/full-movie/full-movie-original/v3-WOSP-12-10-24-sitio1-E1-full-movie-original.mat');
% Protocolo(6).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E2/v2/full-movie/full-movie-original/v2-WOSP-12-10-24-sitio1-E2-full-movie-original.mat');
% Protocolo(7).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-26/sitio2/E1/v3/full-movie/full-movie-original/v3-WOSP-12-10-26-sitio2-E1-full-movie-original.mat');
% % Protocolo(8).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v5/full-movie/full-movie-original/v5-WOSP-12-10-31-sitio1-E3-full-movie-original.mat');
% % Protocolo(9).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v6/full-movie/full-movie-original/v6-WOSP-12-10-31-sitio1-E3-full-movie-original.mat');
% Protocolo(8).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v5/full-movie/full-movie-original/v5-WOSP-12-11-01-sitio1-E1-full-movie-original.mat');
% Protocolo(9).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v6/full-movie/full-movie-original/v6-WOSP-12-11-01-sitio1-E1-full-movie-original.mat');
% Protocolo(10).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v5/full-movie/full-movie-original/v5-WOSP-12-11-06-sitio1-E1-full-movie-original.mat');
% Protocolo(11).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v6/full-movie/full-movie-original/v6-WOSP-12-11-06-sitio1-E1-full-movie-original.mat');
% Protocolo(12).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v5/full-movie/full-movie-original/v5-WOSP-12-11-07-sitio1-E1-full-movie-original.mat');
% Protocolo(13).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v6/full-movie/full-movie-original/v6-WOSP-12-11-07-sitio1-E1-full-movie-original.mat');
% Protocolo(14).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v1/full-movie/full-movie-original/v1-WOSP-12-11-08-sitio1-E1-full-movie-original.mat');
% Protocolo(15).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v4/full-movie/full-movie-original/v4-WOSP-12-11-08-sitio1-E1-full-movie-original.mat');
% Protocolo(16).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v3/full-movie/full-movie-original/v3-WOSP-12-11-09-sitio1-E1-full-movie-original.mat');
% Protocolo(17).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v6/full-movie/full-movie-original/v6-WOSP-12-11-09-sitio1-E1-full-movie-original.mat');
% Protocolo(18).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v1/full-movie/full-movie-original/v1-WOSP-12-11-13-sitio1-E1-full-movie-original.mat');
% Protocolo(19).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v2/full-movie/full-movie-original/v2-WOSP-12-11-13-sitio1-E1-full-movie-original.mat');
% Protocolo(20).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v8/full-movie/full-movie-original/v8-WOSP-12-12-03-sitio1-E3-full-movie-original.mat');
% Protocolo(21).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v11/full-movie/full-movie-original/v11-WOSP-12-12-03-sitio1-E3-full-movie-original.mat');
% Protocolo(22).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v9/full-movie/full-movie-original/v9-WOSP-12-12-04-sitio1-E1-full-movie-original.mat');
% Protocolo(23).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v10/full-movie/full-movie-original/v10-WOSP-12-12-04-sitio1-E1-full-movie-original.mat');
% Protocolo(24).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-05/sitio1/E1/v6/full-movie/full-movie-original/v6-WOSP-12-12-05-sitio1-E1-full-movie-original.mat');
% Protocolo(25).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-05/sitio1/E1/v8/full-movie/full-movie-original/v8-WOSP-12-12-05-sitio1-E1-full-movie-original.mat');
% Protocolo(26).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-07/sitio1/E1/v3/full-movie/full-movie-original/v3-WOSP-12-12-07-sitio1-E1-full-movie-original.mat');
% Protocolo(27).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-07/sitio1/E1/v8/full-movie/full-movie-original/v8-WOSP-12-12-07-sitio1-E1-full-movie-original.mat');
% Protocolo(28).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-10/sitio1/E2/v5/full-movie/full-movie-original/v5-WOSP-12-12-10-sitio1-E2-full-movie-original.mat');
% Protocolo(29).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-10/sitio1/E2/v8/full-movie/full-movie-original/v8-WOSP-12-12-10-sitio1-E2-full-movie-original.mat');
% Protocolo(30).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-12/sitio2/E1/v3/full-movie/full-movie-original/v3-WOSP-12-12-12-sitio2-E1-full-movie-original.mat');
% Protocolo(31).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-12/sitio2/E1/v4/full-movie/full-movie-original/v4-WOSP-12-12-12-sitio2-E1-full-movie-original.mat');
% Protocolo(32).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-13/sitio1/E1/v8/full-movie/full-movie-original/v8-WOSP-12-12-13-sitio1-E1-full-movie-original.mat');
% Protocolo(33).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-13/sitio1/E1/v2/full-movie/full-movie-original/v2-WOSP-12-12-13-sitio1-E1-full-movie-original.mat');
% Protocolo(34).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v511/full-movie/full-movie-original/v511-WOSP-12-12-17-sitio1-E1-full-movie-original.mat');
% Protocolo(35).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v512/full-movie/full-movie-original/v512-WOSP-12-12-17-sitio1-E1-full-movie-original.mat');
% Protocolo(36).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v521/full-movie/full-movie-original/v521-WOSP-12-12-17-sitio1-E1-full-movie-original.mat');
% Protocolo(37).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v811/full-movie/full-movie-original/v811-WOSP-12-12-17-sitio1-E1-full-movie-original.mat');
% Protocolo(38).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v812/full-movie/full-movie-original/v812-WOSP-12-12-17-sitio1-E1-full-movie-original.mat');
% Protocolo(39).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v821/full-movie/full-movie-original/v821-WOSP-12-12-17-sitio1-E1-full-movie-original.mat');
% Protocolo(40).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v822/full-movie/full-movie-original/v822-WOSP-12-12-17-sitio1-E1-full-movie-original.mat');
% Protocolo(41).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v531/full-movie/full-movie-original/v531-WOSP-12-12-17-sitio1-E3-full-movie-original.mat');
% Protocolo(42).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v532/full-movie/full-movie-original/v532-WOSP-12-12-17-sitio1-E3-full-movie-original.mat');
% Protocolo(43).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v831/full-movie/full-movie-original/v831-WOSP-12-12-17-sitio1-E3-full-movie-original.mat');
% Protocolo(44).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-18/sitio1/E2/v8/full-movie/full-movie-original/v8-WOSP-12-12-18-sitio1-E2-full-movie-original.mat');
% Protocolo(45).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-18/sitio1/E2/v11/full-movie/full-movie-original/v11-WOSP-12-12-18-sitio1-E2-full-movie-original.mat');
% Protocolo(46).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-19/sitio1/E2/v1/full-movie/full-movie-original/v1-WOSP-12-12-19-sitio1-E2-full-movie-original.mat');
% Protocolo(47).metric = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-19/sitio1/E2/v8/full-movie/full-movie-original/v8-WOSP-12-12-19-sitio1-E2-full-movie-original.mat');

%%% PROTOCOLOS FILMES DIRETO
% 
% Protocolo(1).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-20/sitio1/E1/v3/full-movie/full-movie-original/v3-WOSP-20-08-12-sitio1-E1-full-movie-original.mat');
% Protocolo(2).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-27/sitio1/E2/v3/full-movie/full-movie-original/v3-WOSP-27-08-12-sitio1-E2-full-movie-original.mat');
% Protocolo(3).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio1/E1/v3/full-movie/full-movie-original/v3-WOSP-29-08-12-sitio1-E1-full-movie-original.mat');
% Protocolo(4).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio1/E1/v4/full-movie/full-movie-original/v4-WOSP-29-08-12-sitio1-E1-full-movie-original.mat');
% Protocolo(5).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio2/E1/v1/full-movie/full-movie-original/v1-WOSP-29-08-12-sitio2-E1-full-movie-original.mat');
% Protocolo(6).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E1/v1/full-movie/full-movie-original/v1-WOSP-30-08-12-sitio1-E1-full-movie-original.mat');
% Protocolo(7).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E1/v4/full-movie/full-movie-original/v4-WOSP-30-08-12-sitio1-E1-full-movie-original.mat');
% Protocolo(8).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v1/full-movie/full-movie-original/v1-WOSP-30-08-12-sitio1-E3-full-movie-original.mat');
% Protocolo(9).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v3/full-movie/full-movie-original/v3-WOSP-30-08-12-sitio1-E3-full-movie-original.mat');
% Protocolo(10).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v4/full-movie/full-movie-original/v4-WOSP-30-08-12-sitio1-E3-full-movie-original.mat');
% Protocolo(11).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v1/full-movie/full-movie-original/v1-WOSP-31-08-12-sitio1-E1-full-movie-original.mat');
% Protocolo(12).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v3/full-movie/full-movie-original/v3-WOSP-31-08-12-sitio1-E1-full-movie-original.mat');
% Protocolo(13).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v4/full-movie/full-movie-original/v4-WOSP-31-08-12-sitio1-E1-full-movie-original.mat');
% Protocolo(14).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v1/full-movie/full-movie-original/v1-WOSP-31-08-12-sitio1-E3-full-movie-original.mat');
% Protocolo(15).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v3/full-movie/full-movie-original/v3-WOSP-31-08-12-sitio1-E3-full-movie-original.mat');
% Protocolo(16).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v4/full-movie/full-movie-original/v4-WOSP-31-08-12-sitio1-E3-full-movie-original.mat');
% Protocolo(17).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-03/sitio1/E3/v2/full-movie/full-movie-original/v2-WOSP-12-09-03-sitio1-E3-full-movie-original.mat');
% Protocolo(18).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-03/sitio2/E1/v4/full-movie/full-movie-original/v4-WOSP-12-09-03-sitio2-E1-full-movie-original.mat');
% Protocolo(19).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-04/sitio1/E1/v2/full-movie/full-movie-original/v2-WOSP-12-09-04-sitio1-E1-full-movie-original.mat');
% Protocolo(20).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-04/sitio1/E1/v4/full-movie/full-movie-original/v4-WOSP-12-09-04-sitio1-E1-full-movie-original.mat');
% Protocolo(21).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-05/sitio1/E1/v3/full-movie/full-movie-original/v3-WOSP-12-09-05-sitio1-E1-full-movie-original.mat');
% Protocolo(22).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-05/sitio1/E1/v4/full-movie/full-movie-original/v4-WOSP-12-09-05-sitio1-E1-full-movie-original.mat');
% Protocolo(23).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-06/sitio1/E1/v3/full-movie/full-movie-original/v3-WOSP-12-09-06-sitio1-E1-full-movie-original.mat');
% Protocolo(24).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-06/sitio1/E1/v4/full-movie/full-movie-original/v4-WOSP-12-09-06-sitio1-E1-full-movie-original.mat');
% Protocolo(25).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-10-16/sitio1/E1/v1/full-movie/full-movie-original/v1-WOSP-12-10-16-sitio1-E1-full-movie-original.mat');
% Protocolo(26).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E1/v311/full-movie/full-movie-original/v311-WOSP-12-12-17-sitio1-E1-full-movie-original.mat');
% Protocolo(27).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E1/v321/full-movie/full-movie-original/v321-WOSP-12-12-17-sitio1-E1-full-movie-original.mat');
% Protocolo(28).metricDireto = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E3/v331/full-movie/full-movie-original/v331-WOSP-12-12-17-sitio1-E3-full-movie-original.mat');
% 
% % %%%% PROTOCOLOS FILMES REVERSO
% % % 
% Protocolo(1).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-20/sitio1/E1/v3/full-movie/full-movie-invertido/v3-WOSP-12-08-20-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(2).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-27/sitio1/E2/v3/full-movie/full-movie-invertido/v3-WOSP-12-08-27-sitio1-E2-full-movie-invertido-bin-size-1.mat');
% Protocolo(3).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio1/E1/v3/full-movie/full-movie-invertido/v3-WOSP-12-08-29-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(4).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio1/E1/v4/full-movie/full-movie-invertido/v4-WOSP-12-08-29-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(5).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-29/sitio2/E1/v1/full-movie/full-movie-invertido/v1-WOSP-12-08-29-sitio2-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(6).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E1/v1/full-movie/full-movie-invertido/v1-WOSP-12-08-30-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(7).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E1/v4/full-movie/full-movie-invertido/v4-WOSP-12-08-30-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(8).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v1/full-movie/full-movie-invertido/v1-WOSP-12-08-30-sitio1-E3-full-movie-invertido-bin-size-1.mat');
% Protocolo(9).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v3/full-movie/full-movie-invertido/v3-WOSP-12-08-30-sitio1-E3-full-movie-invertido-bin-size-1.mat');
% Protocolo(10).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-30/sitio1/E3/v4/full-movie/full-movie-invertido/v4-WOSP-12-08-30-sitio1-E3-full-movie-invertido-bin-size-1.mat');
% Protocolo(11).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v1/full-movie/full-movie-invertido/v1-WOSP-12-08-31-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(12).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v3/full-movie/full-movie-invertido/v3-WOSP-12-08-31-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(13).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E1/v4/full-movie/full-movie-invertido/v4-WOSP-12-08-31-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(14).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v1/full-movie/full-movie-invertido/v1-WOSP-12-08-31-sitio1-E3-full-movie-invertido-bin-size-1.mat');
% Protocolo(15).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v3/full-movie/full-movie-invertido/v3-WOSP-12-08-31-sitio1-E3-full-movie-invertido-bin-size-1.mat');
% Protocolo(16).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-08-31/sitio1/E3/v4/full-movie/full-movie-invertido/v4-WOSP-12-08-31-sitio1-E3-full-movie-invertido-bin-size-1.mat');
% Protocolo(17).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-03/sitio1/E3/v2/full-movie/full-movie-invertido/v2-WOSP-12-09-03-sitio1-E3-full-movie-invertido-bin-size-1.mat');
% Protocolo(18).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-03/sitio2/E1/v4/full-movie/full-movie-invertido/v4-WOSP-12-09-03-sitio2-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(19).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-04/sitio1/E1/v2/full-movie/full-movie-invertido/v2-WOSP-12-09-04-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(20).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-04/sitio1/E1/v4/full-movie/full-movie-invertido/v4-WOSP-12-09-04-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(21).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-05/sitio1/E1/v3/full-movie/full-movie-invertido/v3-WOSP-12-09-05-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(22).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-05/sitio1/E1/v4/full-movie/full-movie-invertido/v4-WOSP-12-09-05-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(23).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-06/sitio1/E1/v3/full-movie/full-movie-invertido/v3-WOSP-12-09-06-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(24).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-09-06/sitio1/E1/v4/full-movie/full-movie-invertido/v4-WOSP-12-09-06-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(25).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-10-16/sitio1/E1/v1/full-movie/full-movie-invertido/v1-WOSP-12-10-16-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(26).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E1/v311/full-movie/full-movie-invertido/v311-WOSP-12-12-17-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(27).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E1/v321/full-movie/full-movie-invertido/v321-WOSP-12-12-17-sitio1-E1-full-movie-invertido-bin-size-1.mat');
% Protocolo(28).metricReverso = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Analise/12-12-17/sitio1/E3/v331/full-movie/full-movie-invertido/v331-WOSP-12-12-17-sitio1-E3-full-movie-invertido-bin-size-1.mat');


% for i=1:size(Protocol,2)
%     
%     
%    diretoDireto(i) = Protocol(i).reliability.reliabilityFullMovieData.DiretoDireto.schreiberData.reliability; 
%    
%    reversoReverso(i) = Protocol(i).reliability.reliabilityFullMovieData.ReversoReverso.schreiberData.reliability; 
%    
%    diretoReverso(i) = Protocol(i).reliability.reliabilityFullMovieData.DiretoReverso.schreiberData.reliability;
%    
%    diretoInvertido(i) = Protocol(i).reliability.reliabilityFullMovieData.DiretoInvertido.schreiberData.reliability;    
%     
%     
% end
% 
% figure;
% 
% for i=1:size(Protocolo,2)
%     
%     
%    QmaxIdxD(i) = Protocolo(i).metricDireto.metric_analysis.max_info.max_info_tpmc_idx_T;
%    
%    HTmaxD(i) = Protocolo(i).metricDireto.metric_analysis.max_info.max_info_tpmc_T;
%    
%    HImaxD(i) = Protocolo(i).metricDireto.metric_analysis.max_info.max_info_tpmc_I;
%     
%    HcountD(i) = Protocolo(i).metricDireto.metric_analysis.max_info.Hcount_tpmc_info_T; 
%     
%    HBias_TD(i) = Protocolo(i).metricDireto.metric_analysis.max_info.HBias_T(QmaxIdxD(i)); 
%     
%    QmaxIdxR(i) = Protocolo(i).metricReverso.metric_analysis.max_info.max_info_tpmc_idx_T;
%    
%    HTmaxR(i) = Protocolo(i).metricReverso.metric_analysis.max_info.max_info_tpmc_T;
%    
%    HImaxR(i) = Protocolo(i).metricReverso.metric_analysis.max_info.max_info_tpmc_I;
%     
%    HcountR(i) = Protocolo(i).metricReverso.metric_analysis.max_info.Hcount_tpmc_info_T; 
%     
%    HBias_TR(i) = Protocolo(i).metricReverso.metric_analysis.max_info.HBias_T(QmaxIdxR(i)); 
%    
% end

% diretoDireto(isnan(diretoDireto)) = 0;
% reversoReverso(isnan(reversoReverso)) = 0;
% diretoReverso(isnan(diretoReverso)) = 0;
% diretoInvertido(isnan(diretoInvertido)) = 0;

% plot(HTmaxD,HTmaxR,'r.','MarkerSize',15);
% %s = fitoptions('Method','NonLinearLeastSquares','Lower',[0,0],'Upper',[Inf,Inf],'StartPoint',[0,0,0]);
% g = fittype('a*x + b','coeff',{'a','b'});
% x = HTmaxD.';
% [X, gofReverso] = fit(x,HTmaxR.',g);
% hold on;
% plot(X);

% 
% [r,p] = pearson(HTmaxD,diretoDireto)
% HT = HTmaxD';
% direto = diretoDireto';
% [RHO,PVAL] = corr(HT,direto,'type','Pearson')
% 
% figure;
% plot(HTmaxD,diretoDireto,'ro');
% xlabel('Hmax[spike]','FontSize',30);
% ylabel('Fidedignidade','FontSize',30);
% title('{Filme Direto x Entropia M\''axima}','interpreter','latex','FontSize',30);
% %s = fitoptions('Method','NonLinearLeastSquares','Lower',[0,0],'Upper',[Inf,Inf],'StartPoint',[0,0,0]);
% g = fittype('a*x + b','coeff',{'a','b'});
% x = (1/28:1/28:1).';
% [X, gofDireto] = fit(x,diretoDireto.',g);
% hold on;
% plot(X);
% xlabel('Hmax[spike]','FontSize',30);
% ylabel('Fidedignidade','FontSize',30);
% 
% [r,p] = pearson(HTmaxD,diretoReverso)
% reverso = diretoReverso';
% [RHO,PVAL] = corr(HT,reverso,'type','Pearson')
% 
% figure
% plot(HTmaxD,diretoReverso,'ro');
% xlabel('Hmax[spike]','FontSize',30);
% ylabel('Fidedignidade','FontSize',30);
% title('{Filme DiretoXReverso x Entropia M\''axima}','interpreter','latex','FontSize',30);
% %s = fitoptions('Method','NonLinearLeastSquares','Lower',[0,0],'Upper',[Inf,Inf],'StartPoint',[0,0,0]);
% g = fittype('a*x + b','coeff',{'a','b'});
% x = (1/28:1/28:1).';
% [X, gofReverso] = fit(x,diretoReverso.',g);
% hold on;
% plot(X);
% xlabel('Hmax[spike]','FontSize',30);
% ylabel('Fidedignidade','FontSize',30);
% 
% [r,p] = pearson(HTmaxD,diretoInvertido)
% invertido = diretoInvertido';
% [RHO,PVAL] = corr(HT,invertido,'type','Pearson')
% 
% figure;
% plot(HTmaxD,diretoInvertido,'ro');
% xlabel('Hmax[spike]','FontSize',30);
% ylabel('Fidedignidade','FontSize',30);
% title('{Filme DiretoXInvertido x Entropia M\''axima}','interpreter', 'latex','FontSize',30);
% %s = fitoptions('Method','NonLinearLeastSquares','Lower',[0,0],'Upper',[Inf,Inf],'StartPoint',[0,0,0]);
% g = fittype('a*x + b','coeff',{'a','b'});
% x = (1/28:1/28:1).';
% [X, gofInvertido] = fit(x,diretoInvertido.',g);
% hold on;
% plot(X);
% xlabel('Hmax[spike]','FontSize',30);
% ylabel('Fidedignidade','FontSize',30);

% Qmax = [0 0.0625 0.125 0.25 0.5 1 2 4 8 16 32 64 128 256 512];
% 
% k = 1;
% j = 1;
% for i=1:size(QmaxIdxD,2)
%     
%     ResolutionD(i) = 1./Qmax(QmaxIdxD(i));
%     
%     
%      if 1./Qmax(QmaxIdxD(i)) >= 1;
%          
%          ResolutionM1(j) =  1./Qmax(QmaxIdxD(i));
%          j = j + 1;
%      else
%          
%          Resolutionm1(k) =  1./Qmax(QmaxIdxD(i));
%          k = k + 1;
%      end
% 
% end

% [unique, nunique] = count_unique(ResolutionM1);
% 
% for i=1:size(unique,1)
%     
%    
%     histResolutionM1(unique(i)) = nunique(i);
%     
% end
% 
% figure
% bar(histResolutionM1);
% xlabel('{Precis\~ao Temporal (s)}','interpreter','latex','FontSize',30);
% ylabel('{Frequ\^encia}','interpreter','latex','FontSize',30);
% title('{Distribui\c{c}\~ao do 1/q(max)}','interpreter','latex','FontSize',30);
% 
% [unique, nunique] = count_unique(ResolutionD);
% 
% figure
% bar(unique,nunique);
% xlabel('{Precis\~ao Temporal (s)}','interpreter','latex','FontSize',30);
% ylabel('{Frequ\^encia}','interpreter','latex','FontSize',30);
% title('{Distribui\c{c}\~ao do 1/q(max) - filme direto}','interpreter','latex','FontSize',30);
% xlim([0 2]);
% 
% 
% figure
% plot(HcountD,HTmaxD,'ro');
% xlabel('H[count]','FontSize',30);
% ylabel('H[spike]','FontSize',30);
% title('{Rela\c{c}\~ao entre H[spike] e H[count] - filme direto}','interpreter','latex','FontSize',30);
% xlim([0 1]);
% hold on;
% plot(0:0.01:1,0:0.01:1);
% 
% figure
% plot(HImaxD,HTmaxD,'ro');
% xlabel('H[interval]','FontSize',30);
% ylabel('H[spike]','FontSize',30);
% title('{Rela\c{c}\~ao entre H[spike] e H[interval] - filme direto}','interpreter','latex','FontSize',30);
% xlim([0 1]);
% hold on;
% plot(0:0.01:1,0:0.01:1);


% 
% k = 1;
% j = 1;
% for i=1:size(QmaxIdxR,2)
%     
%     ResolutionR(i) = 1./Qmax(QmaxIdxR(i));
%     
%      if 1./Qmax(QmaxIdxR(i)) >= 1;
%          
%          ResolutionM1(j) =  1./Qmax(QmaxIdxR(i));
%          j = j + 1;
%      else
%          
%          Resolutionm1(k) =  1./Qmax(QmaxIdxR(i));
%          k = k + 1;
%      end
% 
% end

% [unique, nunique] = count_unique(ResolutionM1);
% 
% for i=1:size(unique,1)
%     
%    
%     histResolutionM1(unique(i)) = nunique(i);
%     
% end
% 
% figure
% bar(histResolutionM1);
% xlabel('{Precis\~ao Temporal (s)}','interpreter','latex','FontSize',30);
% ylabel('{Frequ\^encia}','interpreter','latex','FontSize',30);
% title('{Distribui\c{c}\~ao do 1/q(max)}','interpreter','latex','FontSize',30);
% 
% [unique, nunique] = count_unique(ResolutionR);
% 
% figure
% bar(unique,nunique);
% xlabel('{Precis\~ao Temporal (s)}','interpreter','latex','FontSize',30);
% ylabel('{Frequ\^encia}','interpreter','latex','FontSize',30);
% title('{Distribui\c{c}\~ao do 1/q(max) - filme reverso}','interpreter','latex','FontSize',30);
% xlim([0 2]);
% 
% 
% figure
% plot(HcountR,HTmaxR,'ro');
% xlabel('H[count]','FontSize',30);
% ylabel('H[spike]','FontSize',30);
% title('{Rela\c{c}\~ao entre H[spike] e H[count] - filme reverso}','interpreter','latex','FontSize',30);
% xlim([0 1]);
% hold on;
% plot(0:0.01:1,0:0.01:1);
% 
% figure
% plot(HImaxR,HTmaxR,'ro');
% xlabel('H[interval]','FontSize',30);
% ylabel('H[spike]','FontSize',30);
% title('{Rela\c{c}\~ao entre H[spike] e H[interval] - filme reverso}','interpreter','latex','FontSize',30);
% xlim([0 1]);
% hold on;
% plot(0:0.01:1,0:0.01:1);

% hHTmaxD = lillietest(HTmaxD)
% hHTmaxR = lillietest(HTmaxR)
% 
%  hResolutionD = lillietest(ResolutionD)
%  hResolutionR = lillietest(ResolutionR)
% 
% [p,h] = signrank(HTmaxD,HTmaxR)
% 
%  [p,h] = signrank(ResolutionD,ResolutionR)
% 
% median(HTmaxD)
% median(HTmaxR)
% 
%  median(ResolutionD)
%  median(ResolutionR)

% figure;
% 
% ResolutionAll = [ResolutionD ResolutionR];
% 
% [unique, nunique] = count_unique(ResolutionAll);
% 
% % [uniqueAgain, nuniqueAgain] = count_unique(nunique);
% % 
% % indexAll = [];
% % 
% % for i=1:length(uniqueAgain)
% %    
% %     idxs = [];
% %     
% %     if uniqueAgain(i) ~= 1
% %        
% %         idxs = find(nunique == uniqueAgain(i));
% %         
% %         vals = unique(idxs);
% %         
% %         maxVal = max(vals);
% %         
% %         index = find(unique == maxVal);
% %         
% %         idxs(idxs == index) = [];
% %         
% %     end
% %     
% %     indexAll = [indexAll idxs];
% %     
% % end
% % 
% % unique(indexAll) = [];
% % 
% % nunique(indexAll) = [];
% 
% 
% resolucoes = unique;
% frequencias = nunique; 
% [resolucoes, SortIndex] = sort(resolucoes);
% frequencias = frequencias(SortIndex);
% 
% for i=2:length(frequencias)
%    
%     for k=1:i-1
%     
%         frequencias(i) = frequencias(i) + frequencias(k);
%         
%     end
%         
% end
% 
% frequencias = frequencias./max(frequencias);
%  
%  plot(resolucoes,frequencias);
%  hold on;
%  plot(resolucoes,frequencias,'ro');
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   RELIABILITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   RELIABILITY
% 
% 
% 
% for i=1:size(Protocolo,2)
%     
%    realCRF(i) = Protocolo(i).real.reliabilityFullMovieData.CRFxCRF.schreiberData.reliability;
%    realNonCRF(i) = Protocolo(i).real.reliabilityFullMovieData.NonCRFxNonCRF.schreiberData.reliability;
%    
%    realExtraNonCRF(i) = -1;
%    
%    if i>=20
%        
%        realExtraNonCRF(i) = Protocolo(i).real.reliabilityFullMovieData.ExtraNonCRFxExtraNonCRF.schreiberData.reliability;
%        
%    end
%        
% end
% 
% 
% grupo2CRF = realCRF([8 9 10 11 12 13 14 15 16 17 18 19]);
% grupo2NonCRF = realNonCRF([8 9 10 11 12 13 14 15 16 17 18 19]);
% 
% mean(grupo2CRF)
% mean(grupo2NonCRF)
% 
% median(grupo2CRF)
% median(grupo2NonCRF)
% 
% h2CRF = lillietest(grupo2CRF)
% h2NonCRF = lillietest(grupo2NonCRF)
% 
% [h,p] = ttest(grupo2CRF,grupo2NonCRF)
% 
% grupoAnaCRF = [realNonCRF(22) realNonCRF(23) realNonCRF(24) realNonCRF(25) realCRF(26) realCRF(27) realCRF(26) realCRF(27) realCRF(28) realCRF(29) realCRF(30) realCRF(31) realCRF(30) realCRF(31) realCRF(34) realCRF(35) realCRF(34) realCRF(35) realCRF(36) realCRF(36) realCRF(37) realCRF(38) realCRF(37) realCRF(38) realCRF(39) realCRF(40) realCRF(39) realCRF(40) realCRF(44) realCRF(45) realCRF(44) realCRF(45) realCRF(46) realCRF(47)]; 
% grupoAnaNonCRF = [realExtraNonCRF(22) realExtraNonCRF(23) realExtraNonCRF(24) realExtraNonCRF(25) realNonCRF(26) realNonCRF(27) realExtraNonCRF(26) realExtraNonCRF(27) realExtraNonCRF(28) realExtraNonCRF(29) realNonCRF(30) realNonCRF(31) realExtraNonCRF(30) realExtraNonCRF(31) realNonCRF(34) realNonCRF(35) realExtraNonCRF(34) realExtraNonCRF(35) realNonCRF(36) realExtraNonCRF(36) realNonCRF(37) realNonCRF(38) realExtraNonCRF(37) realExtraNonCRF(38) realNonCRF(39) realNonCRF(40) realExtraNonCRF(39) realExtraNonCRF(40) realNonCRF(44) realNonCRF(44) realExtraNonCRF(44) realExtraNonCRF(45) realExtraNonCRF(46) realExtraNonCRF(47)];                           
% 
% h2CRF = lillietest(grupoAnaCRF)
% h2NonCRF = lillietest(grupoAnaNonCRF)
% 
% [p,h] = signrank(grupoAnaCRF,grupoAnaNonCRF)
% 
% [unique2CRF nunique2CRF] = count_unique(grupo2CRF);
% [unique2NonCRF nunique2NonCRF] = count_unique(grupo2NonCRF);
% 
% [uniqueAnaCRF nuniqueAnaCRF] = count_unique(grupoAnaCRF);
% [uniqueAnaNonCRF nuniqueAnaNonCRF] = count_unique(grupoAnaNonCRF);
% 
% 
% [unique2CRF idx] = sort(unique2CRF);
% nunique2CRF = nunique2CRF(idx);
% 
% [unique2NonCRF idx] = sort(unique2NonCRF);
% nunique2NonCRF = nunique2NonCRF(idx);
% 
% [uniqueAnaCRF idx] = sort(uniqueAnaCRF);
% nuniqueAnaCRF = nuniqueAnaCRF(idx);
% 
% [uniqueAnaNonCRF idx] = sort(uniqueAnaNonCRF);
% nuniqueAnaNonCRF = nuniqueAnaNonCRF(idx);
% 
% freqcum2CRF(1) = nunique2CRF(1);
% for i=1:(length(nunique2CRF)-1)
%     
%     freqcum2CRF(i+1) = freqcum2CRF(i) + nunique2CRF(i+1);
%     
% end
% freqcum2CRF = freqcum2CRF./max(freqcum2CRF);
% 
% freqcum2NonCRF(1) = nunique2NonCRF(1);
% for i=1:(length(nunique2NonCRF)-1)
%     
%     freqcum2NonCRF(i+1) = freqcum2NonCRF(i) + nunique2NonCRF(i+1);
%     
% end
% freqcum2NonCRF = freqcum2NonCRF./max(freqcum2NonCRF);
% 
% freqcumAnaCRF(1) = nuniqueAnaCRF(1);
% for i=1:(length(nuniqueAnaCRF)-1)
%     
%     freqcumAnaCRF(i+1) = freqcumAnaCRF(i) + nuniqueAnaCRF(i+1);
%     
% end
% freqcumAnaCRF = freqcumAnaCRF./max(freqcumAnaCRF);
% 
% freqcumAnaNonCRF(1) = nuniqueAnaNonCRF(1);
% for i=1:(length(nuniqueAnaNonCRF)-1)
%     
%     freqcumAnaNonCRF(i+1) = freqcumAnaNonCRF(i) + nuniqueAnaNonCRF(i+1);
%     
% end
% freqcumAnaNonCRF = freqcumAnaNonCRF./max(freqcumAnaNonCRF);
% 
% % unique2CRF = round(unique2CRF.*100);
% % unique2NonCRF = round(unique2NonCRF.*100);
% % uniqueAnaCRF = round(uniqueAnaCRF.*100);
% % uniqueAnaNonCRF = round(uniqueAnaNonCRF.*100);
% 
% figure;
% plot(unique2CRF,freqcum2CRF,'r--','LineWidth',2);
% hold on;
% plot(unique2NonCRF,freqcum2NonCRF,'r:','LineWidth',2);
% hold on;
% plot(uniqueAnaCRF,freqcumAnaCRF,'b--','LineWidth',2);
% hold on;
% plot(uniqueAnaNonCRF,freqcumAnaNonCRF,'b:','LineWidth',2);
% title('Freq Cum 4 Grupos');
% xlim([0 1]);

% 
% 
% 
% 
% 
% 
% 
% 
% 
% realCRF = realCRF.';
% realNonCRF = realNonCRF.';
% realExtraNonCRF = realExtraNonCRF.';
% 
% matriz = horzcat(realCRF,realNonCRF,realExtraNonCRF);
% 
% g = figure;
% 
% colnamesReal = {'realCRF','realNonCRF','realExtraNonCRF'};
%  
% t = uitable(g,'Data',matriz,'ColumnName',colnamesReal,'Position',[10 10 700 750]);
% 
% stimuliSizeCRF = [2 2 2 2 1 -1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1];
% stimuliSizeNonCRF = [8 8 10 10 5 -1 9 5 5 5 5 5 5 7 7 3 3 4 4 2 2 2 2 2 2 5 5 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 2 2 2 2];
% stimuliSizeExtraNonCRF = [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 4 4 4 4 4 4 9 9 7 7 7 7 6 6 6 6 6 6 6 6 6 8 8 8 4 4 5 5];
% 
% g = 1;
% for k=1:8
%     
%     j = 1;
%     l = 1;
%     h = 1;
%     
%     for i=1:size(stimuliSizeCRF,2)
%        
%         if stimuliSizeCRF(i) == k
%             
%             allSizes(k).realCRF(j) = realCRF(i);
%             
%             j = j + 1;
%             
%             if (k~=1)
%                 
%                 demais(g) = realCRF(i);
%                 g = g + 1;
%             
%             end
%             
%         elseif stimuliSizeNonCRF(i) == k
%             
%             allSizes(k).realNonCRF(h) = realNonCRF(i);
%             
%             h = h + 1;
%            
%              if (k~=1)
%                 
%                 demais(g) = realNonCRF(i);
%                 g = g + 1;
%             
%              end
%             
%         elseif stimuliSizeExtraNonCRF(i) == k
%             
%             allSizes(k).realExtraNonCRF(l) = realExtraNonCRF(i);
%             
%             l = l + 1;
%             
%              if (k~=1)
%                 
%                 demais(g) = realExtraNonCRF(i);
%                 g = g + 1;
%             
%              end
%             
%         end         
%     
%     end
%     
% end
% 
% h = lillietest(demais)
% h = lillietest(allSizes(1).realCRF)
% 
% [p,h] = ranksum(demais,allSizes(1).realCRF)
% 
% mean(demais)
% mean(allSizes(1).realCRF)
% 
% median(demais)
% median(allSizes(1).realCRF)
% 
% figure;
% [unique,nunique] = count_unique(demais);
% bar(unique,nunique);
% xlabel('Fidedignidade','FontSize',30);
% ylabel('Frequencia','FontSize',30);
% xlim([0 1]);
% title('{Distribui\c{c}\~ao para maiores que 1x}','interpreter','latex','FontSize',30);
% 
% [unique, sortIdx] = sort(unique);
% nunique = nunique(sortIdx);
% 
% frequencias = nunique;
% for i=2:length(frequencias)
% 
%     for R=1:i-1
% 
%         frequencias(i) = frequencias(i) + frequencias(R);
% 
%     end
% 
% end
% 
% frequencias = frequencias./max(frequencias);
% 
% figure;
% plot(unique,frequencias)
% xlabel('Fidedignidade');
% ylabel('Frequencia Cumulativa');
% title('Maiores que 1x');
% 
% figure;
% [unique,nunique] = count_unique(allSizes(1).realCRF);
% bar(unique,nunique);
% xlabel('Fidedignidade','FontSize',30);
% ylabel('Frequencia','FontSize',30);
% xlim([0 1]);
% title('{Distribui\c{c}\~ao para 1x}','interpreter','latex','FontSize',30);
% 
% [unique, sortIdx] = sort(unique);
% nunique = nunique(sortIdx);
% 
% frequencias = nunique;
% for i=2:length(frequencias)
% 
%     for R=1:i-1
% 
%         frequencias(i) = frequencias(i) + frequencias(R);
% 
%     end
% 
% end
% 
% frequencias = frequencias./max(frequencias);
% 
% figure;
% plot(unique,frequencias)
% xlabel('Fidedignidade');
% ylabel('Frequencia Cumulativa');
% title('1x');

% figure;
% 
% for i=1:size(Protocolo,2)
%     
%     
%    QmaxIdx(i) = Protocolo(i).metric.metric_analysis.max_info.max_info_tpmc_idx_T;
%    
%    HTmax(i) = Protocolo(i).metric.metric_analysis.max_info.max_info_tpmc_T;
%    
%    HImax(i) = Protocolo(i).metric.metric_analysis.max_info.max_info_tpmc_I;
%     
%    Hcount(i) = Protocolo(i).metric.metric_analysis.max_info.Hcount_tpmc_info_T; 
%     
%    HBias_T(i) = Protocolo(i).metric.metric_analysis.max_info.HBias_T(QmaxIdx(i)); 
%     
%    if Hcount(i) < 0
%        
%        i
%        Hcount
%        
%    end
%    
% end
% 
% Qmax = [0 0.0625 0.125 0.25 0.5 1 2 4 8 16 32 64 128 256 512];
% 
% k = 1;
% j = 1;
% for i=1:size(QmaxIdx,2)
%     
%      ResolutionAll(i) = 1./Qmax(QmaxIdx(i));
%     
%      if 1./Qmax(QmaxIdx(i)) >= 1;
%          
%          ResolutionM1(j) =  1./Qmax(QmaxIdx(i));
%          j = j + 1;
%          
%      else
%          
%          Resolutionm1(k) =  1./Qmax(QmaxIdx(i));
%          k = k + 1;
%          
%      end
% 
% end
% 
% [unique, nunique] = count_unique(ResolutionM1);
% 
% for i=1:size(unique,1)
%     
%    
%     histResolutionM1(unique(i)) = nunique(i);
%     
% end
% 
% figure
% bar(histResolutionM1);
% xlabel('{Precis\~ao Temporal (s)}','interpreter','latex','FontSize',30);
% ylabel('{Frequ\^encia}','interpreter','latex','FontSize',30);
% title('{Distribui\c{c}\~ao do 1/q(max)}','interpreter','latex','FontSize',30);
% 
% [unique, nunique] = count_unique(Resolutionm1);
% 
% figure
% bar(unique,nunique);
% xlabel('{Precis\~ao Temporal (s)}','interpreter','latex','FontSize',30);
% ylabel('{Frequ\^encia}','interpreter','latex','FontSize',30);
% title('{Distribui\c{c}\~ao do 1/q(max)}','interpreter','latex','FontSize',30);
% 
% [unique, nunique] = count_unique(ResolutionAll);
% 
% [unique, sortIdx] = sort(unique);
% nunique = nunique(sortIdx);
% 
%         frequencias = nunique;
%         for i=2:length(frequencias)
% 
%             for R=1:i-1
% 
%                 frequencias(i) = frequencias(i) + frequencias(R);
% 
%             end
% 
%         end
% 
%         frequencias = frequencias./max(frequencias);
% 
% unique
% figure;
% plot(unique,frequencias);
% ylabel('Frequencia Cumulativa');
% xlabel('Resolucao');
% hold on;
% plot(unique,frequencias,'r.','markersize',10);
% 
% figure
% plot(Hcount,HTmax,'ro');
% xlabel('H[count]','FontSize',30);
% ylabel('H[spike]','FontSize',30);
% xlim([0 1]);
% ylim([0 1]);
% title('{Rela\c{c}\~ao entre H[spike] e H[count]}','interpreter','latex','FontSize',30);
% hold on;
% plot(0:0.01:1,0:0.01:1,'k');
% 
% figure
% plot(HImax,HTmax,'ro');
% xlabel('H[interval]','FontSize',30);
% ylabel('H[spike]','FontSize',30);
% xlim([0 1]);
% ylim([0 1]);
% title('{Rela\c{c}\~ao entre H[spike] e H[interval]}','interpreter','latex','FontSize',30);
% hold on;
% plot(0:0.01:1,0:0.01:1,'k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     SPARSENESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     SPARSENESS
% g = figure;
% maior = 0;
% menor = 0;
% 
% for i=1:size(Protocolo,2)
%     
%     sparsenessNonCRF(i) = Protocolo(i).sparseness.sparseData(1,1).gallant2002SelectivityIndex;
%     sparsenessCRF(i) = Protocolo(i).sparseness.sparseData(1,2).gallant2002SelectivityIndex;
%     
%     if i>=10
%         
%         sparsenessCRF(i) = Protocolo(i).sparseness.sparseData(1,1).gallant2002SelectivityIndex;
%         sparsenessNonCRF(i) = Protocolo(i).sparseness.sparseData(1,2).gallant2002SelectivityIndex;
%         
%     end
%     
%     if i>=20
%        
%         sparsenessExtraNonCRF(i) = Protocolo(i).sparseness.sparseData(1,3).gallant2002SelectivityIndex;
%         
%     end
%     
%     if sparsenessNonCRF(i) > sparsenessCRF(i)
%         
%         maior = maior + 1;
%         
%     else
%         
%         menor = menor + 1;
%         
%     end
%     
%     if i>=20
%         
%             if sparsenessExtraNonCRF(i) > sparsenessCRF(i)
% 
%                 maior = maior + 1;
% 
%             else
% 
%                 menor = menor + 1;
% 
%             end
%             
%     end
%     
% end
% 
% stimuliSizeCRF = [2 2 2 2 1 -1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1];
% stimuliSizeNonCRF = [8 8 10 10 5 -1 9 5 5 5 5 5 5 7 7 3 3 4 4 2 2 2 2 2 2 5 5 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 2 2 2 2];
% stimuliSizeExtraNonCRF = [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 4 4 4 4 4 4 9 9 7 7 7 7 6 6 6 6 6 6 6 6 6 8 8 8 4 4 5 5];
% 
% sparsenessCRF = sparsenessCRF.';
% 
% sparsenessNonCRF = sparsenessNonCRF.';
% 
% sparsenessExtraNonCRF = sparsenessExtraNonCRF.';
% 
% stimuliSizeCRF = stimuliSizeCRF.';
% 
% stimuliSizeNonCRF = stimuliSizeNonCRF.';
% 
% stimuliSizeExtraNonCRF = stimuliSizeExtraNonCRF.';
% 
% sparsenessMatrix = horzcat(sparsenessCRF,stimuliSizeCRF,sparsenessNonCRF,stimuliSizeNonCRF,sparsenessExtraNonCRF,stimuliSizeExtraNonCRF);
% 
% colnamesSparse = {'sparsenessCRF','sparsenessNonCRF','sparsenessExtraNonCRF'};
% 
% t = uitable(g,'Data',sparsenessMatrix,'ColumnName',colnamesSparse,'Position',[10 10 700 750]);
% 
% maior
% menor
% 
% grupo2CRF = sparsenessCRF([8 9 10 11 12 13 14 15 16 17 18 19]);
% grupo2NonCRF = sparsenessNonCRF([8 9 10 11 12 13 14 15 16 17 18 19]);
% 
% mean(grupo2CRF)
% mean(grupo2NonCRF)
% 
% median(grupo2CRF)
% median(grupo2NonCRF)
% 
% h2CRF = lillietest(grupo2CRF)
% h2NonCRF = lillietest(grupo2NonCRF)
% 
% [h,p] = ttest(grupo2CRF,grupo2NonCRF)
% 
% grupoAnaCRF = [sparsenessNonCRF(22) sparsenessNonCRF(23) sparsenessNonCRF(24) sparsenessNonCRF(25) sparsenessCRF(26) sparsenessCRF(27) sparsenessCRF(26) sparsenessCRF(27) sparsenessCRF(28) sparsenessCRF(29) sparsenessCRF(30) sparsenessCRF(31) sparsenessCRF(30) sparsenessCRF(31) sparsenessCRF(34) sparsenessCRF(35) sparsenessCRF(34) sparsenessCRF(35) sparsenessCRF(36) sparsenessCRF(36) sparsenessCRF(37) sparsenessCRF(38) sparsenessCRF(37) sparsenessCRF(38) sparsenessCRF(39) sparsenessCRF(40) sparsenessCRF(39) sparsenessCRF(40) sparsenessCRF(44) sparsenessCRF(45) sparsenessCRF(44) sparsenessCRF(45) sparsenessCRF(46) sparsenessCRF(47)]; 
% grupoAnaNonCRF = [sparsenessExtraNonCRF(22) sparsenessExtraNonCRF(23) sparsenessExtraNonCRF(24) sparsenessExtraNonCRF(25) sparsenessNonCRF(26) sparsenessNonCRF(27) sparsenessExtraNonCRF(26) sparsenessExtraNonCRF(27) sparsenessExtraNonCRF(28) sparsenessExtraNonCRF(29) sparsenessNonCRF(30) sparsenessNonCRF(31) sparsenessExtraNonCRF(30) sparsenessExtraNonCRF(31) sparsenessNonCRF(34) sparsenessNonCRF(35) sparsenessExtraNonCRF(34) sparsenessExtraNonCRF(35) sparsenessNonCRF(36) sparsenessExtraNonCRF(36) sparsenessNonCRF(37) sparsenessNonCRF(38) sparsenessExtraNonCRF(37) sparsenessExtraNonCRF(38) sparsenessNonCRF(39) sparsenessNonCRF(40) sparsenessExtraNonCRF(39) sparsenessExtraNonCRF(40) sparsenessNonCRF(44) sparsenessNonCRF(44) sparsenessExtraNonCRF(44) sparsenessExtraNonCRF(45) sparsenessExtraNonCRF(46) sparsenessExtraNonCRF(47)];                           
% 
% median(grupoAnaCRF)
% median(grupoAnaNonCRF)

% [unique2CRF nunique2CRF] = count_unique(grupo2CRF);
% [unique2NonCRF nunique2NonCRF] = count_unique(grupo2NonCRF);
% 
% [uniqueAnaCRF nuniqueAnaCRF] = count_unique(grupoAnaCRF);
% [uniqueAnaNonCRF nuniqueAnaNonCRF] = count_unique(grupoAnaNonCRF);
% 
% 
% [unique2CRF idx] = sort(unique2CRF);
% nunique2CRF = nunique2CRF(idx);
% 
% [unique2NonCRF idx] = sort(unique2NonCRF);
% nunique2NonCRF = nunique2NonCRF(idx);
% 
% [uniqueAnaCRF idx] = sort(uniqueAnaCRF);
% nuniqueAnaCRF = nuniqueAnaCRF(idx);
% 
% [uniqueAnaNonCRF idx] = sort(uniqueAnaNonCRF);
% nuniqueAnaNonCRF = nuniqueAnaNonCRF(idx);
% 
% freqcum2CRF(1) = nunique2CRF(1);
% for i=1:(length(nunique2CRF)-1)
%     
%     freqcum2CRF(i+1) = freqcum2CRF(i) + nunique2CRF(i+1);
%     
% end
% freqcum2CRF = freqcum2CRF./max(freqcum2CRF);
% 
% freqcum2NonCRF(1) = nunique2NonCRF(1);
% for i=1:(length(nunique2NonCRF)-1)
%     
%     freqcum2NonCRF(i+1) = freqcum2NonCRF(i) + nunique2NonCRF(i+1);
%     
% end
% freqcum2NonCRF = freqcum2NonCRF./max(freqcum2NonCRF);
% 
% freqcumAnaCRF(1) = nuniqueAnaCRF(1);
% for i=1:(length(nuniqueAnaCRF)-1)
%     
%     freqcumAnaCRF(i+1) = freqcumAnaCRF(i) + nuniqueAnaCRF(i+1);
%     
% end
% freqcumAnaCRF = freqcumAnaCRF./max(freqcumAnaCRF);
% 
% freqcumAnaNonCRF(1) = nuniqueAnaNonCRF(1);
% for i=1:(length(nuniqueAnaNonCRF)-1)
%     
%     freqcumAnaNonCRF(i+1) = freqcumAnaNonCRF(i) + nuniqueAnaNonCRF(i+1);
%     
% end
% freqcumAnaNonCRF = freqcumAnaNonCRF./max(freqcumAnaNonCRF);
% 
% unique2CRF = round(unique2CRF.*100);
% unique2NonCRF = round(unique2NonCRF.*100);
% uniqueAnaCRF = round(uniqueAnaCRF.*100);
% uniqueAnaNonCRF = round(uniqueAnaNonCRF.*100);
% 
% figure;
% plot(unique2CRF,freqcum2CRF,'r--','LineWidth',1);
% hold on;
% plot(unique2NonCRF,freqcum2NonCRF,'r:','LineWidth',2);
% hold on;
% plot(uniqueAnaCRF,freqcumAnaCRF,'b--','LineWidth',1);
% hold on;
% plot(uniqueAnaNonCRF,freqcumAnaNonCRF,'b:','LineWidth',2);
% title('Freq Cum 4 Grupos');
% xlim([0 100]);
% 
% 
% hAnaCRF = lillietest(grupoAnaCRF)
% hAnaNonCRF = lillietest(grupoAnaNonCRF)
% 
% [h,p] = ttest(grupoAnaCRF,grupoAnaNonCRF)
% 
% size(grupoAnaCRF)

% for k=1:10
%     
%     h = 1;
%     j = 1;
%     for i=1:size(Protocolo,2)
% 
%         if stimuliSizeCRF(i) == k
% 
%             quantasVezes(k).sparseness(j) = sparsenessCRF(i);
%             
%             if k == 1
%                 
%                 um(h) = sparsenessCRF(i);
%                 umPar(h) = sparsenessNonCRF(i);
%                 h = h + 1;
%             end
%             
%           
%             j = j + 1;
%             
%             
%         elseif stimuliSizeNonCRF(i) == k
% 
%             quantasVezes(k).sparseness(j) = sparsenessNonCRF(i);
%             
%             j = j + 1;
%             
%         elseif stimuliSizeExtraNonCRF(i) == k
% 
%             quantasVezes(k).sparseness(j) = sparsenessExtraNonCRF(i);
%             j = j + 1;
%             
%         end
% 
%     end
%     
% end
% 
% 
% tam = size(quantasVezes(2).sparseness,2);
% 
% for k=1:8
%    
%     s = size(quantasVezes(k).sparseness,2);
%     s = s + 1;
%     for i=s:tam
%        
%         quantasVezes(k).sparseness(i) = NaN;
%         
%     end
%     
% end
% 
% matrix = horzcat(quantasVezes(1).sparseness.',quantasVezes(2).sparseness.',quantasVezes(3).sparseness.',quantasVezes(4).sparseness.',quantasVezes(5).sparseness.',quantasVezes(6).sparseness.',quantasVezes(7).sparseness.',quantasVezes(8).sparseness.')
% 
% [p,table,stats] = anova1(matrix)
% 
%           h = 1;
%           for i=1:size(Protocolo,2)
%                 
%                 if stimuliSizeCRF(i) == 2
%                     
%                     dois(h) = sparsenessCRF(i);
%                     doisPar(h) = sparsenessNonCRF(i);
%                     h = h + 1;
%                     
%                 end
%    
%                 if stimuliSizeNonCRF(i) == 2
%                     
%                      if i>=22
% 
%                          dois(h) = sparsenessNonCRF(i);
%                          doisPar(h) = sparsenessExtraNonCRF(i);
%                          h = h + 1;
% 
%                      end
% 
%                 end
%           end
% 
% h1 = lillietest(um)
% hPar1 = lillietest(umPar)
% 
% h2 = lillietest(dois)
% hPar2 = lillietest(doisPar)
% 
% [hum,pum] = ttest(um,umPar)
% 
% [hdois,pdois] = ttest(dois,doisPar)
% 
% mean(um)
% mean(umPar)
% 
% for k=1:10
%     
%     for i=1:size(quantasVezes(k).sparseness,2)
%         
%         quantasVezesRounded(k).sparseness(i) = round(quantasVezes(k).sparseness(i)*100);
% 
%     end
%     
% end
% 
% j = 1;
% for k=1:10
%     
%     if ((k ~= 1) && (k ~= 2) && (k ~= 4))
%         
%         for i=1:size(quantasVezesRounded(k).sparseness,2)
%             
%             demais(j) = quantasVezes(k).sparseness(i);
%             demaisR(j) = quantasVezesRounded(k).sparseness(i);
%             j = j + 1;
%             
%         end
%         
%     end
%     
% end
% 
% 
% for k=1:4
%    
%     if k~=3
%         
%         [unique(k).u numunique(k).n] = count_unique(quantasVezesRounded(k).sparseness);
% 
%         figure;
% 
%         [unique(k).u sortIdx] = sort(unique(k).u);
%         numunique(k).n = numunique(k).n(sortIdx);
% 
%         frequencias = numunique(k).n;
%         for i=2:length(frequencias)
% 
%             for R=1:i-1
% 
%                 frequencias(i) = frequencias(i) + frequencias(R);
% 
%             end
% 
%         end
% 
%         frequencias = frequencias./max(frequencias);
% 
%         plot(unique(k).u,frequencias);
%         xlabel('Sparseness (%)','FontSize',30);
%         ylabel('{Frequ\^encia}','interpreter','latex','FontSize',30);
%         title(strcat('{',int2str(k),'x o Tamanho do Campo Receptivo}'),'interpreter','latex','FontSize',30);
% 
% %         for i=1:size(unique(k).u,1)
% % 
% %             histQuantasVezes(k).count(unique(k).u(i)) = numunique(k).n(i);        
% % 
% %         end
% %         figure;
% %         bar(histQuantasVezes(k).count);
% %         title(strcat('{',int2str(k),'x o Tamanho do Campo Receptivo}'),'interpreter','latex','FontSize',30);
% %         xlabel('Sparseness (%)','FontSize',30);
% %         ylabel('{Frequ\^encia}','interpreter','latex','FontSize',30);
% %         xlim([0 100]);
%         
%     end
%     
% end
% 
%         [unique numunique] = count_unique(demaisR);
% 
%         [unique sortIdx] = sort(unique);
%         numuique = numunique(sortIdx);
% 
%         frequencias = numunique;
%         for i=2:length(frequencias)
% 
%             for R=1:i-1
% 
%                 frequencias(i) = frequencias(i) + frequencias(R);
% 
%             end
% 
%         end
% 
%         frequencias = frequencias./max(frequencias);
% 
%         figure;
% 
%         plot(unique,frequencias);
%         title(strcat('{','Demais Tamanhos de Campo Receptivo}'),'interpreter','latex','FontSize',30);
%         xlabel('Sparseness (%)','FontSize',30);
%         ylabel('{Frequ\^encia}','interpreter','latex','FontSize',30);        
% 
%         for i=1:size(unique,1)
% 
%             histDemaisR(unique(i)) = numunique(i);        
% 
%         end
%         figure;
%         bar(histDemaisR);
%         title(strcat('{','Demais Tamanhos de Campo Receptivo}'),'interpreter','latex','FontSize',30);
%         xlabel('Sparseness (%)','FontSize',30);
%         ylabel('{Frequ\^encia}','interpreter','latex','FontSize',30);
%         xlim([0 100]);
%         
%         
% mean(demais)
% std(demais)
% 
% for k=1:8
%     
%     k
%     mean(quantasVezes(k).sparseness)
%     std(quantasVezes(k).sparseness)
% end
% 
% for k=1:8
%     
%     %1 & 2
%     h1 = lillietest(quantasVezes(1).sparseness)
%     h(k) = lillietest(quantasVezes(k).sparseness)
% 
%     [h(k),p(k)] = ttest2(quantasVezes(1).sparseness,quantasVezes(k).sparseness)
% 
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   TEORIA DA INFORMACAO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   TEORIA DA INFORMACAO
% f = figure;
% 
% for i=1:size(Protocolo,2)
%    
%     condicao(2).InfoPerSecond(i) = Protocolo(i).TI.TI.condicao(1,1).InfoPerSecond;
%     condicao(1).InfoPerSecond(i) = Protocolo(i).TI.TI.condicao(1,2).InfoPerSecond;
%     
%     condicao(2).InfoPerSpike(i) = Protocolo(i).TI.TI.condicao(1,1).InfoPerSpike;
%     condicao(1).InfoPerSpike(i) = Protocolo(i).TI.TI.condicao(1,2).InfoPerSpike;
%     
%     condicao(2).CodingEfficiency(i) = Protocolo(i).TI.TI.condicao(1,1).CodingEfficiency;
%     condicao(1).CodingEfficiency(i) = Protocolo(i).TI.TI.condicao(1,2).CodingEfficiency;
%     
%     if i>=10
%        
%             condicao(1).InfoPerSecond(i) = Protocolo(i).TI.TI.condicao(1,1).InfoPerSecond;
%             condicao(2).InfoPerSecond(i) = Protocolo(i).TI.TI.condicao(1,2).InfoPerSecond;
%     
%             condicao(1).InfoPerSpike(i) = Protocolo(i).TI.TI.condicao(1,1).InfoPerSpike;
%             condicao(2).InfoPerSpike(i) = Protocolo(i).TI.TI.condicao(1,2).InfoPerSpike;
%     
%             condicao(1).CodingEfficiency(i) = Protocolo(i).TI.TI.condicao(1,1).CodingEfficiency;
%             condicao(2).CodingEfficiency(i) = Protocolo(i).TI.TI.condicao(1,2).CodingEfficiency;
%     
%     end
%     
%     if i>=20
%        
%             condicao(3).InfoPerSecond(i) = Protocolo(i).TI.TI.condicao(1,3).InfoPerSecond;
%             
%             condicao(3).InfoPerSpike(i) = Protocolo(i).TI.TI.condicao(1,3).InfoPerSpike;
%             
%             condicao(3).CodingEfficiency(i) = Protocolo(i).TI.TI.condicao(1,3).CodingEfficiency;
%             
%     end
%     
% end
% 
% InfoPerSecond1 = condicao(1).InfoPerSecond;
% InfoPerSecond2 = condicao(2).InfoPerSecond;
% InfoPerSecond3 = condicao(3).InfoPerSecond;
% 
% InfoPerSecond1 = InfoPerSecond1.';
% InfoPerSecond2 = InfoPerSecond2.';
% InfoPerSecond3 = InfoPerSecond3.';
% 
% InfoPerSpike1 = condicao(1).InfoPerSpike;
% InfoPerSpike2 = condicao(2).InfoPerSpike;
% InfoPerSpike3 = condicao(3).InfoPerSpike;
% 
% InfoPerSpike1 = InfoPerSpike1.';
% InfoPerSpike2 = InfoPerSpike2.';
% InfoPerSpike3 = InfoPerSpike3.';
% 
% CodingEfficiency1 = condicao(1).CodingEfficiency;
% CodingEfficiency2 = condicao(2).CodingEfficiency;
% CodingEfficiency3 = condicao(3).CodingEfficiency;
% 
% CodingEfficiency1 = CodingEfficiency1.';
% CodingEfficiency2 = CodingEfficiency2.';
% CodingEfficiency3 = CodingEfficiency3.';
% 
% TIMatrix1 = horzcat(InfoPerSecond1,InfoPerSecond2,InfoPerSecond3);
% TIMatrix2 = horzcat(InfoPerSpike1,InfoPerSpike2,InfoPerSpike3);
% TIMatrix3 = horzcat(CodingEfficiency1,CodingEfficiency2,CodingEfficiency3);
% 
% colnamesTI1 = {'InfoPerSecond-1','InfoPerSecond-2','InfoPerSecond-3'};
% 
% t = uitable(f,'Data',TIMatrix1,'ColumnName',colnamesTI1,'Position',[10 10 350 750]);
% 
% colnamesTI3 = {'CodEff-1','CodEff-2','CodEff-3'};
% 
% t = uitable(f,'Data',TIMatrix3,'ColumnName',colnamesTI3,'Position',[360 10 280 750]);
% 
% colnamesTI2 = {'InfoPerSpike-1','InfoPerSpike-2','InfoPerSpike-3'};
% 
% t = uitable(f,'Data',TIMatrix2,'ColumnName',colnamesTI2,'Position',[650 10 300 750]);
% 
% stimuliSizeCRF = [2 2 2 2 1 -1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1];
% stimuliSizeNonCRF = [8 8 10 10 5 -1 9 5 5 5 5 5 5 7 7 3 3 4 4 2 2 2 2 2 2 5 5 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 2 2 2 2];
% stimuliSizeExtraNonCRF = [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 4 4 4 4 4 4 9 9 7 7 7 7 6 6 6 6 6 6 6 6 6 8 8 8 4 4 5 5];
% 
% stimuliSizeCRF = stimuliSizeCRF.';
% stimuliSizeNonCRF = stimuliSizeNonCRF.';
% stimuliSizeExtraNonCRF = stimuliSizeExtraNonCRF.';
% 
% matrix = horzcat(stimuliSizeCRF,stimuliSizeNonCRF,stimuliSizeExtraNonCRF);
% colnames = {'CRF','NonCRF','ExtraNonCRF'};
% 
% %t = uitable(f,'Data',matrix,'ColumnName',colnames,'Position',[960 10 300 750]);
% 
% 
% for k=1:8
%    
%     j = 1;
%     allSizes(k).infoSecond = [];
%     allSizes(k).infoPerSpike = [];
%     allSizes(k).CodeEfficiency = [];
%    
%     for i=1:size(stimuliSizeCRF,1)
%         
%         if stimuliSizeCRF(i) == k
%             
%             allSizes(k).infoSecond(j) = InfoPerSecond1(i);
%             allSizes(k).infoPerSpike(j) = InfoPerSpike1(i);
%             allSizes(k).CodeEfficiency(j) = CodingEfficiency1(i);
%             j = j + 1;
%             
%         end
%         
%         if stimuliSizeNonCRF(i) == k
%             
%             allSizes(k).infoSecond(j) = InfoPerSecond2(i);
%             allSizes(k).infoPerSpike(j) = InfoPerSpike2(i);
%             allSizes(k).CodeEfficiency(j) = CodingEfficiency2(i);
%             j = j + 1;
%             
%         end
%             
%         if stimuliSizeExtraNonCRF(i) == k
%             
%             allSizes(k).infoSecond(j) = InfoPerSecond3(i);
%             allSizes(k).infoPerSpike(j) = InfoPerSpike3(i);
%             allSizes(k).CodeEfficiency(j) = CodingEfficiency3(i);
%             j = j + 1;
%             
%         end
%             
%         
%     end
%         
% end
% 
% % for k=1:8
% %     
% %    disp(strcat('Teste Lilliefors - ',int2str(k))); 
% %    h = lillietest(allSizes(k).infoSecond)
% %    h = lillietest(allSizes(k).infoPerSpike)
% %    h = lillietest(allSizes(k).CodeEfficiency)    
% %     
% % end
% 
% for k=2:8
%     
%     disp(int2str(k));
%     [h(k-1),p(k-1)] = ttest2(allSizes(k).infoSecond,allSizes(1).infoSecond)
%     [h(k-1),p(k-1)] = ttest2(allSizes(k).infoPerSpike,allSizes(1).infoPerSpike)
%     [h(k-1),p(k-1)] = ttest2(allSizes(k).CodeEfficiency,allSizes(1).CodeEfficiency)
%     
% end
% 
% for k=1:8
%    
%     meanAllinfoSecond(k) = mean(allSizes(k).infoSecond);
%     meanAllinfoPerSpike(k) = mean(allSizes(k).infoPerSpike);
%     meanAllCodeEfficiency(k) = mean(allSizes(k).CodeEfficiency);
%     
% end
% 
% 
% meanAllinfoSecond = meanAllinfoSecond.';
% meanAllinfoPerSpike = meanAllinfoPerSpike.';
% meanAllCodeEfficiency = meanAllCodeEfficiency.';
% 
% matrixAll = horzcat(meanAllinfoSecond,meanAllinfoPerSpike,meanAllCodeEfficiency);
% colnames = {'meanInfoSecond','meanInfoPerSpike','meanInfoCodeEfficiency'};
% t = uitable(f,'Data',matrixAll,'ColumnName',colnames,'Position',[960 10 300 750]);
% 
% disp('InfoPerSecond Group 2');
% 
% grupo2CRFInfoSecond = InfoPerSecond1([8 9 10 11 12 13 14 15 16 17 18 19]);
% grupo2NonCRFInfoSecond = InfoPerSecond2([8 9 10 11 12 13 14 15 16 17 18 19]);
% 
% mean(grupo2CRFInfoSecond)
% mean(grupo2NonCRFInfoSecond)
% 
% h2CRF = lillietest(grupo2CRFInfoSecond)
% h2NonCRF = lillietest(grupo2NonCRFInfoSecond)
% 
% [h,p] = ttest(grupo2CRFInfoSecond,grupo2NonCRFInfoSecond)
% 
% disp('InfoPerSpike Group 2');
% 
% grupo2CRFInfoSpike = InfoPerSpike1([8 9 10 11 12 13 14 15 16 17 18 19]);
% grupo2NonCRFInfoSpike = InfoPerSpike2([8 9 10 11 12 13 14 15 16 17 18 19]);
% 
% mean(grupo2CRFInfoSpike)
% mean(grupo2NonCRFInfoSpike)
% 
% h2CRF = lillietest(grupo2CRFInfoSpike)
% h2NonCRF = lillietest(grupo2NonCRFInfoSpike)
% 
% [h,p] = ttest(grupo2CRFInfoSpike,grupo2NonCRFInfoSpike)
% 
% disp('Code Efficiency Group 2');
% 
% grupo2CRFCodeEfficiency = CodingEfficiency1([8 9 10 11 12 13 14 15 16 17 18 19]);
% grupo2NonCRFCodeEfficiency = CodingEfficiency2([8 9 10 11 12 13 14 15 16 17 18 19]);
% 
% mean(grupo2CRFCodeEfficiency)
% mean(grupo2NonCRFCodeEfficiency)
% 
% h2CRF = lillietest(grupo2CRFCodeEfficiency)
% h2NonCRF = lillietest(grupo2NonCRFCodeEfficiency)
% 
% [h,p] = ttest(grupo2CRFInfoSecond,grupo2NonCRFCodeEfficiency)
% 
% disp('InfoPerSecond Ana');
% 
% grupoAnaCRFInfoPerSecond = [InfoPerSecond2(22) InfoPerSecond2(23) InfoPerSecond2(24) InfoPerSecond2(25) InfoPerSecond1(26) InfoPerSecond1(27) InfoPerSecond1(26) InfoPerSecond1(27) InfoPerSecond1(28) InfoPerSecond1(29) InfoPerSecond1(30) InfoPerSecond1(31) InfoPerSecond1(30) InfoPerSecond1(31) InfoPerSecond1(34) InfoPerSecond1(35) InfoPerSecond1(34) InfoPerSecond1(35) InfoPerSecond1(36) InfoPerSecond1(36) InfoPerSecond1(37) InfoPerSecond1(38) InfoPerSecond1(37) InfoPerSecond1(38) InfoPerSecond1(39) InfoPerSecond1(40) InfoPerSecond1(39) InfoPerSecond1(40) InfoPerSecond1(44) InfoPerSecond1(45) InfoPerSecond1(44) InfoPerSecond1(45) InfoPerSecond1(46) InfoPerSecond1(47)]; 
% grupoAnaNonCRFInfoPerSecond = [InfoPerSecond3(22) InfoPerSecond3(23) InfoPerSecond3(24) InfoPerSecond3(25) InfoPerSecond2(26) InfoPerSecond2(27) InfoPerSecond3(26) InfoPerSecond3(27) InfoPerSecond3(28) InfoPerSecond3(29) InfoPerSecond2(30) InfoPerSecond2(31) InfoPerSecond3(30) InfoPerSecond3(31) InfoPerSecond2(34) InfoPerSecond2(35) InfoPerSecond3(34) InfoPerSecond3(35) InfoPerSecond2(36) InfoPerSecond3(36) InfoPerSecond2(37) InfoPerSecond2(38) InfoPerSecond3(37) InfoPerSecond3(38) InfoPerSecond2(39) InfoPerSecond2(40) InfoPerSecond3(39) InfoPerSecond3(40) InfoPerSecond3(44) InfoPerSecond2(44) InfoPerSecond1(44) InfoPerSecond1(45) InfoPerSecond3(46) InfoPerSecond3(47)];                           
% 
% mean(grupoAnaCRFInfoPerSecond)
% mean(grupoAnaNonCRFInfoPerSecond)
% 
% median(grupoAnaCRFInfoPerSecond)
% median(grupoAnaNonCRFInfoPerSecond)
% 
% hAnaCRF = lillietest(grupoAnaCRFInfoPerSecond)
% hAnaNonCRF = lillietest(grupoAnaNonCRFInfoPerSecond)
% 
% [p,h] = signrank(grupoAnaCRFInfoPerSecond,grupoAnaNonCRFInfoPerSecond)
% 
% disp('InfoPerSpike Ana');
% 
% grupoAnaCRFInfoPerSpike = [InfoPerSpike2(22) InfoPerSpike2(23) InfoPerSpike2(24) InfoPerSpike2(25) InfoPerSpike1(26) InfoPerSpike1(27) InfoPerSpike1(26) InfoPerSpike1(27) InfoPerSpike1(28) InfoPerSpike1(29) InfoPerSpike1(30) InfoPerSpike1(31) InfoPerSpike1(30) InfoPerSpike1(31) InfoPerSpike1(34) InfoPerSpike1(35) InfoPerSpike1(34) InfoPerSpike1(35) InfoPerSpike1(36) InfoPerSpike1(36) InfoPerSpike1(37) InfoPerSpike1(38) InfoPerSpike1(37) InfoPerSpike1(38) InfoPerSpike1(39) InfoPerSpike1(40) InfoPerSpike1(39) InfoPerSpike1(40) InfoPerSpike1(44) InfoPerSpike1(45) InfoPerSpike1(44) InfoPerSpike1(45) InfoPerSpike1(46) InfoPerSpike1(47)]; 
% grupoAnaNonCRFInfoPerSpike = [InfoPerSpike3(22) InfoPerSpike3(23) InfoPerSpike3(24) InfoPerSpike3(25) InfoPerSpike2(26) InfoPerSpike2(27) InfoPerSpike3(26) InfoPerSpike3(27) InfoPerSpike3(28) InfoPerSpike3(29) InfoPerSpike2(30) InfoPerSpike2(31) InfoPerSpike3(30) InfoPerSpike3(31) InfoPerSpike2(34) InfoPerSpike2(35) InfoPerSpike3(34) InfoPerSpike3(35) InfoPerSpike2(36) InfoPerSpike3(36) InfoPerSpike2(37) InfoPerSpike2(38) InfoPerSpike3(37) InfoPerSpike3(38) InfoPerSpike2(39) InfoPerSpike2(40) InfoPerSpike3(39) InfoPerSpike3(40) InfoPerSpike3(44) InfoPerSpike2(44) InfoPerSpike1(44) InfoPerSpike1(45) InfoPerSpike3(46) InfoPerSpike3(47)];                           
%  
% mean(grupoAnaCRFInfoPerSpike)
% mean(grupoAnaNonCRFInfoPerSpike)
% 
% median(grupoAnaCRFInfoPerSpike)
% median(grupoAnaNonCRFInfoPerSpike)
% 
% hAnaCRF = lillietest(grupoAnaCRFInfoPerSpike)
% hAnaNonCRF = lillietest(grupoAnaNonCRFInfoPerSpike)
%  
% [p,h] = signrank(grupoAnaCRFInfoPerSpike,grupoAnaNonCRFInfoPerSpike)
% 
% disp('Coding Efficiency Ana');
% 
% grupoAnaCRFCodingEfficiency = [CodingEfficiency2(22) CodingEfficiency2(23) CodingEfficiency2(24) CodingEfficiency2(25) CodingEfficiency1(26) CodingEfficiency1(27) CodingEfficiency1(26) CodingEfficiency1(27) CodingEfficiency1(28) CodingEfficiency1(29) CodingEfficiency1(30) CodingEfficiency1(31) CodingEfficiency1(30) CodingEfficiency1(31) CodingEfficiency1(34) CodingEfficiency1(35) CodingEfficiency1(34) CodingEfficiency1(35) CodingEfficiency1(36) CodingEfficiency1(36) CodingEfficiency1(37) CodingEfficiency1(38) CodingEfficiency1(37) CodingEfficiency1(38) CodingEfficiency1(39) CodingEfficiency1(40) CodingEfficiency1(39) CodingEfficiency1(40) CodingEfficiency1(44) CodingEfficiency1(45) CodingEfficiency1(44) CodingEfficiency1(45) CodingEfficiency1(46) CodingEfficiency1(47)]; 
% grupoAnaNonCRFCodingEfficiency = [CodingEfficiency3(22) CodingEfficiency3(23) CodingEfficiency3(24) CodingEfficiency3(25) CodingEfficiency2(26) CodingEfficiency2(27) CodingEfficiency3(26) CodingEfficiency3(27) CodingEfficiency3(28) CodingEfficiency3(29) CodingEfficiency2(30) CodingEfficiency2(31) CodingEfficiency3(30) CodingEfficiency3(31) CodingEfficiency2(34) CodingEfficiency2(35) CodingEfficiency3(34) CodingEfficiency3(35) CodingEfficiency2(36) CodingEfficiency3(36) CodingEfficiency2(37) CodingEfficiency2(38) CodingEfficiency3(37) CodingEfficiency3(38) CodingEfficiency2(39) CodingEfficiency2(40) CodingEfficiency3(39) CodingEfficiency3(40) CodingEfficiency3(44) CodingEfficiency2(44) CodingEfficiency1(44) CodingEfficiency1(45) CodingEfficiency3(46) CodingEfficiency3(47)];                           
%  
% mean(grupoAnaCRFCodingEfficiency)
% mean(grupoAnaNonCRFCodingEfficiency)
% 
% median(grupoAnaCRFCodingEfficiency)
% median(grupoAnaNonCRFCodingEfficiency)
% 
% hAnaCRF = lillietest(grupoAnaCRFCodingEfficiency)
% hAnaNonCRF = lillietest(grupoAnaNonCRFCodingEfficiency)
%  
% [p,h] = signrank(grupoAnaCRFCodingEfficiency,grupoAnaNonCRFCodingEfficiency)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        SUPRESSAO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        SUPRESSAO
f = figure

for i=1:size(Protocolo,2)
    
    meanFRNonCRF(i) = Protocolo(i).rate.suppressionAnalysis.conditions(1,1).meanFiringRate;
    meanFRCRF(i) = Protocolo(i).rate.suppressionAnalysis.conditions(1,2).meanFiringRate;
    
    meanRateModulationNonCRF(i) = Protocolo(i).rate.suppressionAnalysis.meanRateModulationNonCRF;
    
    for k=1:300

        pro(i).meanRateBinCRF(k) = Protocolo(i).rate.suppressionAnalysis.meanRateBinCRF(k);
        pro(i).meanRateBinNon(k) = Protocolo(i).rate.suppressionAnalysis.meanRateBinNon(k);
        pro(i).meanRateBinExtraNon(k) = 0;

    end


    if i>=10
    
        meanFRCRF(i) = Protocolo(i).rate.suppressionAnalysis.conditions(1,1).meanFiringRate;
        meanFRNonCRF(i) = Protocolo(i).rate.suppressionAnalysis.conditions(1,2).meanFiringRate;
        
    end

    if i>=20
       
        meanFRExtraNonCRF(i) = Protocolo(i).rate.suppressionAnalysis.conditions(1,3).meanFiringRate;;
        
        meanRateModulationExtraNonCRF(i) = Protocolo(i).rate.suppressionAnalysis.meanRateModulationExtraNonCRF;
        
        for k=1:300

            pro(i).meanRateBinExtraNon(k) = Protocolo(i).rate.suppressionAnalysis.meanRateBinExtraNon(k);

        end

    end
    
end

meanFRCRF = meanFRCRF.';
meanFRNonCRF = meanFRNonCRF.';
meanFRExtraNonCRF = meanFRExtraNonCRF.';

meanRateModulationNonCRF = meanRateModulationNonCRF.';
meanRateModulationExtraNonCRF = meanRateModulationExtraNonCRF.';

h = lillietest(meanFRCRF)
h = lillietest(meanFRNonCRF)
h = lillietest(meanFRExtraNonCRF)

for i=1:size(Protocolo,2)
   
    [pFRNon(i),hFRNon(i)] = signrank(meanFRCRF,meanFRNonCRF(i));
    
    pFRExtraNon(i) = -1;
    hFRExtraNon(i) = -1;
    
    if i>=20
       
        [pFRExtraNon(i),hFRExtraNon(i)] = signrank(meanFRCRF,meanFRExtraNonCRF(i));
        
    end
    
end

pFRNon = pFRNon.';
hFRNon = hFRNon.';
pFRExtraNon = pFRExtraNon.';
hFRExtraNon = hFRExtraNon.';

mean(meanFRCRF)


for i=1:size(Protocolo,2)

hlillieRMCRF(i) = lillietest(pro(i).meanRateBinCRF);
hlillieRMNon(i) = lillietest(pro(i).meanRateBinNon);

[pRMNon(i),hRMNon(i)] = signrank(pro(i).meanRateBinCRF,pro(i).meanRateBinNon);
    
    hlillieRMExtraNon(i) = -1;
    pRMExtra(i) = -1;
    hRMExtra(i) = -1;

    if i>=20

        hlillieRMExtraNon(i) = lillietest(pro(i).meanRateBinExtraNon);

        [pRMExtra(i),hRMExtra(i)] = signrank(pro(i).meanRateBinCRF,pro(i).meanRateBinExtraNon);
    
    end

end

hlillieRMCRF = hlillieRMCRF.';
hlillieRMNon = hlillieRMNon.';
hlillieRMExtraNon = hlillieRMExtraNon.';
hRMNon = hRMNon.';
pRMNon = pRMNon.';
pRMExtra = pRMExtra.';
hRMExtra = hRMExtra.';

FiringRateMatrix = horzcat(meanFRCRF,meanFRNonCRF,pFRNon,hFRNon,meanFRExtraNonCRF,pFRExtraNon,hFRExtraNon);

RateModulationMatrix = horzcat(hlillieRMCRF,hlillieRMNon,hlillieRMExtraNon,hRMNon,pRMNon,pRMExtra,hRMExtra,meanRateModulationNonCRF,meanRateModulationExtraNonCRF);

colnamesMeanFR = {'meanFRCRF','meanFRNonCRF','pFRNon','hFRNon','meanFRExtraNonCRF','pFRExtraNon','hFRExtraNon'};

colnamesMeanRM = {'hlillieRMCRF','hlillieRMNon','hlillieRMExtraNon','hRMNon','pRMNon','pRMExtra','hRMExtra','RateMod-NonCRF','RateMod-ExtraNonCRF'};

t = uitable(f,'Data',FiringRateMatrix,'ColumnName',colnamesMeanFR,'Position',[10 10 500 700]);

t = uitable(f,'Data',RateModulationMatrix,'ColumnName',colnamesMeanRM,'Position',[520 10 700 700]);


RateModulation = [meanRateModulationNonCRF meanRateModulationExtraNonCRF];

figure;
bar(RateModulation);

% video(1) = 1;
% video(2) = 3;
% video(3) = 2;
% video(4) = 4;
% video(5) = 3;
% video(6) = 2;
% video(7) = 3;
% video(8) = 5;
% video(9) = 6;
% video(10) = 5;
% video(11) = 6;
% video(12) = 5;
% video(13) = 6;
% video(14) = 5;
% video(15) = 6;
% video(16) = 1;
% video(17) = 4;
% video(18) = 3;
% video(19) = 6;
% video(20) = 2;
% video(21) = 1;
% video(22) = 8;
% video(23) = 11;
% video(24) = 9;
% video(25) = 10;
% video(26) = 6;
% video(27) = 8;
% video(28) = 3;
% video(29) = 8;
% video(30) = 5;
% video(31) = 8;
% video(32) = 3;
% video(33) = 4;
% video(34) = 8;
% video(35) = 2;
% 
% contrast1(1) = 1.1585;
% contrast1(2) = 1.7122;
% contrast1(3) = 1.5632;
% contrast1(4) = 1.9120;
% contrast1(5) = 1.7122;
% contrast1(6) = 1.5632;
% contrast1(7) = 1.7122;
% contrast1(8) = 1.7620;
% contrast1(9) = 1.8462;
% contrast1(10) = 1.7620;
% contrast1(11) = 1.5313;
% contrast1(12) = 1.5014;
% contrast1(13) = 1.5313;
% contrast1(14) = 1.5014;
% contrast1(15) = 1.5313;
% contrast1(16) = 0.9072;
% contrast1(17) = 1.6924;
% contrast1(18) = 1.4985;
% contrast1(19) = 1.5313;
% contrast1(20) = 1.1631;
% contrast1(21) = 0.9072;
% contrast1(22) = 1.6485;
% contrast1(23) = 1.2140;
% contrast1(24) = 1.1043;
% contrast1(25) = 1.5110;
% contrast1(26) = 1.5313;
% contrast1(27) = 1.6485;
% contrast1(28) = 1.6322;
% contrast1(29) = 1.7270;
% contrast1(30) = 1.6337;
% contrast1(31) = 1.7753;
% contrast1(32) = 1.6273;
% contrast1(33) = 1.7823;
% contrast1(34) = 1.6953;
% contrast1(35) = 1.3035;
% 
% contrast2(1) = 1.9442;
% contrast2(2) = 1.9212;
% contrast2(3) = 1.9333;
% contrast2(4) = 2;
% contrast2(5) = 1.9212;
% contrast2(6) = 1.9333;
% contrast2(7) = 1.9212;
% contrast2(8) = 1.8881;
% contrast2(9) = 1.9235;
% contrast2(10) = 1.8765;
% contrast2(11) = 1.9281;
% contrast2(12) = 1.8765;
% contrast2(13) = 1.9281;
% contrast2(14) = 1.8765;
% contrast2(15) = 1.9281;
% contrast2(16) = 1.8585;
% contrast2(17) = 1.9960;
% contrast2(18) = 1.8727;
% contrast2(19) = 1.9281;
% contrast2(20) = 1.8037;
% contrast2(21) = 1.8383;
% contrast2(22) = 1.9066;
% contrast2(23) = 1.6203;
% contrast2(24) = 1.9163;
% contrast2(25) = 1.8636;
% contrast2(26) = 1.9281;
% contrast2(27) = 1.9042;
% contrast2(28) = 0;
% contrast2(29) = 0;
% contrast2(30) = 1.8431;
% contrast2(31) = 1.8388;
% contrast2(32) = 1.8715;
% contrast2(33) = 1.9879;
% contrast2(34) = 1.8740;
% contrast2(35) = 1.7276;
% 
% contrast3(22) = 2;
% contrast3(23) = 2;
% contrast3(24) = 2;
% contrast3(25) = 2;
% contrast3(26) = 2;
% contrast3(27) = 2;
% contrast3(28) = 1.8727;
% contrast3(29) = 1.9066;
% contrast3(30) = 1.8765;
% contrast3(31) = 1.9066;
% contrast3(32) = 1.9212;
% contrast3(33) = 2;
% contrast3(34) = 1.9392;
% contrast3(35) = 1.9333;
% 
% video = video.';
% contrast1 = contrast1.';
% contrast2 = contrast2.';
% contrast3 = contrast3.';
% 
% videoMatrix = horzcat(video,contrast1,contrast2,contrast3);
% 
% colnamesVideo = {'Video','Contrast-1','Contrast-2','Contrast-3'};
% 
% %t = uitable(f,'Data',videoMatrix,'ColumnName',colnamesVideo,'Position',[10 510 350 270]);
% 





end