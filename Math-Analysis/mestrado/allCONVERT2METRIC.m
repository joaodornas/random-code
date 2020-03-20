function h = allCONVERT2METRIC(invert_movie,bin_size,start_time,end_time)

%%%  primeiro registro

convertFullMovieSpass2Metric('12-08-20',1,'E1','nsp003a01_1b',3,0,4000,invert_movie,bin_size,40);

%%%%%%%%%%%%%%%%%%


%%% primeira semana

convertFullMovieSpass2Metric('12-08-27',1,'E1','nsp005a01_1b',3,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-08-27',1,'E2','nsp005a01_2b',3,start_time,end_time,invert_movie,bin_size,60);

convertFullMovieSpass2Metric('12-08-29',1,'E1','nsp006a01_1b',3,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-08-29',1,'E1','nsp006a02_1b',4,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-08-29',2,'E1','nsp006b01_1b',1,start_time,end_time,invert_movie,bin_size,60);

convertFullMovieSpass2Metric('12-08-30',1,'E1','nsp007a01_1b',4,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-08-30',1,'E3','nsp007a01_2b',4,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-08-30',1,'E1','nsp007a02_1b',1,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-08-30',1,'E3','nsp007a02_2b',1,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-08-30',1,'E3','nsp007a03_2a',3,start_time,end_time,invert_movie,bin_size,60);

convertFullMovieSpass2Metric('12-08-31',1,'E1','nsp008a01_1b',1,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-08-31',1,'E3','nsp008a01_2b',1,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-08-31',1,'E1','nsp008a02_1a',3,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-08-31',1,'E3','nsp008a02_2b',3,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-08-31',1,'E1','nsp008a03_1b',4,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-08-31',1,'E3','nsp008a03_2a',4,start_time,end_time,invert_movie,bin_size,60);

%%%%%%%%%%%%%%%%%%%%%%


%%% segunda semana

convertFullMovieSpass2Metric('12-09-03',1,'E3','nsp009a01_2b',2,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-09-03',2,'E1','nsp009b01_1a',4,start_time,end_time,invert_movie,bin_size,60);

convertFullMovieSpass2Metric('12-09-04',1,'E1','nsp010a02_1b',4,start_time,end_time,invert_movie,bin_size,60);

convertFullMovieSpass2Metric('12-09-05',1,'E1','nsp011a01_1a',3,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-09-05',1,'E1','nsp011a02_1b',4,start_time,end_time,invert_movie,bin_size,60);

convertFullMovieSpass2Metric('12-09-06',1,'E1','nsp012a01_1b',3,start_time,end_time,invert_movie,bin_size,60);
convertFullMovieSpass2Metric('12-09-06',1,'E1','nsp012a02_1a',4,start_time,end_time,invert_movie,bin_size,60);

%%%%%%%%%%%%%%%%%%%%

convertFullMovieSpass2Metric('12-10-16',1,'E1','nps013a01_1b',1,start_time,end_time,invert_movie,bin_size,60);



end