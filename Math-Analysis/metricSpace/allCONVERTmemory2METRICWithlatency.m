function allCONVERTmemory2METRICWithlatency(start_time,end_time)

%convertFullMovieSpass2MetricWithlatency(date,site_index,channel,registro,video_index,start_time,end_time,invert_movie,bin_size,nRepetitions,ae)

for i=0:1
    
    %%%  primeiro registro

    convertFullMovieSpass2MetricWithlatency('12-08-20',1,'E1','nsp003a01_1b',3,0,4000,i,1,40,0);

    %%%%%%%%%%%%%%%%%%

    %%% primeira semana

    convertFullMovieSpass2MetricWithlatency('12-08-27',1,'E1','nsp005a01_1b',3,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-08-27',1,'E2','nsp005a01_2b',3,start_time,end_time,i,1,60,0);

    convertFullMovieSpass2MetricWithlatency('12-08-29',1,'E1','nsp006a01_1b',3,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-08-29',1,'E1','nsp006a02_1b',4,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-08-29',2,'E1','nsp006b01_1b',1,start_time,end_time,i,1,60,0);

    convertFullMovieSpass2MetricWithlatency('12-08-30',1,'E1','nsp007a01_1b',4,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-08-30',1,'E3','nsp007a01_2b',4,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-08-30',1,'E1','nsp007a02_1b',1,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-08-30',1,'E3','nsp007a02_2b',1,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-08-30',1,'E3','nsp007a03_2a',3,start_time,end_time,i,1,60,0);

    convertFullMovieSpass2MetricWithlatency('12-08-31',1,'E1','nsp008a01_1b',1,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-08-31',1,'E3','nsp008a01_2b',1,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-08-31',1,'E1','nsp008a02_1a',3,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-08-31',1,'E3','nsp008a02_2b',3,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-08-31',1,'E1','nsp008a03_1b',4,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-08-31',1,'E3','nsp008a03_2a',4,start_time,end_time,i,1,60,0);

    %%%%%%%%%%%%%%%%%%%%%%


    %%% segunda semana

    convertFullMovieSpass2MetricWithlatency('12-09-03',1,'E3','nsp009a01_2b',2,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-09-03',2,'E1','nsp009b01_1a',4,start_time,end_time,i,1,60,0);

    convertFullMovieSpass2MetricWithlatency('12-09-04',1,'E1','nsp010a1_1b',2,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-09-04',1,'E1','nsp010a02_1b',4,start_time,end_time,i,1,60,0);

    convertFullMovieSpass2MetricWithlatency('12-09-05',1,'E1','nsp011a01_1a',3,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-09-05',1,'E1','nsp011a02_1b',4,start_time,end_time,i,1,60,0);

    convertFullMovieSpass2MetricWithlatency('12-09-06',1,'E1','nsp012a01_1b',3,start_time,end_time,i,1,60,0);
    convertFullMovieSpass2MetricWithlatency('12-09-06',1,'E1','nsp012a02_1a',4,start_time,end_time,i,1,60,0);

    %%%%%%%%%%%%%%%%%%%%

    convertFullMovieSpass2MetricWithlatency('12-10-16',1,'E1','nps013a01_1b',1,start_time,end_time,i,1,60,0);
    
    convertFullMovieSpass2MetricWithlatency('12-10-17',1,'E1','nsp033a09_1b',311,start_time,end_time,i,1,60,0);

    convertFullMovieSpass2MetricWithlatency('12-10-17',1,'E1','nsp033a09_1c',321,start_time,end_time,i,1,60,0);

    convertFullMovieSpass2MetricWithlatency('12-10-17',1,'E3','nsp033a09_3b',331,start_time,end_time,i,1,60,0);
    
end

 %%%  primeiro registro

    convertFullMovieSpass2MetricWithlatency('12-08-20',1,'E1','nsp003a01_1b',3,0,4000,0,1,40,1);

    %%%%%%%%%%%%%%%%%%

    %%% primeira semana

    convertFullMovieSpass2MetricWithlatency('12-08-27',1,'E1','nsp005a01_1b',3,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-08-27',1,'E2','nsp005a01_2b',3,start_time,end_time,0,1,60,1);

    convertFullMovieSpass2MetricWithlatency('12-08-29',1,'E1','nsp006a01_1b',3,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-08-29',1,'E1','nsp006a02_1b',4,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-08-29',2,'E1','nsp006b01_1b',1,start_time,end_time,0,1,60,1);

    convertFullMovieSpass2MetricWithlatency('12-08-30',1,'E1','nsp007a01_1b',4,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-08-30',1,'E3','nsp007a01_2b',4,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-08-30',1,'E1','nsp007a02_1b',1,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-08-30',1,'E3','nsp007a02_2b',1,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-08-30',1,'E3','nsp007a03_2a',3,start_time,end_time,0,1,60,1);

    convertFullMovieSpass2MetricWithlatency('12-08-31',1,'E1','nsp008a01_1b',1,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-08-31',1,'E3','nsp008a01_2b',1,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-08-31',1,'E1','nsp008a02_1a',3,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-08-31',1,'E3','nsp008a02_2b',3,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-08-31',1,'E1','nsp008a03_1b',4,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-08-31',1,'E3','nsp008a03_2a',4,start_time,end_time,0,1,60,1);

    %%%%%%%%%%%%%%%%%%%%%%


    %%% segunda semana

    convertFullMovieSpass2MetricWithlatency('12-09-03',1,'E3','nsp009a01_2b',2,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-09-03',2,'E1','nsp009b01_1a',4,start_time,end_time,0,1,60,1);

    convertFullMovieSpass2MetricWithlatency('12-09-04',1,'E1','nsp010a1_1b',2,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-09-04',1,'E1','nsp010a02_1b',4,start_time,end_time,0,1,60,1);

    convertFullMovieSpass2MetricWithlatency('12-09-05',1,'E1','nsp011a01_1a',3,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-09-05',1,'E1','nsp011a02_1b',4,start_time,end_time,0,1,60,1);

    convertFullMovieSpass2MetricWithlatency('12-09-06',1,'E1','nsp012a01_1b',3,start_time,end_time,0,1,60,1);
    convertFullMovieSpass2MetricWithlatency('12-09-06',1,'E1','nsp012a02_1a',4,start_time,end_time,0,1,60,1);

    %%%%%%%%%%%%%%%%%%%%

    convertFullMovieSpass2MetricWithlatency('12-10-16',1,'E1','nps013a01_1b',1,start_time,end_time,0,1,60,1);
    
    convertFullMovieSpass2MetricWithlatency('12-10-17',1,'E1','nsp033a09_1b',311,start_time,end_time,0,1,60,1);

    convertFullMovieSpass2MetricWithlatency('12-10-17',1,'E1','nsp033a09_1c',321,start_time,end_time,0,1,60,1);

    convertFullMovieSpass2MetricWithlatency('12-10-17',1,'E3','nsp033a09_3b',331,start_time,end_time,0,1,60,1);

end