function metric1xManyx

% gallant2002Sparseness('12-11-01',1,'E1','nsp020a01_1b',5,500,9500,bin_size,2);
% gallant2002Sparseness('12-11-01',1,'E1','nsp020a02_1b',6,500,9500,bin_size,2);
% 
% gallant2002Sparseness('12-11-06',1,'E1','nsp021a01_1a',5,500,9500,bin_size,2);
% gallant2002Sparseness('12-11-06',1,'E1','nsp021a02_1a',6,500,9500,bin_size,2);
% 
% gallant2002Sparseness('12-11-07',1,'E1','nsp022a01_1b',5,500,9500,bin_size,2);
% gallant2002Sparseness('12-11-07',1,'E1','nsp022a02_1b',6,500,9500,bin_size,2);
% 
% gallant2002Sparseness('12-11-08',1,'E1','nsp023a01_1a',1,500,9500,bin_size,2);
% gallant2002Sparseness('12-11-08',1,'E1','nsp023a02_1a',4,500,9500,bin_size,2);
% 
% gallant2002Sparseness('12-11-09',1,'E1','nsp024a01_1a',6,500,9500,bin_size,2);
% gallant2002Sparseness('12-11-09',1,'E1','nsp024a02_1a',3,500,9500,bin_size,2);
% 
% gallant2002Sparseness('12-11-13',1,'E1','nsp025a01_1b',2,500,9500,bin_size,2);
% gallant2002Sparseness('12-11-13',1,'E1','nsp025a02_1a',1,500,9500,bin_size,2);

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

% reliabilityCC('12-12-03',1,'E3','nsp026a01_3b',8,500,9500,30,15,3);
% reliabilityCC('12-12-03',1,'E3','nsp026a02_3a',11,500,9500,30,15,3);
% 
% reliabilityCC('12-12-04',1,'E1','nsp027a03_1a',10,500,9500,30,15,3);
% reliabilityCC('12-12-04',1,'E1','nsp027a04_1b',9,500,9500,30,15,3);
% 
% reliabilityCC('12-12-05',1,'E1','nsp028a03_1a',6,500,9500,30,15,3);
% reliabilityCC('12-12-05',1,'E1','nsp028a04_1b',8,500,9500,30,15,3);
% 
% reliabilityCC('12-12-07',1,'E1','nsp029a03_1a',8,500,9500,30,15,3);
% reliabilityCC('12-12-07',1,'E1','nsp029a04_1b',3,500,9500,30,15,3);
% 
% reliabilityCC('12-12-10',1,'E2','nsp030a03_2a',5,500,9500,30,15,3);
% reliabilityCC('12-12-10',1,'E2','nsp030a04_2b',8,500,9500,30,15,3);
%  
% reliabilityCC('12-12-12',2,'E1','nsp031b03_1b',3,500,9500,30,15,3);
% reliabilityCC('12-12-12',2,'E1','nsp031b04_1b',4,500,9500,30,15,3);
% 
% reliabilityCC('12-12-13',1,'E1','nsp032a03_1a',8,500,9500,30,15,3);
% reliabilityCC('12-12-13',1,'E1','nsp032a04_1a',2,500,9500,30,15,3);
%  
% reliabilityCC('12-12-17',1,'E1','nsp033a04_1c',811,500,9500,30,15,3);
% reliabilityCC('12-12-17',1,'E1','nsp033a06_1b',812,500,9500,30,15,3);
% 
% reliabilityCC('12-12-17',1,'E1','nsp033a04_1b',821,500,9500,30,15,3);
% reliabilityCC('12-12-17',1,'E1','nsp033a06_1c',822,500,9500,30,15,3);
% 
% reliabilityCC('12-12-17',1,'E1','nsp033a05_1b',511,500,9500,30,15,3);
% reliabilityCC('12-12-17',1,'E1','nsp033a07_1b',512,500,9500,30,15,3);
% 
% reliabilityCC('12-12-17',1,'E1','nsp033a05_1c',521,500,9500,30,15,3);
% 
% reliabilityCC('12-12-17',1,'E3','nsp033a04_3a',831,500,9500,30,15,3);
% reliabilityCC('12-12-17',1,'E3','nsp033a05_3a',531,500,9500,30,15,3);
% reliabilityCC('12-12-17',1,'E3','nsp033a07_3b',532,500,9500,30,15,3);
% 
% reliabilityCC('12-12-18',1,'E2','nsp034a04_2b',8,500,9500,30,15,3);
% reliabilityCC('12-12-18',1,'E2','nsp034a05_2b',11,500,9500,30,15,3);
% 
% reliabilityCC('12-12-19',1,'E2','nsp035a03_2b',8,500,9500,30,15,3);
% reliabilityCC('12-12-19',1,'E2','nsp035a04_2b',1,500,9500,30,15,3);

%registrosGrupo2 = {'nsp020a01_1b' 'nsp020a02_1b' 'nsp021a01_1a' 'nsp021a02_1a' 'nsp022a01_1b' 'nsp022a02_1b' 'nsp023a01_1a' 'nsp023a02_1a' 'nsp024a01_1a' 'nsp024a02_1a' 'nsp025a01_1b' 'nsp025a02_1a'};
%video_index = [5 6 5 6 5 6 1 4 6 3 2 1];

registrosGrupoAna = {'nsp026a01_3b' 'nsp026a02_3a' 'nsp027a03_1a' 'nsp027a04_1b' 'nsp028a03_1a' 'nsp028a04_1b' 'nsp029a03_1a' 'nsp029a04_1b' 'nsp030a03_2a' 'nsp030a04_2b' 'nsp031b03_1b' 'nsp031b04_1b' 'nsp032a03_1a' 'nsp032a04_1a' 'nsp033a04_1c' 'nsp033a06_1b' 'nsp033a04_1b' 'nsp033a06_1c' 'nsp033a05_1b' 'nsp033a07_1b' 'nsp033a05_1c' 'nsp033a04_3a' 'nsp033a05_3a' 'nsp033a07_3b' 'nsp034a04_2b' 'nsp034a05_2b' 'nsp035a03_2b' 'nsp035a04_2b'};
video_index = [8 11 10 9 6 8 8 3 5 8 3 4 8 2 811 812 821 822 511 512 521 831 531 532 8 11 8 1];
condition1 = [2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; 
condition2 = [3 3 3 3 2 2 3 3 3 3 2 2 3 3 2 2 3 3 2 3 2 2 3 3 2 2 3 3 2 2 1 1 3 3];                           

%%%   SET METRIC SPACE PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Set Metric Space parameters...');

%Escala de tempo em segundos
time_scale = 1;

%Escala da resolu??o temporal em segundos
time_resolution = 0.0001;

%Sistema Internacional de Medidas
si_prefix = 1;

%N?mero de classes de est?mulos
M = size(registrosGrupoAna,2);

%N?mero de s?tios (canais)
N = 1;

%N?mero de trials (repeti??es por condi??o)
P = 60;

%Tempos inicial e final do trial
start_time = 500;
end_time = 9500;

%N?mero de pontos no vetor list
%Q = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

categories1x(1:M) = struct('label','','P',P,'trials',zeros(P));
categoriesMx(1:M) = struct('label','','P',P,'trials',zeros(P));

for i=1:size(registrosGrupoAna,2)
    
    %%%   LOAD TRIALS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Load trials...');

    Spass = load(strcat('_',registrosGrupoAna{i},'-','v',int2str(video_index(i)),'.mat'));

    for j=1:3

        trials_label(j).labels = find(Spass.stimIds == j); 

    end

    Spass.spike_times = Spass.spike_times./ 32000;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%   BEGIN CATEGORIES LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Begin categories loop...');

    for j=1:2

        if j == 1
            
            c = condition1(i);
        
        else
            
            c = condition2(i);
            
        end
        
        disp(strcat('Category .',int2str(c)));

        disp('...get trials');

        trials_spikes = Spass.spike_times(trials_label(c).labels(:),:);

        nTrials = size(trials_spikes,1);

        trials(N,1:nTrials) = struct('start_time',start_time/1000,'end_time',end_time/1000,'Q',0,'list',zeros(size(trials_spikes,2)));

        for k=1:nTrials

            spikes = trials_spikes(k,:);
            spikes = spikes(spikes>0);

            spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));
            spikes = sort(spikes);
            
            disp(strcat('...set trial .',int2str(k)));

            trials(N,k).Q = size(spikes,2);
            trials(N,k).list = spikes;

        end

        trials = trials';

        if j == 1
            
            categories1x(i).label = registrosGrupoAna(i);
            categories1x(i).P = nTrials;
            categories1x(i).trials = trials
            
        else
            
            categoriesMx(i).label = registrosGrupoAna{i};
            categoriesMx(i).P = nTrials;
            categoriesMx(i).trials = trials
            
        end

        clear trials;

    end

    disp('End categories loop...');

end

%%% CREATE METRIC SPACE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Create Metric Space data...');

categories1x = categories1x';
categoriesMx = categoriesMx';

sites = struct('label','unit_001','recording_tag','episodic','time_scale',time_scale,'time_resolution',time_resolution,'si_unit','none','si_prefix',si_prefix);

X1x = struct('M',M,'N',N,'sites',sites,'categories',categories1x);
XMx = struct('M',M,'N',N,'sites',sites,'categories',categoriesMx);

%X = struct('M',M,'N',N,'sites',struct('label','unit_001','recording_tag','episodic','time_scale',time_scale,'time_resolution',time_resolution,'si_unit','none','si_prefix',si_prefix),'categories',struct('label','','P',P,'trials',struct('start_time',start_time,'end_time',end_time,'Q',Q,'list',list)));

filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/');

filepath1x = strcat(filepath,'metric1x');
filepathMx = strcat(filepath,'metricMx');

stawrite(X1x,filepath1x);
stawrite(XMx,filepathMx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END');

end

