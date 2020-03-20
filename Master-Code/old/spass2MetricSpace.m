function spass2MetricSpace(nome_do_registro,inverter_movie,latency_file,WOlatency,Conditions,nome_metric_file)

%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function loads a data set of electrophysiology coming from Spass
% software format and transfom the data set to the format used as input
% for the STATTOOLKIT Toolbox[1].
%
% INPUT:
%      1. spassDataSetName: the name of Spass data set file.      
%
% OUTPUT:
%      1.sparseness: structure with three variables:
%
% [1] Vinje, W., & Gallant, J. (2000). Sparse coding and decorrelation in 
% primary visual cortex during natural vision. Science (New York, NY), 287(5456), 1273.
%
% [2] Vinje, W., & Gallant, J. (2002). Natural stimulation of the nonclassical receptive 
% field increases information transmission efficiency in V1. Journal of Neuroscience, 22(7), 2904.
%
%
% Author: Joao V. Dornas, joaodornas@gmail.com, 03/2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% PARAMETROS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nome_do_registro = colocar o nome do registro sem a extencao '.mat'

% inverter_movie = define se quer inverter a segunda condicao ou nao
%
%    0 - nao inverte a segunda condicao
%    1 - inverte a segunda condicao
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

latencies = load(char(latency_file));
latency = latencies.latency_time;

disp('BEGIN - spass2metric');

%%%   ABRE REGISTRO E L? PAR?METROS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Abre arquivo do registro
Spass = load(strcat(nome_do_registro,'.mat'));
 
%L? n?mero de condi??es
nConditions = Spass.parameters.nconditions;

%Separa os Trials Labels para cada condi??o
for i=1:nConditions
    
    trials_label(i).labels = find(Spass.stimIds == i); 
 
end

%L? o n?mero de trials
nTrials = length(trials_label(1).labels);

%Ajusta os valores do tempo de cada spike
spike_times = Spass.spike_times./ 32000;

%L? tempo inciail e final do Trial
start_time = Spass.parameters.baseline_duration;
end_time = Spass.parameters.trialsize - Spass.parameters.baseline_duration;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   DEFINE PAR?METROS DO METRIC SPACE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Escala de tempo em segundos
time_scale = 1;

%Escala da resolucao temporal em segundos
time_resolution = 0.0001;

%Sistema Internacional de Medidas
si_prefix = 1;

%Numero de classes de estimulos
%M = nConditions;
M = length(Conditions);

if nConditions < M
    
    M = nConditions;
    
end

%Numero de sitios (canais)
N = 1;

%Numero de trials (repeticoes por condicao)
P = nTrials;

%Numero de pontos no vetor list
%Q = 0;

%Cria vetor de categorias para as condicoes
%categories(1:nConditions) = struct('label','','P',P,'trials',zeros(P));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   PROCESSA AS CATEGORIAS/CONDI??ES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:nConditions
  
    %L? conjunto de trials especifico para uma condicao
    trials_spikes = spike_times(trials_label(i).labels(:),:);
    
    nTrials = size(trials_spikes,1);

    %Define o vetor trials no formato do Metric Space
    trials(N,1:nTrials) = struct('start_time',start_time/1000,'end_time',end_time/1000,'Q',0,'list',zeros(size(trials_spikes,2)));

    for k=1:nTrials

        spikes = trials_spikes(k,:);
        
        if WOlatency
               
            spikes = spikes - latency(i);
               
        end

        spikes = spikes(spikes>0);

        spikes = sort(spikes);

        spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));

        %Se for a segunda condicao, e for solicitado que se inverta as
        %trilhas, inverte a trilha
        if (i == 2) && (inverter_movie == 1)

            nBins = (end_time - start_time);

            spikes = inverteTrilha(spikes,start_time,nBins);
            
        end

        %Salva a trilha no formato do Metric Space
        trials(N,k).Q = size(spikes,2);
        trials(N,k).list = spikes;

    end

    trials = trials';

    %Define nome da Categoria
    category_name{i} = int2str(i);
    
    %Salva as informacoes da condicao como categoria para o Metric Space
    the_category.label = category_name{i};
    the_category.P = nTrials;
    the_category.trials = trials;

    clear trials;
    
    categories(i) = the_category;
    
    clear the_category;
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CRIA ARQUIVO DE DADOS DO METRIC SPACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

categories = categories';

sites = struct('label','unit_001','recording_tag','episodic','time_scale',time_scale,'time_resolution',time_resolution,'si_unit','none','si_prefix',si_prefix);

X = struct('M',M,'N',N,'sites',sites,'categories',categories);

%Salva as categorias em um arquivo '.stam' e '.stad' no formato Metric
%Space
stawrite(X,nome_metric_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END - spass2metric');

end

