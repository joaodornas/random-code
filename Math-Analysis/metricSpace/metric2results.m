function metric2results(nome_do_registro)

%%%%%% PARAMETROS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nome_do_registro = colocar o nome do registro sem a extencao '.mat'
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

disp('BEGIN - metric2results');

%%%   L? NUMERO DE CONDICOES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Abre arquivo do registro
Spass = load(strcat(nome_do_registro,'.mat'));
 
%L? n?mero de condi??es
nConditions = Spass.parameters.nconditions;

%L? tempos inicial e final da estimulacao
start_time = Spass.parameters.baseline_duration;
end_time = Spass.parameters.trialsize - Spass.parameters.baseline_duration;

clear Spass;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% ABRE DADOS DO METRIC SPACE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = staread(strrep(strcat(nome_do_registro,'.stam'),'/',filesep));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE OPCOES DO METRIC SPACE    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define tipos de entropia a calcular de acordo com o metodo de correcao
opts.entropy_estimation_method = {'plugin','tpmc','jack'};
%opts.variance_estimation_method = {'jack'};

%Sobre bins vazios
%opts.unoccupied_bins_strategy = -1; % Ignore unoccupied bins
opts.unoccupied_bins_strategy = 0; % Use an unoccupied bin only if its row and column are occupied
%opts.unoccupied_bins_strategy = 1; % Use all bins

%Calcula as dist?ncias para todos os 'shift_cost' ao mesmo tempo
opts.parallel = 1;

%Estrat?gia para calcular o n?mero de palavras poss?veis
opts.possible_words = 'unique';

%Tempos inicial e final da estimulacao
opts.start_time = start_time / 1000;
opts.end_time = end_time / 1000;

%Resolucoes a serem calculadas
opts.shift_cost = [0 2.^(-4:9)];

%Expoente da fun??o que calcula a m?dia da dist?ncias de uma trilha
opts.clustering_exponent = -2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DEFINE NUMERO S PARA CONDICOES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = nConditions;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CALCULA METRIC SPACE - CODIGO TIMING   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define c?digo Timing
opts.metric_family = 0;

%Calcula Metric Space
[out_T,opts_used] = metric(X,opts);

%L? valores de Informa??o para cada tipo de entropia
for q_idx=1:length(opts.shift_cost)

    info_plugin_T(q_idx) = out_T(q_idx).table.information(1).value;
    info_tpmc_T(q_idx) = out_T(q_idx).table.information(2).value;
    info_jack_T(q_idx) = out_T(q_idx).table.information(3).value;

end

%L? valores m?ximos de Informa??o para cada tipo de entropia
[max_info_plugin_T,max_info_plugin_idx_T] = max(info_plugin_T);
[max_info_tpmc_T,max_info_tpmc_idx_T] = max(info_tpmc_T);
[max_info_jack_T,max_info_jack_idx_T] = max(info_jack_T);

%L? valores de Informa??o para cada tipo de entropia para c?digo Taxa
Hcount_plugin_T = info_plugin_T(1);
Hcount_tpmc_T = info_tpmc_T(1);
Hcount_jack_T = info_jack_T(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% CALCULA METRIC SPACE - CODIGO INTERVAL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define c?digo Interval
opts.metric_family = 1;

%Calcula Metric Space
[out_I,opts_used] = metric(X,opts);

%L? valores de Informa??o para cada tipo de entropia e resolu??o
for q_idx=1:length(opts.shift_cost)

    info_plugin_I(q_idx) = out_I(q_idx).table.information(1).value;
    info_tpmc_I(q_idx) = out_I(q_idx).table.information(2).value;
    info_jack_I(q_idx) = out_I(q_idx).table.information(3).value;

end

%L? valores m?ximos de Informa??o para cada tipo de entropia
[max_info_plugin_I,max_info_plugin_idx_I] = max(info_plugin_I);
[max_info_tpmc_I,max_info_tpmc_idx_I] = max(info_tpmc_I);
[max_info_jack_I,max_info_jack_idx_I] = max(info_jack_I);

%L? valores de Informa??o para cada tipo de entropia para c?digo Taxa
Hcount_plugin_I = info_plugin_I(1);
Hcount_tpmc_I = info_tpmc_I(1);
Hcount_jack_I = info_jack_I(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CALCULA H-BIAS VIA SHUFFLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Tipo de entropia usada
opts.entropy_estimation_method = {'plugin'};

%Calcula Shuffle para c?digo Timing
rand('state',0);
opts.metric_family = 0;
[Y_T,SHUF_T,opts_used] = metric_shuf(X,opts,S);

SHUF_T = SHUF_T';

%Calcula Shuffle para c?digo Interval
rand('state',0);
opts.metric_family = 1;
[Y_I,SHUF_I,opts_used] = metric_shuf(X,opts,S);

SHUF_I = SHUF_I';

%L? valores de Informa??o para cada c?digo
for q_idx=1:length(opts.shift_cost)
    
    for i=1:S
        
        HShuffle_T(i,q_idx) = SHUF_T(i,q_idx).table.information.value;
        
    end

    for i=1:S
        
        HShuffle_I(i,q_idx) = SHUF_I(i,q_idx).table.information.value;
        
    end

end

%L? a m?dia do Shuffle para cada resolu??o
HBias_T = mean(HShuffle_T,1);
HBias_std_T = std(HShuffle_T,[],1);
HBias_I = mean(HShuffle_I,1);
HBias_std_I = std(HShuffle_I,[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  CALCULA JACK-KNIFE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Seleciona tipo de entropia
opts.entropy_estimation_method = {'plugin'};

%Seleciona tipo de codigo
opts.metric_family = 0;

%Calcula Metric Space Jack Knife
[out_unjk,jk,opts_used] = metric_jack(X,opts);

%L? as Informa??es para cada resolu??o
P_total = size(jk,1);

temp_info_jk = zeros(P_total,length(opts.shift_cost));

for q_idx=1:length(opts.shift_cost)
    
  info_unjk(q_idx)= out_unjk(q_idx).table.information.value;
  
  for p=1:P_total
      
    temp_info_jk(p,q_idx) = jk(p,q_idx).table.information.value;
    
  end
  
end

info_jk_sem = sqrt((P_total-1)*var(temp_info_jk,1,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% SALVA DADOS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info = struct('info_plugin_T',info_plugin_T,'info_tpmc_T',info_tpmc_T,'info_jack_T',info_jack_T,'info_unjk',info_unjk,'info_jk_sem',info_jk_sem,'info_plugin_I',info_plugin_I,'info_tpmc_I',info_tpmc_I,'info_jack_I',info_jack_I);

max_info = struct('max_info_plugin_idx_T',max_info_plugin_idx_T,'max_info_plugin_T',max_info_plugin_T,'max_info_plugin_idx_I',max_info_plugin_idx_I,'max_info_plugin_I',max_info_plugin_I,'Hcount_plugin_info_T',Hcount_plugin_T,'Hcount_plugin_info_I',Hcount_plugin_I,'max_info_tpmc_idx_T',max_info_tpmc_idx_T,'max_info_tpmc_T',max_info_tpmc_T,'max_info_tpmc_idx_I',max_info_tpmc_idx_I,'max_info_tpmc_I',max_info_tpmc_I,'Hcount_tpmc_info_T',Hcount_tpmc_T,'Hcount_tpmc_info_I',Hcount_tpmc_I,'max_info_jack_idx_T',max_info_jack_idx_T,'max_info_jack_T',max_info_jack_T,'max_info_jack_idx_I',max_info_jack_idx_I,'max_info_jack_I',max_info_jack_I,'Hcount_jack_info_T',Hcount_jack_T,'Hcount_jack_info_I',Hcount_jack_I,'HBias_T',HBias_T,'HBias_std_T',HBias_std_T,'HBias_I',HBias_I,'HBias_std_I',HBias_std_I);

metric_analysis = struct('X', X, 'out_T', out_T, 'out_I', out_I, 'info', info, 'max_info', max_info, 'opts', opts, 'Y_I', Y_I, 'Y_T', Y_T, 'SHUF_I', SHUF_I, 'SHUF_T', SHUF_T);

save(nome_do_registro,'metric_analysis');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END - metric2results');

end

