function metricjack(dataset, start_time, end_time)

ct = cputime;

disp('READING DATA');
% READ DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%
X = staread(strrep(dataset,'/',filesep));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('SETTING OPTIONS');
% OPTIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.entropy_estimation_method = {'plugin'};
%opts.variance_estimation_method = {'jack'};

%opts.unoccupied_bins_strategy = -1; % Ignore unoccupied bins
opts.unoccupied_bins_strategy = 0; % Use an unoccupied bin only if its row and column are occupied
%opts.unoccupied_bins_strategy = 1; % Use all bins

opts.parallel = 1;
opts.possible_words = 'unique';

opts.start_time = start_time / 1000;
opts.end_time = end_time / 1000;
opts.shift_cost = [0 2.^(-4:9)];
%opts.label_cost = [0 1 2];
opts.clustering_exponent = -2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('METRIC TIMING JACKKNIFE');
% METRIC TIMING   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.metric_family = 0;

%%% leave-one-out Jackknife 
disp('leave-one-out Jackknife');

[out_unjk,jk,opts_used] = metric_jack(X,opts);

P_total = size(jk,1);

temp_info_jk = zeros(P_total,length(opts.shift_cost));

for q_idx=1:length(opts.shift_cost)
    
  info_unjk(q_idx)= out_unjk(q_idx).table.information.value;
  
  for p=1:P_total
      
    temp_info_jk(p,q_idx) = jk(p,q_idx).table.information.value;
    
  end
  
end

info_jk_sem = sqrt((P_total-1)*var(temp_info_jk,1,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING INFO TIMING-JACKKNIFE');
% INFO + SHUFFLE INFO PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f7 = figure;
errorbar(1:length(opts.shift_cost),info_unjk,2*info_jk_sem);
set(gca,'xtick',1:length(opts.shift_cost));
set(gca,'xticklabel',opts.shift_cost);
set(gca,'xlim',[1 length(opts.shift_cost)]);
set(gca,'ylim',[-0.5 2.5]);
xlabel('Temporal precision (1/sec)');
ylabel('Information (bits)');
legend('Original data (\pm 2 SE via Jackknife)');

scalefig(gcf,1.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cte = cputime - ct;

disp(strcat('CPUTIME: ',num2str(cte)));


end

