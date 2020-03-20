function metricshuffle(dataset, start_time, end_time, nRepetitions)


ct = cputime;

disp('READING DATA');
% READ DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%
X = staread(strrep(dataset,'/',filesep));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('SETTING OPTIONS');
% OPTIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.entropy_estimation_method = {'plugin','tpmc','jack'};
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

S = nRepetitions;

disp('CALCULATE H-BIAS VIA SHUFFLE');
% CALCULATE H-BIAS VIA SHUFFLE %%%%%%%%%%%%%%%

disp('...computing metric space timing shuffle');

opts.entropy_estimation_method = {'plugin'};

rand('state',0);

opts.metric_family = 0;

[Y_T,SHUF_T,opts_used] = metric_shuf(X,opts,S);

SHUF_T = SHUF_T';

disp('...computing metric space interval shuffle');

rand('state',0);

opts.metric_family = 1;

[Y_I,SHUF_I,opts_used] = metric_shuf(X,opts,S);

SHUF_I = SHUF_I';

for q_idx=1:length(opts.shift_cost)
    
    for i=1:S
        
        HShuffle_T(i,q_idx) = SHUF_T(i,q_idx).table.information.value;
        
    end

    for i=1:S
        
        HShuffle_I(i,q_idx) = SHUF_I(i,q_idx).table.information.value;
        
    end

end

HBias_T = mean(HShuffle_T,1);
HBias_std_T = std(HShuffle_T,[],1);
HBias_I = mean(HShuffle_I,1);
HBias_std_I = std(HShuffle_I,[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PLOTTING INFO TIMING-SHUFFLE');
% INFO + SHUFFLE INFO PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f7 = figure;
hold on;
errorbar(1:length(opts.shift_cost),HBias_T,2*HBias_std_T,'r');
hold off;
set(gca,'xtick',1:length(opts.shift_cost));
set(gca,'xticklabel',opts.shift_cost);
set(gca,'xlim',[1 length(opts.shift_cost)]);
set(gca,'ylim',[-0.5 2.5]);
xlabel('Temporal precision (1/sec)');
ylabel('Information (bits)');
legend('Shuffled data (\pm 2 SD)');

scalefig(gcf,1.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cte = cputime - ct;

disp(strcat('CPUTIME: ',num2str(cte)));

end
