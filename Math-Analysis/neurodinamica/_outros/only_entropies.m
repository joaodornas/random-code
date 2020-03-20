
function [X, out_T, opts] = only_entropies(dataset, start_time, end_time, entropy)

echo on;

opts.clustering_exponent = -2;
opts.unoccupied_bins_strategy = 0;
opts.parallel = 1;
opts.possible_words = 'unique';
opts.shift_cost = [0 2.^(-4:9)];

opts.metric_family = 0;

opts.start_time = start_time;
opts.end_time = end_time;

switch entropy
    
    case 0
        
        opts.entropy_estimation_method = {'plugin','tpmc','jack','ma','bub','chaoshen','ww','nsb'};

    case 1
        
        opts.entropy_estimation_method = {'plugin'};

    case 2
        
        opts.entropy_estimation_method = {'tpmc'};

    case 3
        
        opts.entropy_estimation_method = {'jack'};

    case 4
        
        opts.entropy_estimation_method = {'ma'};

    case 5
        
        opts.entropy_estimation_method = {'bub'};
    
    case 6
        
        opts.entropy_estimation_method = {'chaoshen'};

    case 7
        
        opts.entropy_estimation_method = {'ww'};

    case 8
        
        opts.entropy_estimation_method = {'nsb'};

end

X = staread(strrep(dataset,'/',filesep));

clear out out_unshuf shuf out_unjk jk;
clear info_plugin info_tpmc info_jack info_unshuf info_unjk;
clear temp_info_shuf temp_info_jk;

[out_T,opts_used] = metric(X,opts);

for q_idx=1:length(opts.shift_cost)

    if entropy == 0
        
        info_plugin(q_idx) = out_T(q_idx).table.information(1).value;
        info_tpmc(q_idx) = out_T(q_idx).table.information(2).value;
        info_jack(q_idx) = out_T(q_idx).table.information(3).value;
        info_ma(q_idx) = out_T(q_idx).table.information(4).value;
        info_bub(q_idx) = out_T(q_idx).table.information(5).value;
        info_chaoshen(q_idx) = out_T(q_idx).table.information(6).value;
        info_ww(q_idx) = out_T(q_idx).table.information(7).value;
        info_nsb(q_idx) = out_T(q_idx).table.information(8).value;
        
    else
        
        info_whatever(q_idx) = out_T(q_idx).table.information(1).value;
        
    end

end

figure;
set(gcf,'name',['Metric - Timing - ' dataset]); 

if entropy == 0
    
    plot(1:length(opts.shift_cost),info_plugin,'b');
    hold on;
    plot(1:length(opts.shift_cost),info_tpmc,'b--');
    plot(1:length(opts.shift_cost),info_jack,'b-.');
    plot(1:length(opts.shift_cost),info_ma,'b:');
    plot(1:length(opts.shift_cost),info_bub,'r');
    plot(1:length(opts.shift_cost),info_chaoshen,'r--');
    plot(1:length(opts.shift_cost),info_ww,'r-.');
    plot(1:length(opts.shift_cost),info_nsb,'r:');
    hold off;
    set(gca,'xtick',1:length(opts.shift_cost));
    set(gca,'xticklabel',opts.shift_cost);
    set(gca,'xlim',[1 length(opts.shift_cost)]);
    set(gca,'ylim',[-0.5 2.5]);
    xlabel('Temporal precision (1/sec)');
    ylabel('Information (bits)');
    legend('No correction','TPMC','Jackknife','Ma bound','Best upper bound','Chao-Shen','Wolpert-Wolf Bayesian with a Dirichlet prior','Nemenman-Shafee-Bialek');

else
    
    switch entropy
        
        case 1
            
            entropy_name = 'plugin';
            
        case 2
            
            entropy_name = 'TPMC';
            
        case 3
            
            entropy_name = 'Jackknife';
            
        case 4
            
            entropy_name = 'Ma bound';
            
        case 5
            
            entropy_name = 'Best upper bound';
            
        case 6
            
            entropy_name = 'Chao-Shen';
            
        case 7
            
            entropy_name = 'Wolpert-Wolf Bayesian with a Dirichlet prior';
            
        case 8
            
            entropy_name = 'Nemenman-Shafee-Bialek';
            
    end
    
    plot(1:length(opts.shift_cost),info_whatever,'b');
    set(gca,'xtick',1:length(opts.shift_cost));
    set(gca,'xticklabel',opts.shift_cost);
    set(gca,'xlim',[1 length(opts.shift_cost)]);
    set(gca,'ylim',[-0.5 2.5]);
    xlabel('Temporal precision (1/sec)');
    ylabel('Information (bits)');
    legend(entropy_name);
    
end

end

