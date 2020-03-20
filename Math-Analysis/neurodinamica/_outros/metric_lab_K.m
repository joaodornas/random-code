
function [X, out, opts] = metric_lab_K(dataset, start_time, end_time)

echo on;

opts.clustering_exponent = -2;
opts.unoccupied_bins_strategy = 0;
opts.parallel = 1;
opts.possible_words = 'unique';
opts.shift_cost = [0 2.^(-4:9)];

opts.metric_family = 0;

opts.start_time = start_time;
opts.end_time = end_time;

opts.entropy_estimation_method = {'plugin','tpmc','jack','ma','bub','chaoshen','ww','nsb'};

X = staread(strrep(dataset,'/',filesep));

clear out out_unshuf shuf out_unjk jk;
clear info_plugin info_tpmc info_jack info_unshuf info_unjk;
clear temp_info_shuf temp_info_jk;

[out,opts_used] = metric(X,opts);

for q_idx=1:length(opts.shift_cost)

    info_plugin(q_idx) = out(q_idx).table.information(1).value;
    info_tpmc(q_idx) = out(q_idx).table.information(2).value;
    info_jack(q_idx) = out(q_idx).table.information(3).value;
    info_ma(q_idx) = out(q_idx).table.information(4).value;
    info_bub(q_idx) = out(q_idx).table.information(5).value;
    info_chaoshen(q_idx) = out(q_idx).table.information(6).value;
    info_ww(q_idx) = out(q_idx).table.information(7).value;
    info_nsb(q_idx) = out(q_idx).table.information(8).value;

end

figure;
set(gcf,'name',['Metric - Timing - ' dataset]); 
    
subplot(131);
[max_info,max_info_idx]=max(info_plugin);
imagesc(out(max_info_idx).d);
xlabel('Spike train index');
ylabel('Spike train index');
title('Distance matrix at maximum information');

subplot(132);
[max_info,max_info_idx]=max(info_plugin);
imagesc(out(max_info_idx).cm);
xlabel('Spike train index');
ylabel('Spike train index');
title('Confusion matrix from clustering of distances');

subplot(133);

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

end

