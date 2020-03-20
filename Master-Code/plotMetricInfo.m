
f = figure;

load('C:\Users\João V. Dornas\Hard Quale .com\DornasLab.org - Documents\_CODE\TP-RANDOM\Master-Code\MetricSpace\500-9500\nsp008a02_1a-v3\full-movie\full-movie-original\nsp008a02_1a-v3-full-movie-original-WO-latency.mat');

for i=1:15
    info(i) = metric_analysis.Y_T(i).table.information(1).value;
end

for j=1:15
    for k=1:60
        info_std(i,k) = metric_analysis.SHUF_T(k,i).table.information.value;
    end
end

istd = std(info_std,0,2);

errorbar(1:length(opts.shift_cost),info,istd,'b');

hold on

load('C:\Users\João V. Dornas\Hard Quale .com\DornasLab.org - Documents\_CODE\TP-RANDOM\Master-Code\MetricSpace\500-9500\nsp008a02_1a-v3\full-movie\full-movie-invertido\nsp008a02_1a-v3-bin-size-1-full-movie-invertido-WO-latency.mat');

for i=1:15
    info(i) = metric_analysis.Y_T(i).table.information(1).value;
end

for j=1:15
    for k=1:60
        info_std(i,k) = metric_analysis.SHUF_T(k,i).table.information.value;
    end
end

istd = std(info_std,0,2);

errorbar(1:length(opts.shift_cost),info,istd,'r');

resolution = 1 ./ opts.shift_cost ;

set(gca,'xtick',1:length(opts.shift_cost));
set(gca,'xticklabel',resolution);
set(gca,'xlim',[1 length(opts.shift_cost)]);
set(gca,'ylim',[0 2]);
xlabel('Resolution (sec)');
ylabel('Information (bits)');
legend({'Before Inverting','After Inverting'});


