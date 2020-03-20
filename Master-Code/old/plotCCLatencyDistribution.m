

function plotCCLatencyDistribution


load('CC-latency-distribution.mat');

latencies = [];

latencies = [latencies; latency_distribution(:,1)];

latencies = [latencies; latency_distribution(:,2)];

N = 10;

vector = linspace(min(latencies),max(latencies),N);

[Y, X] = hist(latencies,vector);

f = figure;

bar(X,Y);
title('Center-Surround');
xlabel('Latency (s)');
ylabel('Number of registers');
print(f,'-depsc','CC-Latency-Distribution');


end

