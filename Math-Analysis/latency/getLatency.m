function latency = getLatency(registro,condition)


database = load(strcat(char(registro),'-latencies-for-conditions.mat'));

latency = database.results.latency(condition);


end

