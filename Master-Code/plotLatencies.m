%% plot Latency distribution

files = dir('Latency\*.mat');

nFiles = length(files);

for iFile=1:nFiles
    
   load(strcat('Latency\',files(iFile).name));
   
   latencies(iFile) = latency_time;
   
   clear latency_time
    
end


