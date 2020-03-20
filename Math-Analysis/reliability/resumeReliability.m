function resumeReliability

registro = importdata('memoryBackwardProtocols.txt');

for r=1:length(registro)
    
    protocol = registro{r};
    
    if r < 27
        
        protocol = protocol(2:end-7);
        
    else
        
        protocol = protocol(2:end-9);
        
    end
    
    getDataFor = load(strcat(protocol,'-reliabilityOptimal-Forward.mat'));
    getDataBack = load(strcat(protocol,'-reliabilityOptimal-Backward.mat'));
    
    results(r).name = getDataFor.data.name;
    
    results(r).max_for_real = max(getDataFor.data.max_for_real);
    results(r).max_resolution_for = max(getDataFor.data.max_resolution_for);

    max_for_real(r) = max(getDataFor.data.max_for_real);
    max_resolution_for(r) = max(getDataFor.data.max_resolution_for);
    
    results(r).max_back_real = max(getDataBack.data.max_back_real);
    results(r).max_resolution_back = max(getDataBack.data.max_resolution_back);
    
end

resultados = struct('results',results);

save('/Volumes/Data/DATA/Forward-Backward/reliabilityOptimal/resultados','resultados');

f = figure;
bar(max_for_real)
xlim([0 length(max_for_real)]);
ylim([0 max(max_for_real)]);

print(f,'-depsc',strcat('/Volumes/Data/DATA/Forward-Backward/reliabilityOptimal/','reliability-distribution-max-realiability-forward.eps'));

g = figure;
bar(max_resolution_for)
xlim([0 length(max_resolution_for)]);
ylim([0 max(max_resolution_for)]);

print(g,'-depsc',strcat('/Volumes/Data/DATA/Forward-Backward/reliabilityOptimal/','resolution-distribution-max-resolution-forward.eps'));




end

