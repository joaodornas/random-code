function rasterPlotForBack

registro = importdata('memoryBackwardProtocols.txt');

nConditions = 3;

for r=1:length(registro)
    
    protocols = char(registro(r));
    
    rasterAll(protocols,nConditions);

end

end

