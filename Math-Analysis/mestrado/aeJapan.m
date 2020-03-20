function aeJapan

  
    protocol(1).japan = load('v3-12-08-31-sitio1-E1-psth-kernel.mat');
    protocol(1).cell_label = 'v3-12-08-31-sitio1-E1';
    
    protocol(2).japan = load('v311-12-12-17-sitio1-E1-psth-kernel.mat');
    protocol(2).cell_label = 'v311-12-12-17-sitio1-E1';
    
    protocol(3).japan = load('v3-12-08-29-sitio1-E1-psth-kernel.mat');
    protocol(3).cell_label = 'v3-12-08-29-sitio1-E1';
    
    protocol(4).japan = load('v2-12-09-04-sitio1-E1-psth-kernel.mat');
    protocol(4).cell_label = 'v2-12-09-04-sitio1-E1';


for i=1:length(protocol)

       
       density = protocol(i).japan.histData(1).condition(3).kernel.density./max(protocol(i).japan.histData(1).condition(3).kernel.density);
       maxY = max(protocol(i).japan.histData(1).condition(3).psth.binHeight);
       density = density.*maxY;
       
       w(i) = figure;
       plot(protocol(i).japan.histData(1).condition(3).kernel.timePoints,density,'g');
       
       ylim([0 maxY]);
       ylabel('Taxa de Disparo por Segundo');
       xlabel('Tempo (s)');
    
       print(w(i),'-depsc',strcat(protocol(i).cell_label,'-condition-3-ae'));
       
end



end

