function allTIgroupCC


Protocol(1) = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v5/TI/v5-12-11-01-sitio1-E1-bin_size-30-TI.mat');
Protocol(2) = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v6/TI/v6-12-11-01-sitio1-E1-bin_size-30-TI.mat');
Protocol(3) = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v5/TI/v5-12-11-06-sitio1-E1-bin_size-30-TI.mat');
Protocol(4) = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v6/TI/v6-12-11-06-sitio1-E1-bin_size-30-TI.mat');
Protocol(5) = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v5/TI/v5-12-11-07-sitio1-E1-bin_size-30-TI.mat');
Protocol(6) = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v6/TI/v6-12-11-07-sitio1-E1-bin_size-30-TI.mat');
Protocol(7) = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v1/TI/v1-12-11-08-sitio1-E1-bin_size-30-TI.mat');
Protocol(8) = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v4/TI/v4-12-11-08-sitio1-E1-bin_size-30-TI.mat');
Protocol(9) = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v3/TI/v3-12-11-09-sitio1-E1-bin_size-30-TI.mat');
Protocol(10) = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v6/TI/v6-12-11-09-sitio1-E1-bin_size-30-TI.mat');
Protocol(11) = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v1/TI/v1-12-11-13-sitio1-E1-bin_size-30-TI.mat');
Protocol(12) = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v2/TI/v2-12-11-13-sitio1-E1-bin_size-30-TI.mat');



p = size(Protocol,2);

for i=1:p
   
    condicao(1).InfoPerSecond(i) = Protocol(i).TI.condicao(1).InfoPerSecond;
    condicao(2).InfoPerSecond(i) = Protocol(i).TI.condicao(2).InfoPerSecond;
    
    condicao(1).InfoPerSpike(i) = Protocol(i).TI.condicao(1).InfoPerSpike;
    condicao(2).InfoPerSpike(i) = Protocol(i).TI.condicao(2).InfoPerSpike;
    
    condicao(1).CodingEfficiency(i) = Protocol(i).TI.condicao(1).CodingEfficiency;
    condicao(2).CodingEfficiency(i) = Protocol(i).TI.condicao(2).CodingEfficiency;
    
end


condicao(1).meanInfoPerSecond = mean(condicao(1).InfoPerSecond);
condicao(2).meanInfoPerSecond = mean(condicao(2).InfoPerSecond);


condicao(1).meanInfoPerSpike = mean(condicao(1).InfoPerSpike);
condicao(2).meanInfoPerSpike = mean(condicao(2).InfoPerSpike);


condicao(1).meanCodingEfficiency = mean(condicao(1).CodingEfficiency);
condicao(2).meanCodingEfficiency = mean(condicao(2).CodingEfficiency);

save('TI-group','condicao');

end