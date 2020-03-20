function allSparsenessGroupCC

Protocolo(1).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v1/sparseness/v1-12-10-22-sitio1-E2-bin_size-30-sparseness.mat');
Protocolo(2).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v3/sparseness/v3-12-10-22-sitio1-E2-bin_size-30-sparseness.mat');
Protocolo(3).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v2/sparseness/v2-12-10-23-sitio1-E2-bin_size-30-sparseness.mat');
Protocolo(4).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v4/sparseness/v4-12-10-23-sitio1-E2-bin_size-30-sparseness.mat');
Protocolo(5).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E1/v3/sparseness/v3-12-10-24-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(6).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E2/v2/sparseness/v2-12-10-24-sitio1-E2-bin_size-30-sparseness.mat');
Protocolo(7).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-26/sitio2/E1/v3/sparseness/v3-12-10-26-sitio2-E1-bin_size-30-sparseness.mat');
Protocolo(8).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v5/sparseness/v5-12-10-31-sitio1-E3-bin_size-30-sparseness.mat');
Protocolo(9).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v6/sparseness/v6-12-10-31-sitio1-E3-bin_size-30-sparseness.mat');
Protocolo(10).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v5/sparseness/v5-12-11-01-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(11).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v6/sparseness/v6-12-11-01-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(12).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v5/sparseness/v5-12-11-06-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(13).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v6/sparseness/v6-12-11-06-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(14).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v5/sparseness/v5-12-11-07-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(15).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v6/sparseness/v6-12-11-07-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(16).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v1/sparseness/v1-12-11-08-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(17).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v4/sparseness/v4-12-11-08-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(18).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v3/sparseness/v3-12-11-09-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(19).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v6/sparseness/v6-12-11-09-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(20).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v2/sparseness/v2-12-11-13-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(21).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v1/sparseness/v1-12-11-13-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(22).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v8/sparseness/v8-12-12-03-sitio1-E3-bin_size-30-sparseness.mat');
Protocolo(23).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v11/sparseness/v11-12-12-03-sitio1-E3-bin_size-30-sparseness.mat');
Protocolo(24).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v9/sparseness/v9-12-12-04-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(25).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v10/sparseness/v10-12-12-04-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(26).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-05/sitio1/E1/v6/sparseness/v6-12-12-05-sitio1-E1-bin_size-30-sparseness.mat');
Protocolo(27).sparseness = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-05/sitio1/E1/v8/sparseness/v8-12-12-05-sitio1-E1-bin_size-30-sparseness.mat');

for i=1:size(Protocolo,2)
    
    sparsenessNonCRF(i) = Protocolo(i).sparseness.sparseData(1,1).gallant2002Sparseness;
    sparsenessCRF(i) = Protocolo(i).sparseness.sparseData(1,2).gallant2002Sparseness;
    
    if i>=10
        
        sparsenessCRF(i) = Protocolo(i).sparseness.sparseData(1,1).gallant2002Sparseness;
        sparsenessNonCRF(i) = Protocolo(i).sparseness.sparseData(1,2).gallant2002Sparseness;
        
    end
    
end

sparsenessCRF = sparsenessCRF';
sparsenessNonCRF = sparsenessNonCRF';

sparsenessMatrix = horzcat(sparsenessCRF, sparsenessNonCRF);

j = 0;
for i=10:size(Protocolo,2)
    
    j = j + 1;
    
    increase(j) = sparsenessNonCRF(i) / sparsenessCRF(i) * 100 - 100;
    shift(j) = ( sparsenessNonCRF(i) - sparsenessCRF(i) ) / ( 1 - sparsenessCRF(i) ) * 100;
    
    sparseCRF(j) = sparsenessCRF(i);
    sparseNonCRF(j) = sparsenessNonCRF(i);
    
end

meanIncrease = mean(increase);
stdIncrease = std(increase);

meanShift = mean(shift);
stdShift = std(shift);

meanSparsenessCRF = mean(sparseCRF)*100;
meanSparsenessNonCRF = mean(sparseNonCRF)*100;

h = lillietest(sparseCRF)
h = lillietest(sparseNonCRF)

[H p] = ttest(sparseCRF,sparseNonCRF)

figure;
bar(sparsenessMatrix);

text(1,0.65,strcat('Increase - ',' Mean:',int2str(meanIncrease),'  /  ','Std:',int2str(stdIncrease)));

text(10,0.65,strcat('SHIFT - ',' Mean:',int2str(meanShift),'  /  ','Std:',int2str(stdShift)));

text(1,0.6,strcat('MeanCRF:',int2str(meanSparsenessCRF),'MeanNonCRF:',int2str(meanSparsenessNonCRF)));

j = 0;
for i=22:size(Protocolo,2)
   
    j = j + 1;
    
    sparseness3CCRF(j) = Protocolo(i).sparseness.sparseData(1,1).gallant2002Sparseness;
    sparseness3CNonCRF(j) = Protocolo(i).sparseness.sparseData(1,2).gallant2002Sparseness;
    sparseness3CExtraNonCRF(j) = Protocolo(i).sparseness.sparseData(1,3).gallant2002Sparseness;   
    
end


sparseness3CCRF = sparseness3CCRF.';
sparseness3CNonCRF = sparseness3CNonCRF.';
sparseness3CExtraNonCRF = sparseness3CExtraNonCRF.';

sparseness3CMatrix = horzcat(sparseness3CCRF,sparseness3CNonCRF,sparseness3CExtraNonCRF);

figure;
bar(sparseness3CMatrix);

end
