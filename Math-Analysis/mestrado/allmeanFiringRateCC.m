function allmeanFiringRateCC


Protocolo(1).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v1/suppression/v1-12-10-22-sitio1-E2-bin_size-30-suppression');
Protocolo(2).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v3/suppression/v3-12-10-22-sitio1-E2-bin_size-30-suppression');
Protocolo(3).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v2/suppression/v2-12-10-23-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(4).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v4/suppression/v4-12-10-23-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(5).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E1/v3/suppression/v3-12-10-24-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(6).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E2/v2/suppression/v2-12-10-24-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(7).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-26/sitio2/E1/v3/suppression/v3-12-10-26-sitio2-E1-bin_size-30-suppression.mat');
Protocolo(8).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v5/suppression/v5-12-10-31-sitio1-E3-bin_size-30-suppression.mat');
Protocolo(9).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v6/suppression/v6-12-10-31-sitio1-E3-bin_size-30-suppression.mat');
Protocolo(10).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v5/suppression/v5-12-11-01-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(11).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v6/suppression/v6-12-11-01-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(12).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v5/suppression/v5-12-11-06-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(13).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v6/suppression/v6-12-11-06-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(14).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v5/suppression/v5-12-11-07-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(15).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v6/suppression/v6-12-11-07-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(16).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v1/suppression/v1-12-11-08-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(17).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v4/suppression/v4-12-11-08-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(18).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v3/suppression/v3-12-11-09-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(19).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v6/suppression/v6-12-11-09-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(20).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v2/suppression/v2-12-11-13-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(21).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v1/suppression/v1-12-11-13-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(22).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v8/suppression/v8-12-12-03-sitio1-E3-bin_size-30-suppression');
Protocolo(23).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v11/suppression/v11-12-12-03-sitio1-E3-bin_size-30-suppression');
Protocolo(24).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v9/suppression/v9-12-12-04-sitio1-E1-bin_size-30-suppression');
Protocolo(25).suppression = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v10/suppression/v10-12-12-04-sitio1-E1-bin_size-30-suppression');

for i=1:size(Protocolo,2)
    
   meanFiringRateNonCRF(i) = Protocolo(i).suppression.suppressionAnalysis.conditions(1).meanFiringRate;
   meanFiringRateCRF(i) = Protocolo(i).suppression.suppressionAnalysis.conditions(2).meanFiringRate;
   
   if i >= 10
       
          meanFiringRateCRF(i) = Protocolo(i).suppression.suppressionAnalysis.conditions(1).meanFiringRate;
          meanFiringRateNonCRF(i) = Protocolo(i).suppression.suppressionAnalysis.conditions(2).meanFiringRate;

   end
   
   if i > 21
       
       meanFiringRateExtraNonCRF(i) = Protocolo(i).suppression.suppressionAnalysis.conditions(3).meanFiringRate;
       
   end
   
end



meanFiringRateCRF = meanFiringRateCRF.';
meanFiringRateNonCRF = meanFiringRateNonCRF.';
meanFiringRateExtraNonCRF = meanFiringRateExtraNonCRF.';

meanFiringRateMatrix = horzcat(meanFiringRateCRF,meanFiringRateNonCRF);

figure;
bar(meanFiringRateMatrix);

meanFiringRateExtraMatrix = horzcat(meanFiringRateCRF,meanFiringRateExtraNonCRF);

j = 0;
for i=10:size(Protocolo,2)
    
    j = j + 1;
    
    increase(j) = meanFiringRateNonCRF(i) / meanFiringRateCRF(i) * 100 - 100;
    
    meanFiringRateCRF(j) = meanFiringRateCRF(i);
    meanFiringRateNonCRF(j) = meanFiringRateNonCRF(i);
    
end

meanIncrease = mean(increase);
stdIncrease = std(increase);

meanFRCRF = mean(meanFiringRateCRF)*100;
meanFRNonCRF = mean(meanFiringRateNonCRF)*100;

h = lillietest(meanFiringRateCRF)
h = lillietest(meanFiringRateNonCRF)

[H p] = ttest(meanFiringRateCRF,meanFiringRateNonCRF)

text(1,24,strcat('Increase - ',' Mean:',int2str(meanIncrease),'  /  ','Std:',int2str(stdIncrease)));

text(1,23,strcat('MeanFRCRF:',int2str(meanFRCRF),'MeanFRNonCRF:',int2str(meanFRNonCRF)));

figure;

bar(meanFiringRateExtraMatrix);

end