function histModulationCC


Protocolo(1).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v1/suppression/v1-12-10-22-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(2).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-22/sitio1/E2/v3/suppression/v3-12-10-22-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(3).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v2/suppression/v2-12-10-23-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(4).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-23/sitio1/E2/v4/suppression/v4-12-10-23-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(5).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E1/v3/suppression/v3-12-10-24-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(6).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-24/sitio1/E2/v2/suppression/v2-12-10-24-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(7).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-26/sitio2/E1/v3/suppression/v3-12-10-26-sitio2-E1-bin_size-30-suppression.mat');
% Protocolo(8).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v5/suppression/v5-12-10-31-sitio1-E3-bin_size-30-suppression.mat');
% Protocolo(9).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-10-31/sitio1/E3/v6/suppression/v6-12-10-31-sitio1-E3-bin_size-30-suppression.mat');
Protocolo(8).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v5/suppression/v5-12-11-01-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(9).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-01/sitio1/E1/v6/suppression/v6-12-11-01-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(10).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v5/suppression/v5-12-11-06-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(11).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-06/sitio1/E1/v6/suppression/v6-12-11-06-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(12).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v5/suppression/v5-12-11-07-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(13).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-07/sitio1/E1/v6/suppression/v6-12-11-07-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(14).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v1/suppression/v1-12-11-08-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(15).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-08/sitio1/E1/v4/suppression/v4-12-11-08-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(16).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v3/suppression/v3-12-11-09-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(17).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-09/sitio1/E1/v6/suppression/v6-12-11-09-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(18).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v1/suppression/v1-12-11-13-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(19).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-11-13/sitio1/E1/v2/suppression/v2-12-11-13-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(20).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v8/suppression/v8-12-12-03-sitio1-E3-bin_size-30-suppression.mat');
Protocolo(21).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-03/sitio1/E3/v11/suppression/v11-12-12-03-sitio1-E3-bin_size-30-suppression.mat');
Protocolo(22).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v9/suppression/v9-12-12-04-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(23).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-04/sitio1/E1/v10/suppression/v10-12-12-04-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(24).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-05/sitio1/E1/v6/suppression/v6-12-12-05-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(25).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-05/sitio1/E1/v8/suppression/v8-12-12-05-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(26).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-07/sitio1/E1/v3/suppression/v3-12-12-07-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(27).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-07/sitio1/E1/v8/suppression/v8-12-12-07-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(28).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-10/sitio1/E2/v5/suppression/v5-12-12-10-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(29).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-10/sitio1/E2/v8/suppression/v8-12-12-10-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(30).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-12/sitio2/E1/v3/suppression/v3-12-12-12-sitio2-E1-bin_size-30-suppression.mat');
Protocolo(31).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-12/sitio2/E1/v4/suppression/v4-12-12-12-sitio2-E1-bin_size-30-suppression.mat');
Protocolo(32).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-13/sitio1/E1/v8/suppression/v8-12-12-13-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(33).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-13/sitio1/E1/v2/suppression/v2-12-12-13-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(34).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v511/suppression/v511-12-12-17-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(35).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v512/suppression/v512-12-12-17-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(36).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v521/suppression/v521-12-12-17-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(37).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v811/suppression/v811-12-12-17-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(38).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v812/suppression/v812-12-12-17-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(39).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v821/suppression/v821-12-12-17-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(40).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E1/v822/suppression/v822-12-12-17-sitio1-E1-bin_size-30-suppression.mat');
Protocolo(41).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v531/suppression/v531-12-12-17-sitio1-E3-bin_size-30-suppression.mat');
Protocolo(42).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v532/suppression/v532-12-12-17-sitio1-E3-bin_size-30-suppression.mat');
Protocolo(43).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-17/sitio1/E3/v831/suppression/v831-12-12-17-sitio1-E3-bin_size-30-suppression.mat');
Protocolo(44).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-18/sitio1/E2/v8/suppression/v8-12-12-18-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(45).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-18/sitio1/E2/v11/suppression/v11-12-12-18-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(46).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-19/sitio1/E2/v1/suppression/v1-12-12-19-sitio1-E2-bin_size-30-suppression.mat');
Protocolo(47).rate = load('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/12-12-19/sitio1/E2/v8/suppression/v8-12-12-19-sitio1-E2-bin_size-30-suppression.mat');

s = 1;
r = 1;
for i=1:size(Protocolo,2)
    
   modulation(i,:) = Protocolo(i).rate.suppressionAnalysis.NonCRFxCRF(:);
   
   meanRateModulationNonCRF(i) = Protocolo(i).rate.suppressionAnalysis.meanRateModulationNonCRF;
   
   if i>=20
      
      modulationExtra(s,:) = Protocolo(i).rate.suppressionAnalysis.ExtraNonCRFxCRF(:); 
      
      meanRateModulationExtraNonCRF(s) = Protocolo(i).rate.suppressionAnalysis.meanRateModulationExtraNonCRF;
       
      s = s + 1;
      
   end
    
   
   if i<=10
   
        meanRateModulationAlong(i) = Protocolo(i).rate.suppressionAnalysis.condition(1).meanRate / Protocolo(i).rate.suppressionAnalysis.condition(2).meanRate;
   
   else
       
       meanRateModulationAlong(i) = Protocolo(i).rate.suppressionAnalysis.condition(2).meanRate / Protocolo(i).rate.suppressionAnalysis.condition(1).meanRate;
       
   end
   
   if i>=20
       
       meanRateModulationAlongExtra(r) = Protocolo(i).rate.suppressionAnalysis.condition(3).meanRate / Protocolo(i).rate.suppressionAnalysis.condition(1).meanRate;
      
       r = r + 1;
       
   end
       
end

for i=1:size(Protocolo,2)
    
   
    meanMod(i) = mean(modulation(i,:));
    
end


for i=1:(s-1)
    
   meanModExtra(i) = mean(modulationExtra(i,:));    
    
end

Mod = [meanMod meanModExtra];
RateModulation = [meanRateModulationNonCRF meanRateModulationExtraNonCRF];
RateModulationAlong = [meanRateModulationAlong meanRateModulationAlongExtra];

figure;
bar(Mod,'b');

hold on;
bar(RateModulation,'g');

hold on;
p1 = [1 1];
p2 = [0 size(RateModulation,2)];
plot([p2(1) p2(2)],[p1(1) p1(1)],'Color','r','LineWidth',2);

figure;
bar(RateModulationAlong);

positivo = [];
negativo = [];
s = 1;
r = 1;
for i=1:length(Mod)
   
    if Mod(i) > 0
        
        positivo(s) = Mod(i);
        s = s + 1;
        
    elseif Mod(i) < 0
        
        negativo(r) = Mod(i);
        r = r + 1;
        
    end
    
end

mean(Mod)

nPositivo = s-1
nNegativo = r-1

nM1 = length(RateModulationAlong(RateModulationAlong > 1))
nm1 = length(RateModulationAlong(RateModulationAlong < 1))

end

