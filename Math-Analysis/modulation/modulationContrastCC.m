function modulationradiusCC

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

video(1) = 1;
video(2) = 3;
video(3) = 2;
video(4) = 4;
video(5) = 3;
video(6) = 2;
video(7) = 3;
video(8) = 5;
video(9) = 6;
video(10) = 5;
video(11) = 6;
video(12) = 5;
video(13) = 6;
video(14) = 1;
video(15) = 4;
video(16) = 3;
video(17) = 6;
video(18) = 1;
video(19) = 2;
video(20) = 8;
video(21) = 11;
video(22) = 9;
video(23) = 10;
video(24) = 6;
video(25) = 8;
video(26) = 3;
video(27) = 8;
video(28) = 5;
video(29) = 8;
video(30) = 3;
video(31) = 4;
video(32) = 8;
video(33) = 2;
video(34) = 5;
video(35) = 5;
video(36) = 5;
video(37) = 8;
video(38) = 8;
video(39) = 8;
video(40) = 8;
video(41) = 5;
video(42) = 5;
video(43) = 8;
video(44) = 8;
video(45) = 11;
video(46) = 1;
video(47) = 8;

radius(1).CRF = 68;
radius(1).nCRF = 272;

radius(2).CRF = 68;
radius(2).nCRF = 272;

radius(3).CRF = 68;
radius(3).nCRF = 272;

radius(4).CRF = 68;
radius(4).nCRF = 272;

radius(5).CRF = 68;
radius(5).nCRF = 272;

radius(6).CRF = 68;
radius(6).nCRF = 272;

radius(7).CRF = 68;
radius(7).nCRF = 272;

radius(8).CRF = 44;
radius(8).nCRF = 176;

radius(9).CRF = 44;
radius(9).nCRF = 176;

radius(10).CRF = 44;
radius(10).nCRF = 176;

radius(11).CRF = 44;
radius(11).nCRF = 176;

radius(12).CRF = 44;
radius(12).nCRF = 176;

radius(13).CRF = 44;
radius(13).nCRF = 176;

radius(14).CRF = 44;
radius(14).nCRF = 176;

radius(15).CRF = 44;
radius(15).nCRF = 176;

radius(16).CRF = 44;
radius(16).nCRF = 176;

radius(17).CRF = 44;
radius(17).nCRF = 176;

radius(18).CRF = 44;
radius(18).nCRF = 176;

radius(19).CRF = 44;
radius(19).nCRF = 176;

radius(20).CRF = 44;
radius(20).nCRF = 88;

radius(21).CRF = 44;
radius(21).nCRF = 88;

radius(22).CRF = 44;
radius(22).nCRF = 88;

radius(23).CRF = 44;
radius(23).nCRF = 88;

radius(24).CRF = 44;
radius(24).nCRF = 88;

radius(25).CRF = 44;
radius(25).nCRF = 88;

radius(26).CRF = 68;
radius(26).nCRF = 136;

radius(27).CRF = 68;
radius(27).nCRF = 136;

radius(28).CRF = 88;
radius(28).nCRF = 176;

radius(29).CRF = 88;
radius(29).nCRF = 176;

radius(30).CRF = 68;
radius(30).nCRF = 136;

radius(31).CRF = 68;
radius(31).nCRF = 136;

radius(32).CRF = 68;
radius(32).nCRF = 136;

radius(33).CRF = 68; 
radius(33).nCRF = 136;

radius(34).CRF = 68;
radius(34).nCRF = 136;

radius(35).CRF = 68;
radius(35).nCRF = 136;

radius(36).CRF = 68;
radius(36).nCRF = 136;

radius(37).CRF = 68;
radius(37).nCRF = 136;

radius(38).CRF = 68;
radius(38).nCRF = 136;

radius(39).CRF = 68;
radius(39).nCRF = 136;

radius(40).CRF = 68;
radius(40).nCRF = 136;

radius(41).CRF = 68;
radius(41).nCRF = 136;

radius(42).CRF = 68;
radius(42).nCRF = 136;

radius(43).CRF = 68;
radius(43).nCRF = 136;

radius(44).CRF = 44;
radius(44).nCRF = 88;

radius(45).CRF = 44;
radius(45).nCRF = 88;

radius(46).CRF = 44;
radius(46).nCRF = 88;

radius(47).CRF = 44;
radius(47).nCRF = 88;

s = 1;
for i=1:size(Protocolo,2)
     
    meanRateModulationNonCRF(i) = Protocolo(i).rate.suppressionAnalysis.meanRateModulationNonCRF;   

%     if i>=20
%         
%         meanRateModulationExtraNonCRF(s) = Protocolo(i).rate.suppressionAnalysis.meanRateModulationExtraNonCRF;
%         
%         s = s + 1;
%         
%     end

end


%rateModulation = [meanRateModulationNonCRF meanRateModulationExtraNonCRF];

rateModulation = meanRateModulationNonCRF;

for i=1:size(Protocolo,2)
    
   spatialContrast = load(strcat('v',int2str(video(i)),'-radius-',int2str(radius(i).CRF),'-spatialContrast.mat'));
   
   contrastCRF(i) = spatialContrast.spatialContrast.meanMichelson;
   
   spatialContrast = load(strcat('v',int2str(video(i)),'-radius-',int2str(radius(i).nCRF),'-spatialContrast.mat'));
   
   contrastnCRF(i) = spatialContrast.spatialContrast.meanMichelson; 
    
end

idxFacilitacao = find(rateModulation>1)
idxSupressao = find(rateModulation<1)

facilitacao = rateModulation(rateModulation>1);
supressao = rateModulation(rateModulation<1);

contrastFacilitacaoCRF = contrastCRF(idxFacilitacao);
contrastSupressaoCRF = contrastCRF(idxSupressao);

mean(contrastFacilitacaoCRF)
mean(contrastSupressaoCRF)

end

