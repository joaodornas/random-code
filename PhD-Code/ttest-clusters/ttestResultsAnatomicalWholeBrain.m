
function ttestResultsAnatomicalWholeBrain

load('CorrelationDensities.mat');

[lFCD_AND, lFCD_ANI, lFCD_APD, lFCD_API, lFCD_PND, lFCD_PNI, lFCD_PPD, lFCD_PPI] = removeNaNs(lFCD_AND, lFCD_ANI, lFCD_APD, lFCD_API, lFCD_PND, lFCD_PNI, lFCD_PPD, lFCD_PPI);

check(lFCD_ANI,lFCD_PNI,'Negative - Increase');
check(lFCD_AND,lFCD_PND,'Negative - Decrease');
check(lFCD_API,lFCD_PPI,'Positive - Increase');
check(lFCD_APD,lFCD_PPD,'Positive - Decrease');

end

function check(Att,Pass,label)

pcriterion = 0.05;

[P,H] = ranksum(Att,Pass);

disp(strcat(label,':',num2str(H),':',num2str(P)));


end

function [lFCD_AND, lFCD_ANI, lFCD_APD, lFCD_API, lFCD_PND, lFCD_PNI, lFCD_PPD, lFCD_PPI] = removeNaNs(lFCD_AND, lFCD_ANI, lFCD_APD, lFCD_API, lFCD_PND, lFCD_PNI, lFCD_PPD, lFCD_PPI);

lFCD_AND(isnan(lFCD_AND)) = [];
lFCD_ANI(isnan(lFCD_ANI)) = [];
lFCD_APD(isnan(lFCD_APD)) = [];
lFCD_API(isnan(lFCD_API)) = [];


lFCD_PND(isnan(lFCD_PND)) = [];
lFCD_PNI(isnan(lFCD_PNI)) = [];
lFCD_PPD(isnan(lFCD_PPD)) = [];
lFCD_PPI(isnan(lFCD_PPI)) = [];

end

