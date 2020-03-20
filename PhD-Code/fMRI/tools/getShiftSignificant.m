function shift = getShiftSignificant(rho1,pval1,rho2,pval2)

NComponent = size(rho1,1);

pcriterion = 0.01;

idx1 = find(pval1>pcriterion);
rho1(idx1) = 0;
rho1(rho1==0) = [];

idx2 = find(pval2>pcriterion);
rho2(idx2) = 0;
rho2(rho2==0) = [];

percent1 = length(rho1)/(NComponent*NComponent);
percent2 = length(rho2)/(NComponent*NComponent);

shift = percent1 - percent2;

end
