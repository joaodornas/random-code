%fitted coeficients
A=10;
bestSF=0.5;
tuningWidth=1;
skew=0;
lowest_SF=0.25;

%calculate LSFS from fitted curve
sf=lowest_SF;
lowest_SF_resp=A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))));
sf=bestSF;
best_SF_resp=A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))));
LSFS=lowest_SF_resp/best_SF_resp;
clear sf

%calculate LSFV
syms sf 
a=int((A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))))).*(log(sf)/log(16)-log(bestSF)/log(16)).^2,sf,bestSF/16,bestSF);
b=int((A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))))),sf,bestSF/16,bestSF);
LSFV=double(a)/double(b)

%plot tuning curve from model function

sf=0.1:0.1:64;
response=A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))));
plot(log(sf),response);
clear sf