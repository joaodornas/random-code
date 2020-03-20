function [LSFV,LSFS]=calculate_lsfv(fittedParameters,lowest_SF,plotTuningCurve)

%calculates low spatial frequency variance (lsfv) and low spatial frequency
%suppression based on method described in Xing et al (2002) J Physiol
%
%[LSFV,LSFS]=calculate_lsfv(fittedParameters,lowest_SF,plotTuningCurve)
%where fittedParamaters is a cftool object containing fitted model parameters
%(default matlab name fittedmodel1) from Priebe et al (2006)'s 1D fitting 
%equation, lowest_SF is lowest tested SF and plotTuningCurve is a string 
%argument which should be set to 'plot' in case one wants to visualize the tuning curve.

warning off symbolic:sym:int:warnmsg1

%fitted coeficients
A=fittedParameters.A;
bestSF=fittedParameters.ps;
tuningWidth=fittedParameters.tuningWidth;
skew=fittedParameters.skew;

%calculate LSFS from fitted curve
sf=lowest_SF;
lowest_SF_resp=A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))));
sf=bestSF;
best_SF_resp=A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))));
LSFS=lowest_SF_resp/best_SF_resp;
clear sf

%calculate LSFV
syms sf
a=int((A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))))).*(log2(sf)/log2(16)-log2(bestSF)/log2(16)).^2,bestSF/16,bestSF);
b=int((A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))))),bestSF/16,bestSF);
LSFV=double(a)/double(b);

%plot tuning curve from model function
switch plotTuningCurve
    case 'plot'
        sfs=0.1:0.1:64;
        response=A*((exp((-((log2(sfs)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sfs)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))));
        plot(log(sfs),response);
end

end

