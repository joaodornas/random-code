responseMatrix=responseMatrix';
meanResponse=mean(responseMatrix);
sst=sum((responseMatrix-meanResponse).^2);
n=size(responseMatrix,1);
normResponseMatrix=responseMatrix/max(responseMatrix);
normMeanResponse=mean(normResponseMatrix);
sstNorm=sum((normResponseMatrix-normMeanResponse).^2);

%PriebeFit
initialValuesPriebe=[scaling;prefSf;sigmaSf;prefTf;sigmaTf;skew;dependence];
[priebeParameters,priebeResiduals,priebeJacobian]=nlinfit(conditionMatrix,responseMatrix,@priebefit,initialValuesPriebe,statset('MaxIter',1000,'Display','iter'));
priebeParametersCI=nlparci(priebeParameters,priebeResiduals,priebeJacobian);
[priebePred,deltaPriebeCI] = nlpredci(@priebefit,conditionMatrix,priebeParameters,priebeResiduals,priebeJacobian);
nlintool(conditionMatrix,responseMatrix,@priebefit,initialValuesPriebe);
priebePredBounds=[responseMatrix priebePred deltaPriebeCI];
kPriebe=size(priebeParameters,1);
ssePriebe=sum(priebeResiduals.^2);
rsquaredPriebe=1-(ssePriebe/sst);
adj_rsquaredPriebe=1-((ssePriebe/(n-kPriebe-1))/(sst/(n-1)));
priebeAIC=n*log(ssePriebe/n)+2*kPriebe;
priebeAICc=priebeAIC+((2*(kPriebe+2)*(kPriebe+3))/(n-kPriebe-3));
Fpriebe=(rsquaredPriebe/(kPriebe-1))/((1-rsquaredPriebe)/(n-kPriebe));
p_valuePriebe=1-fcdf(Fpriebe,(kPriebe-1),(n-kPriebe));

%PerroneFit
theta=deg2rad(theta);
initialValuesPerrone=[scaling;prefSf;sigmaSf;prefTf;sigmaTf;theta];
[perroneParameters,perroneResiduals,perroneJacobian]=nlinfit(conditionMatrix,responseMatrix,@perronefit,initialValuesPerrone,statset('MaxIter',1000,'Display','iter'));
perroneParametersCI=nlparci(perroneParameters,perroneResiduals,perroneJacobian);
[perronePred,deltaPerroneCI] = nlpredci(@perronefit,conditionMatrix,perroneParameters,perroneResiduals,perroneJacobian);
nlintool(conditionMatrix,responseMatrix,@perronefit,initialValuesPerrone);
perronePredBounds=[responseMatrix perronePred deltaPerroneCI];
kPerrone=size(perroneParameters,1);
ssePerrone=sum(perroneResiduals.^2);
rsquaredPerrone=1-(ssePerrone/sst);
adj_rsquaredPerrone=1-((ssePerrone/(n-kPerrone-1))/(sst/(n-1)));
perroneAIC=n*log(ssePerrone/n)+2*kPerrone;
perroneAICc=perroneAIC+((2*(kPerrone+2)*(kPerrone+3))/(n-kPerrone-3));
Fperrone=(rsquaredPerrone/(kPerrone-1))/((1-rsquaredPerrone)/(n-kPerrone));
p_valuePerrone=1-fcdf(Fperrone,(kPerrone-1),(n-kPerrone));

%PerroneWIMFit
initialValuesPerroneWIM=[sigma;phase;k;A1;A2;A3;xc1;xs1;xc2;xs2;g;S;v;prefSf;alpha;delta];
[perroneWIMParameters,perroneWIMResiduals,perroneWIMJacobian]=nlinfit(conditionMatrix,normResponseMatrix,@perroneWIMfit,initialValuesPerroneWIM,statset('MaxIter',1000,'Display','iter'));
perroneWIMParametersCI=nlparci(perroneWIMParameters,perroneWIMResiduals,perroneWIMJacobian);
[perroneWIMPred,deltaPerroneWIMCI] = nlpredci(@perroneWIMfit,conditionMatrix,perroneWIMParameters,perroneWIMResiduals,perroneWIMJacobian);
nlintool(conditionMatrix,normResponseMatrix,@perroneWIMfit,initialValuesPerroneWIM);
perroneWIMPredBounds=[normResponseMatrix perroneWIMPred deltaPerroneWIMCI];
ssePerroneWIM=sum(perroneWIMResiduals.^2);
kPerroneWIM=size(perroneWIMParameters,1);
rsquaredPerroneWIM=1-(ssePerroneWIM/sstNorm);
adj_rsquaredPerroneWIM=1-((ssePerroneWIM/(n-kPerroneWIM-1))/(sstNorm/(n-1)));
perroneWIMAIC=n*log(ssePerroneWIM/n)+2*kPerroneWIM;
perroneWIMAICc=perroneWIMAIC+((2*(kPerroneWIM+2)*(kPerroneWIM+3))/(n-kPerroneWIM-3));
FperroneWIM=(rsquaredPerroneWIM/(kPerroneWIM-1))/((1-rsquaredPerroneWIM)/(n-kPerroneWIM));
p_valuePerroneWIM=1-fcdf(FperroneWIM,(kPerroneWIM-1),(n-kPerroneWIM));