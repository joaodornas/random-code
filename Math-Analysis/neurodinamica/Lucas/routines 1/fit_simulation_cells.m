epsilonVector=zeros(1440,1);
epsilonVector_CI=zeros(1440,2);
fitpVector=zeros(1440,1);
rsquareVector=zeros(1440,1);
scaling=1;
sigmaSf=2;
sigmaTf=2.5;
skew=-0.2;
dependence=1;

for i=1:1440
    load('simulation_results',strcat('WIMvector_',int2str(i)),strcat('conditionVector_',int2str(i)));
    response=eval(strcat('WIMvector_',int2str(i)));
    conditions=eval(strcat('conditionVector_',int2str(i)));
    clear(strcat('WIMvector_',int2str(i)));
    clear(strcat('conditionVector_',int2str(i)));
    conditions=2.^conditions;
    t=isreal(response);
    if t==0
        response=abs(response);
    end
    maxresp=max(response);
    a=find(response==maxresp,1,'first');
    prefSf=conditions(a,1);
    prefTf=conditions(a,2);
    initialValuesPriebe=[scaling;prefSf;sigmaSf;prefTf;sigmaTf;skew;dependence];
    [priebeParameters,priebeResiduals,priebeJacobian]=nlinfit(conditions,response,@priebefit,initialValuesPriebe,statset('MaxIter',1000));
    t1=isreal(priebeParameters);
    if t1==0
        priebeParameters=abs(priebeParameters);
        priebeResiduals=abs(priebeResiduals);
        priebeJacobian=abs(priebeJacobian);
    end
    priebeParametersCI=nlparci(priebeParameters,priebeResiduals,priebeJacobian);
    meanResponse=mean(response);
    sst=sum((response-meanResponse).^2);
    kPriebe=size(priebeParameters,1);
    n=size(response,1);
    ssePriebe=sum(priebeResiduals.^2);
    rsquaredPriebe=1-(ssePriebe/sst);
    Fpriebe=(rsquaredPriebe/(kPriebe-1))/((1-rsquaredPriebe)/(n-kPriebe));
    p_valuePriebe=1-fcdf(Fpriebe,(kPriebe-1),(n-kPriebe));
    assignin('base',strcat('fitParameters',int2str(i)),priebeParameters);
    rsquareVector(i)=rsquaredPriebe;
    fitpVector(i)=p_valuePriebe;
    if p_valuePriebe<=0.05
        epsilonVector(i)=priebeParameters(7);
        epsilonVector_CI(i,:)=priebeParametersCI(7,:);
        if priebeParametersCI(7,1)<1 && priebeParametersCI(7,2)<0
            if i==1
                classifVector=cellstr('negative');
            else classifVector=[classifVector;'negative'];
            end
        elseif priebeParametersCI(7,1)<0 && (priebeParametersCI(7,2)>=0 && priebeParametersCI(7,2)<1)
            if i==1
                classifVector=cellstr('independent');
            else classifVector=[classifVector;'independent'];
            end
        elseif (priebeParametersCI(7,1)>0 && priebeParametersCI(7,1)<1) && (priebeParametersCI(7,2)>0 && priebeParametersCI(7,2)<1)
            if i==1
                classifVector=cellstr('intermediate');
            else classifVector=[classifVector;'intermediate'];
            end
        elseif priebeParametersCI(7,1)>0 && priebeParametersCI(7,1)<1 && priebeParametersCI(7,2)>=1
            if i==1
                classifVector=cellstr('speed');
            else classifVector=[classifVector;'speed'];
            end
        else if i==1
                classifVector=cellstr('undetermined');
            else classifVector=[classifVector;'undetermined'];
            end
        end
    else epsilonVector(i)=NaN;
        epsilonVector_CI(i,:)=NaN;
        if i==1
           classifVector=cellstr('NaN');
        else classifVector=[classifVector;'NaN'];
        end
    end
end

speedi=strmatch('speed',classifVector);
indepi=strmatch('independent',classifVector);
intermediatei=strmatch('intermediate',classifVector);
negativei=strmatch('negative',classifVector);
selectedi=[speedi;indepi;intermediatei;negativei];
selectedi=sort(selectedi,'ascend');
epsilonVector_selected=epsilonVector(selectedi);
speed_proportion=size(speedi,1)/size(selectedi,1);
negative_proportion=size(negativei,1)/size(selectedi,1);
indep_proportion=size(indepi,1)/size(selectedi,1);
intermediate_proportion=size(intermediatei,1)/size(selectedi,1);

save('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\PerroneWIM simulation\simulation_fits')      