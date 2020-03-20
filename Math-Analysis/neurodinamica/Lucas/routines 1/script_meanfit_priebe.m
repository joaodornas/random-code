% load ('population_data.mat','cells');
% ncells=size(cell_list,1);
% rsquaredVector=zeros(ncells,1);
% fitpVector=zeros(ncells,1);

for i=8:10%:ncells
    unit=cell_list(i);
    filename=strcat('2Dfits_',unit,'.mat');
    filename=filename{1};
    load(filename,'initialValuesPriebe','conditionMatrix');
    ind=strmatch(unit,cells);
    vector=strcat('responseVector',int2str(ind));
    load('partial_corr_bydata.mat',vector);
    response=eval(vector);
    nconditions=size(response,1);
    firstSF=conditionMatrix(1,1);
    firstTF=conditionMatrix(1,2);
    if nconditions==35
        nSF=5; nTF=7;
    elseif nconditions==36
        nSF=6; nTF=6;
    end
    conditions=zeros(nconditions,2);
    for sf=log2(firstSF):log2(firstSF)+nSF-1
        for tf=log2(firstTF):log2(firstTF)+nTF-1
            conditions((abs(log2(firstSF))+sf)*nTF+(abs(log2(firstTF))+tf+1),:)=[2^sf 2^tf];
        end
    end
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
    rsquaredVector(i)=rsquaredPriebe;
    fitpVector(i)=p_valuePriebe;
end

R2mean=nanmean(rsquaredVector);
R2sem=nanstd(rsquaredVector)./ncells;

  