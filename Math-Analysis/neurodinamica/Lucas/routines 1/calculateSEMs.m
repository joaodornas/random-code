for i=1:6
    SEM=zeros(1,6);
    sf=2^(i-3);
    for j=1:6
        tf=2^(j-3);
        a=find(speed_conditionVector(:,1)==sf & speed_conditionVector(:,2)==tf);
        vector=speed_responseVector(a);
        SEM(j)=nanstd(vector)./sqrt(size(a,1));
    end
    assignin('base',strcat('speed_Sf',int2str(i),'_SEM'),SEM);
end