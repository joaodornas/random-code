for i=1:6
    sf=2^(i-3);
    a=find(neg_conditionVector(:,1)==sf);
    trials=neg_responseVector(a);
    conditions=neg_conditionVector(a,2)./sf;
    assignin('base',strcat('neg_SP',int2str(i),'_trials'),trials);
    assignin('base',strcat('neg_SP',int2str(i),'_conditions'),conditions);
end