function [p_value]=fitp(goodness1,output1)

%calculates fit p value using the F test. Inputs are goodness1 and output1,
%which are saved to workspace from cftool.
%[p_value]=fitp(goodness1,output1)

F=(goodness1.rsquare/(output1.numparam-1))/((1-goodness1.rsquare)/(output1.numobs-output1.numparam));
p_value=1-fcdf(F,(output1.numparam-1),(output1.numobs-output1.numparam));