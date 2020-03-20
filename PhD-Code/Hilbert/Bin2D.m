function [Count, BinMid]= Bin2D(DataSource, BinsN, Limits)
 
%% preparing limits, if not given
if (~exist('Limits', 'var') || isempty(Limits))
  for iDim= 1:2,
    Limits(iDim, 1:2)= [min(DataSource(iDim, :)) max(DataSource(iDim, :))];
  end;
end;


%% preparing bins
BinLo= cell(1, 2);
BinHi= cell(1, 2);
BinMid= cell(1, 2);
for iDim= 1:2,
  BinStep= (Limits(iDim, 2)-Limits(iDim, 1))/(BinsN(iDim)-1);
  BinLo{iDim}= Limits(iDim, 1):BinStep:Limits(iDim, 2);
  BinHi{iDim}= BinLo{iDim}+BinStep;
  BinHi{iDim}(end)= BinHi{iDim}(end)+0.001; %% so we can do >=Lo & <Hi
  BinMid{iDim}= BinLo{iDim}+BinStep/2;
end;
  
%% binning
Count= nan(numel(BinLo{1}), numel(BinLo{2}));
for iBin1= 1:numel(BinLo{1}),
  for iBin2= 1:numel(BinLo{2}),
    iCurrent= find(DataSource(1, :)>=BinLo{1}(iBin1) & DataSource(1, :)<BinHi{1}(iBin1) & DataSource(2, :)>=BinLo{2}(iBin2) & DataSource(2, :)<BinHi{2}(iBin2));
    Count(iBin1, iBin2)= numel(iCurrent);
  end;
end;

end