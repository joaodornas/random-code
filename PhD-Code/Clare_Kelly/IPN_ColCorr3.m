%% Compute Pair-correlation between two sample matrices.
function R = IPN_ColCorr3(X, Y)
%
% XiNian.Zuo@nyumc.org

[numSamp1 numVar1] = size(X); % n*p1
[numSamp2 numVar2] = size(Y); % n*p2
numVar = min(numVar1, numVar2);

if (numSamp1 ~= numSamp2)
    disp('The two matices must have the same size of rows!')
else
    X = (X - repmat(mean(X), numSamp1, 1))./repmat(std(X, 0, 1), numSamp1, 1);
    Y = (Y - repmat(mean(Y), numSamp1, 1))./repmat(std(Y, 0, 1), numSamp1, 1);
    R = X' * Y / (numSamp1 - 1);
end
%if numVar1*numVar2 ~= 1
%    R0 = R(1:numVar, 1:numVar);
%    R0 = R0 - diag(diag(R0));
%    R(1:numVar, 1:numVar) = R0;
%    clear R0
%end