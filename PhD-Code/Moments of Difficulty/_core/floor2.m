function z = floor2(x,y)
%floor2 rounds number to nearest multiple of arbitrary precision.
%   Z = floor2(X,Y) rounds X to nearest multiple of Y.
%
% See also ROUND.

%% defensive programming
error(nargchk(2,2,nargin))
error(nargoutchk(0,1,nargout))
if numel(y)>1
  error('Y must be scalar')
end

%%
z = floor(x/y)*y;
