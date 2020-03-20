function g = perrone(beta,x)

%the model is a two-dimensional gaussian of the form:
%)g(sf,tf)=a*(exp(-(power(((sf-sf0)cos(theta)+(tf-tf0)sin(theta)),2))./(power(sigmaSf
%,2)).*(exp(-(power(-(sf-sf0)sin(theta)+(tf-tf0)cos(theta)),2)./(power(sigm
%aTf,2)))
%independent variables: sf,tf
%free parameters:c (spont. act., constant); sf0 (preferred spatial frequency);
%sigmaSf; tf0 (preferred temporal frequency); sigmaTf; theta (orientation
%angle)
%Perrone & Thiele (2001) Nat Neurosci

scaling=beta(1)
prefSf=beta(2)
sigmaSf=beta(3)
prefTf=beta(4)
sigmaTf=beta(5)
theta=beta(6)
sf=x(:,1)
tf=x(:,2)

g=scaling.*(exp(-(power(((sf-prefSf).*cos(theta)+(tf-prefTf).*sin(theta)),2))./(power(sigmaSf,2)))).*(exp(-(power((-(sf-prefSf).*sin(theta)+(tf-prefTf).*cos(theta)),2))./(power(sigmaTf,2))))