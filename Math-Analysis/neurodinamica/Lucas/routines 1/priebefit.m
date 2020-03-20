function [response] = priebefit(beta,x)

%the model is a two-dimensional gaussian of the form:
%a*exp(-((log2(sf)-log2(sf0))^2)/(2*(sigmaSf^2)))*(exp(-((log2(tf)-(epsilon
%*(log2(sf)-log2(sf0))+log2(tf0)))^2)/(2*((sigmaTf+skew*(log2(tf)-(epsilon*
%(log2(sf)-log2(sf0))+log2(tf0))))^2)))-exp(-1/skew^2))+c
%independent variables: sf,tf
%free parameters:a (peak response); sf0 (preferred spatial frequency); sigmaSf; epsilon (speed
%tuning index); tf0 (preferred temporal frequency); sigmaTf; skew
%Priebe et al (2006) J Neuroscience
%
%arguments: [response] = priebe(beta,x)
%where inputs are vectors beta (function parameters) and x (independent variables
%sf and tf, respectively 1st and 2nd columns) and output is unit response

scaling=beta(1);
prefSf=beta(2);
sigmaSf=beta(3);
prefTf=beta(4);
sigmaTf=beta(5);
skew=beta(6);
dependence=beta(7);
sf=x(:,1);
tf=x(:,2);

response=scaling.*exp(-(power((log2(sf)-log2(prefSf)),2))/(2*sigmaSf^2)).*(exp(-(power((log2(tf)-(dependence.*(log2(sf)-log2(prefSf))+log2(prefTf))),2))./(2.*(power((sigmaTf+skew.*(log2(tf)-(dependence.*(log2(sf)-log2(prefSf))+log2(prefTf)))),2))))-exp(-1/skew^2));