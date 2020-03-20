function [pvalue, pA_significance, pB_significance] = Binostats( NA, OnA, NB, OnB )

%% given are two event counts, NA and NB, and two numbers of On events, 0 <= OnA <= NA and 0 <= OnB <= NB

%% if On events are distributed identically in both sets, then the number of On events will be binomially distributed

Ntotal = NA + NB;
Ontotal = OnA + OnB;

p = Ontotal / Ntotal;

%% check binomial probabilities and cumulative binomial probabilties

POn = binopdf( 0:Ntotal, Ntotal, p );

cumulative_Pon = betainc( 1-p, Ntotal:-1:0, 1:Ntotal+1);

% figure;
% 
% hold on;
% stairs( 0:Ntotal, POn, 'r' );
% stairs( 0:Ntotal, cumulative_Pon, 'b' );
% hold off;

%% compute probability of deviation from null hypothesis

pA = OnA / NA;   pB = OnB / NB;

if pA < p    
    pA_significance = betainc( 1-p, Ntotal*(1-pA), Ntotal*pA + 1, 'lower' );   
else
    pA_significance = betainc( 1-p, Ntotal*(1-pA), Ntotal*pA + 1, 'upper' );     
end

if pB < p    
    pB_significance = betainc( 1-p, Ntotal*(1-pB), Ntotal*pB + 1, 'lower' );   
else
    pB_significance = betainc( 1-p, Ntotal*(1-pB), Ntotal*pB + 1, 'upper' );     
end
    
pvalue = max( pA_significance, pB_significance );

return;