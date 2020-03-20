function [AincreaseFraction, RincreaseFraction, AdecreaseFraction, RdecreaseFraction, Achanges, Rchanges, Aincrease, Rincrease, Adecrease, Rdecrease, Dfisher, Afisher, Rfisher, Asorted, Rsorted, target_order, A_source_order, R_source_order, source_order, Ncomponent, cluster_assignment, Ncluster] = getFCDifferencePositiveNegativeFraction(rho_att,pval_att,rho_rest,pval_rest,shouldIDoCluster)

    pcriterion = 0.01;
    threshold = 0.35;

    rho_att( pval_att > pcriterion ) = 0;
    rho_rest( pval_rest > pcriterion ) = 0;
    
    rho_att( rho_att == 1 ) = 0.99;
    
    rho_rest( rho_rest == 1 ) = 0.99;
    
    [Ncomponent, ~] = size( rho_rest );
    
    if shouldIDoCluster
        
        [cluster_assignment, Tidx, Ncluster] = ClusterWithKmeans( rho_rest, pval_rest );  % kmeans cluster original correlation matrix for resting)
 
        source_order = 1:Ncomponent;
        [target_order, cluster_idx_sorted] = TargetOrder( cluster_assignment, Ncluster );
    
        [Rsorted, R_source_order] = SortMatrix( rho_rest, source_order, target_order, Ncomponent );
    
        [Asorted, A_source_order] = SortMatrix( rho_att, source_order, target_order, Ncomponent );
    
    else
        
        cluster_assignment = [];
        Ncluster = [];
        cluster_idx_sorted = [];
        source_order = 1:Ncomponent;
        target_order = source_order;
        Rsorted = rho_rest;
        Asorted = rho_att;
        R_source_order = source_order;
        A_source_order = source_order;
        
    end
        
    Afisher = 0.5 * log( ( 1 + Asorted ) ./ ( 1 - Asorted ) );
    
    Rfisher = 0.5 * log( ( 1 + Rsorted ) ./ ( 1 - Rsorted ) );
    
    Dfisher = Afisher - Rfisher;
    
    Adecrease = Asorted;
    Adecrease(  Dfisher > -threshold ) = 0;
    Adecrease(  Adecrease > 0 ) = 0;         % limit to negative
    
    Aincrease = Asorted;
    Aincrease(  Dfisher <  threshold ) = 0;
    Aincrease(  Aincrease < 0 ) = 0;         % limit to positive
    
    Rdecrease = Rsorted;
    Rdecrease(  Dfisher > -threshold ) = 0;
    Rdecrease(  Rdecrease < 0 ) = 0;         % limit to positive
       
    Rincrease = Rsorted;
    Rincrease(  Dfisher <  threshold ) = 0;
    Rincrease(  Rincrease > 0 ) = 0;         % limit to negative
    
    Achanges = Asorted;
    Achanges( (Dfisher < threshold) & (Dfisher > -threshold)) = 0;
    
    Rchanges = Rsorted;
    Rchanges( (Dfisher < threshold) & (Dfisher > -threshold)) = 0;
    
    significant = Adecrease ~= 0 & ~isnan(Adecrease);
    AdecreaseFraction = sum(significant,1)./Ncomponent;

    significant = Aincrease ~= 0  & ~isnan(Aincrease);
    AincreaseFraction = sum(significant,1)./Ncomponent;
    
    significant = Rdecrease ~= 0  & ~isnan(Rdecrease);
    RdecreaseFraction = sum(significant,1)./Ncomponent;

    significant = Rincrease ~= 0  & ~isnan(Rincrease);
    RincreaseFraction = sum(significant,1)./Ncomponent;
    
end
