function [IdxClusters, Ncluster, distance] = ClusterWithKmeans_no_symmetrical( rho_input, pval_input, nClusters )

tic

%% keep only significant correlations
    
    pcriterion = 0.01;

    rho_input( pval_input > pcriterion ) = 0;
    
    rho_input( isnan(rho_input) ) = 0;
    
    distance = 'sqEuclidean';
    %distance = 'correlation';

%% compute the clusters

    %[Tidx,C,sumd,D] = kmeans(rho_input, num_Clust , 'distance', 'correlation', 'display', 'final','replicate',20);
    [Tidx,C,sumd,D] = kmeans(rho_input, nClusters , 'distance', distance, 'display', 'off','replicate',20);

    Ncluster = max(Tidx);
    
    IdxClusters = Tidx;

toc

end