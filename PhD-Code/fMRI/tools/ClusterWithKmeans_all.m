function [IdxClusters, Tidx, Ncluster] = ClusterWithKmeans_all( rho_input, nClusters )

tic

nVoxels = size(rho_input,1);

% %% keep only significant correlations
%     
%     pcriterion = 0.01;
% 
%     rho_input( pval_input > pcriterion ) = 0;
    
%% find NaNs

    kk = find(~any(rho_input,2));
    
    rho_input(kk,:) = [];
    rho_input(:,kk) = [];

%% compute the clusters

    num_Clust = nClusters;

    %[Tidx,C,sumd,D] = kmeans(rho_input, num_Clust , 'distance', 'correlation', 'display', 'final','replicate',20);
    [Tidx,C,sumd,D] = kmeans(rho_input, num_Clust , 'distance', 'correlation', 'display', 'off','replicate',20);

    Ncluster = max(Tidx);

%% repositions voxels clusters with NaNs

IdxClusters = zeros(nVoxels,1);
ik = 0;

for iVoxel=1:nVoxels
    
    if find(kk==iVoxel)
        
        IdxClusters(iVoxel) = NaN;
        
    else
        
        ik = ik + 1;
        
        IdxClusters(iVoxel) = Tidx(ik);

    end
    
end

toc

end