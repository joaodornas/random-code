function [nROI, ROI_labels, ROI_labels_inverted, nCluster, Clusters_per_ROI, ROI_of_Cluster, First_fcl_of_ROI, Last_fcl_of_ROI] = LoadClusterInfo()

load 'Functional-Clusters-All-Info.mat';

nROI = numel( ROI );

Clusters_per_ROI = nan( 1, nROI );

for iROI = 1 : nROI
    
    ROI_labels{ iROI } = [num2str( iROI ) '   ' ROI( iROI ).label];
    ROI_labels_inverted{ nROI - iROI + 1 } = [num2str( iROI ) '   ' ROI( iROI ).label];
    Clusters_per_ROI( iROI ) = ROI( iROI ).nClusters;
    
end

char( ROI_labels )

First_fcl_of_ROI = 1 + [0 cumsum(Clusters_per_ROI(1:nROI-1))];

Last_fcl_of_ROI  = cumsum( Clusters_per_ROI( 1 : nROI ) );

% store ROI number for each cluster

nCluster = sum( Clusters_per_ROI );

ROI_of_Cluster = zeros( 1,nCluster );

for iCluster = 1 : nCluster
    
    iROI = find( First_fcl_of_ROI <= iCluster & iCluster <= Last_fcl_of_ROI );
    
    ROI_of_Cluster(iCluster) = iROI;
end


