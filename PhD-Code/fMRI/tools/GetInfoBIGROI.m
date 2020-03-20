function [Clusters_per_BIGROI, BIGROI_labels, BIGROI_of_Cluster] = GetInfoBIGROI()

Clusters_per_ROI = [...
    17    16    17    20     4     4    24    25     4     5     5     6  ...
    12    10     8     8     4     6    10    11     1     1    14    10  ...
     3     4     4     3     9     8     7     6     9    11     2     1  ...
     4     4     4     5     1     1    11     9     7     7    10    11  ...
     6     7    16    10     4     4    11    12    19    19    10    11  ...
    12     6     6     9     5     8    17    16     6     4     4     4  ...
     5     5     1     1     5     5     1     1    11    15     6     6  ...
    24    22     3     5    16    17];

nROI = numel( Clusters_per_ROI );

First_fcl_of_ROI = 1 + [0 cumsum(Clusters_per_ROI(1:nROI-1))];

Last_fcl_of_ROI  = cumsum( Clusters_per_ROI( 1 : nROI ) );

nCluster = sum( Clusters_per_ROI );

ROI_of_Cluster = zeros( 1,nCluster );

for iCluster = 1 : nCluster
    
    iROI = find( First_fcl_of_ROI <= iCluster & iCluster <= Last_fcl_of_ROI );
    
    ROI_of_Cluster(iCluster) = iROI;
end

Clusters_per_BIGROI = [...
   (17  + 16)  (17 +  20  +  4  +  4)  (24 +  25 +   4  +  5)   (5  +  6+ ...
    12  + 10  +  8  +  8)   (4  +  6)  (10 +  11)   (1  +  1)  (14  + 10+ ...
     3  +  4)   (4  +  3)   (9  +  8)   (7 +   6  +  9  + 11 +   2  +  1) ...
    (4  +  4  +  4  +  5)   (1  +  1)  (11 +   9)   (7  +  7)  (10  + 11) ...
    (6  +  7  + 16  + 10  +  4  +  4)  (11 +  12)  (19  + 19)  (10  + 11+ ...
    12  +  6)   (6  +  9)   (5  +  8)  (17 +  16)   (6  +  4)   (4  +  4) ...
    (5  +  5)   (1  +  1)   (5  +  5)   (1 +   1)  (11  + 15 +   6  +  6) ...
   (24  + 22  +  3  +  5)  (16  + 17)];
    
BIGROI_labels = {'Precentral', 'Frontal-Sup', 'Frontal-Mid', 'Frontal-Inf', 'Rolandic-Oper', 'Supp-Motor-Area', 'Olfactory', 'Frontal-Medial', ...
               'Rectus', 'Insula', 'Cingulum', 'Hippocampus', 'Amygdala', 'Calcarine', 'Cuneus', 'Lingual', 'Occipital', 'Fusiform', ...
               'Postcentral', 'Parietal', 'SupraMarginal', 'Angular', 'Precuneus', 'Paracentral-Lobule', 'Caudate', 'Putamen', 'Pallidum', ...
               'Thalamus', 'Heschl', 'Temporal-Sup', 'Temporal-Mid', 'Temporal-Inf'};
               

nBIGROI = numel( Clusters_per_BIGROI );
           
First_fcl_of_BIGROI = 1 + [0 cumsum(Clusters_per_BIGROI(1:nBIGROI-1))];

Last_fcl_of_BIGROI  = cumsum( Clusters_per_BIGROI( 1 : nBIGROI ) );

BIGROI_of_Cluster = zeros( 1,nCluster );

for iCluster = 1 : nCluster
    
    iBIGROI = find( First_fcl_of_BIGROI <= iCluster & iCluster <= Last_fcl_of_BIGROI );
    
    BIGROI_of_Cluster(iCluster) = iBIGROI;
end          


return;