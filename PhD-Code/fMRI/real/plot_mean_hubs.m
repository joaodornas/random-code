
clear all

load('FC-Voxels-AAL-ROI-corr-KMeans-MeaningFul-FC-Att-Stim-Only.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Functional-Parcels.mat');

FC_strength = zeros(1,length(MeaningFul_Attention)-1);
GC_FROM_strength = zeros(1,length(MeaningFul_Attention)-1);
GC_TO_strength = zeros(1,length(MeaningFul_Attention)-1); 

for iMean=2:length(MeaningFul_Attention)
    
    if ~isempty(length(MeaningFul_Attention{iMean,7}))
        FC_strength(iMean-1) = length(MeaningFul_Attention{iMean,7}); % FC - direct
    else
        FC_strength(iMean-1) = 0;
    end
    
    if ~isempty(length(MeaningFul_Attention{iMean,8}))
        FC_strength(iMean-1) = FC_strength(iMean-1) + length(MeaningFul_Attention{iMean,8}); % FC - indirect
    else
        FC_strength(iMean-1) = FC_strength(iMean-1) + 0;
    end

    if ~isempty(length(MeaningFul_Attention{iMean,11}))
        GC_FROM_strength(iMean-1) = length(MeaningFul_Attention{iMean,11}); % GC - direct - FROM
    else
        GC_FROM_strength(iMean-1) = 0;
    end
    
    if ~isempty(length(MeaningFul_Attention{iMean,12}))
        GC_FROM_strength(iMean-1) = GC_FROM_strength(iMean-1) + length(MeaningFul_Attention{iMean,12}); % GC - indirect - FROM
    else
        GC_FROM_strength(iMean-1) = GC_FROM_strength(iMean-1) + 0;
    end
    
    if ~isempty(length(MeaningFul_Attention{iMean,15}))
        GC_TO_strength(iMean-1) = length(MeaningFul_Attention{iMean,15}); % GC - direct - TO
    else
        GC_TO_strength(iMean-1) = 0;
    end
    
    if ~isempty(length(MeaningFul_Attention{iMean,16}))
        GC_TO_strength(iMean-1) = GC_TO_strength(iMean-1) + length(MeaningFul_Attention{iMean,16}); % GC - indirect - TO
    else
        GC_TO_strength(iMean-1) = GC_TO_strength(iMean-1) + 0;
    end
    
end

[s i] = sort(FC_strength,'descend');
idx_cluster(1) = i(1)+1;

[s i] = sort(GC_FROM_strength,'descend');
idx_cluster(2) = i(1)+1;

[s i] = sort(GC_TO_strength,'descend');
idx_cluster(3) = i(1)+1;

all_clusters = 0;
for iCluster=1:length(idx_cluster)
    
    clusters.FC_direct = MeaningFul_Attention{idx_cluster(iCluster),7};
    clusters.FC_indirect = MeaningFul_Attention{idx_cluster(iCluster),8};
    
    clusters.GC_FROM_direct = MeaningFul_Attention{idx_cluster(iCluster),11};
    clusters.GC_FROM_indirect = MeaningFul_Attention{idx_cluster(iCluster),12};
    
    clusters.GC_TO_direct = MeaningFul_Attention{idx_cluster(iCluster),15};
    clusters.GC_TO_indirect = MeaningFul_Attention{idx_cluster(iCluster),16};
    
    all_clusters = [clusters.FC_direct, clusters.FC_indirect, clusters.GC_FROM_direct, clusters.GC_FROM_indirect, clusters.GC_TO_direct, clusters.GC_TO_indirect];
    
end