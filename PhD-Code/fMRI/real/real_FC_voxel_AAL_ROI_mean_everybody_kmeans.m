function real_FC_voxel_AAL_ROI_mean_everybody_kmeans

%doWithOnlyOneAmountOfClusters;

doWithSeveralAmountOfClusters;

end

function doWithOnlyOneAmountOfClusters
        
%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_ROI = 1:90;
nROI = length(idx_ROI);

for iROI=1:nROI
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'.mat'));
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'.mat'));
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'.mat'));
    
    disp('Cluster Track');
    [FC_Track_KMeans.IdxClusters, FC_Track_KMeans.Tidx, FC_Track_KMeans.Ncluster] = ClusterWithKmeans_all( all_Track );
    
    disp('Cluster Passive');
    [FC_Passive_KMeans.IdxClusters, FC_Passive_KMeans.Tidx, FC_Passive_KMeans.Ncluster] = ClusterWithKmeans_all( all_Passive );
    
    disp('Cluster RestingState');
    [FC_RestingState_KMeans.IdxClusters, FC_RestingState_KMeans.Tidx, FC_RestingState_KMeans.Ncluster] = ClusterWithKmeans_all( all_RestingState );
                
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-','KMeans','.mat'),'FC_Track_KMeans');
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-','KMeans','.mat'),'FC_Passive_KMeans');
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-','KMeans','.mat'),'FC_RestingState_KMeans');
    
end

end

function doWithSeveralAmountOfClusters
        
%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_ROI = 1:90;
nROI = length(idx_ROI);

nClusters = 2:2:100;

for iROI=1:nROI
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'.mat'));
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'.mat'));
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'.mat'));
    
    for isizeCluster=1:length(nClusters)
        
        sizeCluster = nClusters(isizeCluster);
        
        disp(strcat('sizeCluster:',int2str(sizeCluster)));
        
        disp('Cluster Track');
        FC_Track_KMeans.cluster(isizeCluster).sizeCluster = sizeCluster;
        [FC_Track_KMeans.cluster(isizeCluster).IdxClusters, FC_Track_KMeans.cluster(isizeCluster).Tidx, FC_Track_KMeans.cluster(isizeCluster).Ncluster] = ClusterWithKmeans_all_size( all_Track, sizeCluster );

        disp('Cluster Passive');
        FC_Passive_KMeans.cluster(isizeCluster).sizeCluster = sizeCluster;
        [FC_Passive_KMeans.cluster(isizeCluster).IdxClusters, FC_Passive_KMeans.cluster(isizeCluster).Tidx, FC_Passive_KMeans.cluster(isizeCluster).Ncluster] = ClusterWithKmeans_all_size( all_Passive, sizeCluster );

        disp('Cluster RestingState');
        FC_RestingState_KMeans.cluster(isizeCluster).sizeCluster = sizeCluster;
        [FC_RestingState_KMeans.cluster(isizeCluster).IdxClusters, FC_RestingState_KMeans.cluster(isizeCluster).Tidx, FC_RestingState_KMeans.cluster(isizeCluster).Ncluster] = ClusterWithKmeans_all_size( all_RestingState, sizeCluster );
                
    end
    
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-','KMeans-size','.mat'),'FC_Track_KMeans');
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-','KMeans-size','.mat'),'FC_Passive_KMeans');
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-','KMeans-size','.mat'),'FC_RestingState_KMeans');
    
end

end


