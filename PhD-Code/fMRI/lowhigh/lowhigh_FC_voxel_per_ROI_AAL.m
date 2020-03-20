function lowhigh_FC_voxel_per_ROI_AAL


settings_jan_0805;
%setting_elena_2905;

%doTheMath(settings);

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26 27 28];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [17 18 19 20 57 58 59 60 61 62 63 64 65 66 67 68 69 70];

% doRegionByRegion(settings,idx_frontal);
% doRegionByRegion(settings,idx_occipital);
% doRegionByRegion(settings,idx_parietal);

% doRegionByRegion_half(settings,idx_frontal);
% doRegionByRegion_half(settings,idx_occipital);
% doRegionByRegion_half(settings,idx_parietal);

% computeClustersRegionByRegion(settings,idx_frontal);
% computeClustersRegionByRegion(settings,idx_occipital);
% computeClustersRegionByRegion(settings,idx_parietal);

% plotFCClusteredRegionByRegion(settings,idx_frontal);
% plotFCClusteredRegionByRegion(settings,idx_occipital);
% plotFCClusteredRegionByRegion(settings,idx_parietal);

idx_ROI = 66; %% ANGULAR R
%computeClustersRegionByRegion_half(settings,idx_ROI);
%plotFCClusteredRegionByRegion_half(settings,idx_ROI);
plotFCClusteredRegionByRegion_half_SameOrder(settings,idx_ROI)

%idx_ROI = [idx_frontal,idx_parietal,idx_occipital];
%computeFCMeanSQClustersAllROIs(settings,idx_ROI);
%plotFCMeanSQClustersAllROIs(settings,idx_ROI);

% computeClustersRegionByRegion_half(settings,idx_frontal);
% computeClustersRegionByRegion_half(settings,idx_occipital);
% computeClustersRegionByRegion_half(settings,idx_parietal);

%idx_ROI = [idx_frontal,idx_parietal,idx_occipital];
% computeFCMeanSQClustersAllROIs_half(settings,idx_ROI);
% plotFCMeanSQClustersAllROIs_half(settings,idx_ROI);

%computePearsonBtwROI(settings);

% area1_ID = 13; %% FRONTAL INF TRI L
% area2_ID = 66; %% ANGULAR R
% plotFCClusterTwoRegions_all_rhos(settings,area1_ID,area2_ID);

% area1_ID = 13; %% FRONTAL INF TRI L
% area2_ID = 66; %% ANGULAR R
% plotFCBtwConditionFractionDistance(settings,area1_ID,area2_ID);

% area1_ID = 13; %% FRONTAL INF TRI L
% area2_ID = 66; %% ANGULAR R
% plotFCDistributionAndDifference(settings,area1_ID,area2_ID);

% area_ID = 66; %% ANGULAR R
% plotFCDifferencePositiveNegativeFraction(settings,area_ID);
% area_ID = 13; %% FRONTAL INF TRI L
% plotFCDifferencePositiveNegativeFraction(settings,area_ID);

% idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26 27 28];
% idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
% idx_parietal = [17 18 19 20 57 58 59 60 61 62 63 64 65 66 67 68 69 70];
% idx_areas = [idx_frontal,idx_occipital,idx_parietal];
% computeAllFractionsFCDifferencePositiveNegative(settings,idx_areas);

%generate3DimgFractions(settings);



end


function doTheMath(settings)

subject = settings.folders.subject;

%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

lowhigh_load_all_data_FSL;

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

idx_ROI = [59 60 61 62 65 66];

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    disp(label_ROI{iROI});
    
    disp('High Attention');
    [FC_Voxels.ROI(iROI).run(1).rho_MOT4, FC_Voxels.ROI(iROI).run(1).pval_MOT4] = getVoxelCorrelations(MOT4Run1,AAL_img,AAL_ROI(idx_ROI(iROI)).ID);
    [FC_Voxels.ROI(iROI).run(2).rho_MOT4, FC_Voxels.ROI(iROI).run(2).pval_MOT4] = getVoxelCorrelations(MOT4Run2,AAL_img,AAL_ROI(idx_ROI(iROI)).ID);

%     disp('Low Attention');
%     [FC_Voxels.ROI(iROI).run(1).rho_MOT2, FC_Voxels.ROI(iROI).run(1).pval_MOT2] = getCorrelations(MOT2Run1,AAL_img,AAL_ROI(idx_ROI(iROI)).ID);
%     [FC_Voxels.ROI(iROI).run(2).rho_MOT2, FC_Voxels.ROI(iROI).run(2).pval_MOT2] = getCorrelations(MOT2Run2,AAL_img,AAL_ROI(idx_ROI(iROI)).ID);
    
    disp('Resting State');
    [FC_Voxels.ROI(iROI).run(1).rho_RestingState, FC_Voxels.ROI(iROI).run(1).pval_RestingState] = getVoxelCorrelations(RestingStateRun1,AAL_img,AAL_ROI(idx_ROI(iROI)).ID);
    [FC_Voxels.ROI(iROI).run(2).rho_RestingState, FC_Voxels.ROI(iROI).run(2).pval_RestingState] = getVoxelCorrelations(RestingStateRun2,AAL_img,AAL_ROI(idx_ROI(iROI)).ID);
    
end

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL','.mat'),'subject','file','idx_ROI','label_ROI','FC_Voxels');

end

function doRegionByRegion_half(settings,idx_ROI)

subject = settings.folders.subject;

%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

lowhigh_load_all_data_FSL;

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
%load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(label_ROI{iROI});
    
    disp('High Attention - FH');
    [FC_Voxels.run(1).rho_MOT4_FH, FC_Voxels.run(1).pval_MOT4_FH] = getVoxelCorrelations_half(MOT4Run1,AAL_img,AAL_ROI(idx_ROI(iROI)).ID,'FH');
    [FC_Voxels.run(2).rho_MOT4_FH, FC_Voxels.run(2).pval_MOT4_FH] = getVoxelCorrelations_half(MOT4Run2,AAL_img,AAL_ROI(idx_ROI(iROI)).ID,'FH');
 
    disp('Resting State - FH');
    [FC_Voxels.run(1).rho_RestingState_FH, FC_Voxels.run(1).pval_RestingState_FH] = getVoxelCorrelations_half(RestingStateRun1,AAL_img,AAL_ROI(idx_ROI(iROI)).ID,'FH');
    [FC_Voxels.run(2).rho_RestingState_FH, FC_Voxels.run(2).pval_RestingState_FH] = getVoxelCorrelations_half(RestingStateRun2,AAL_img,AAL_ROI(idx_ROI(iROI)).ID,'FH');
    
    disp('High Attention - SH');
    [FC_Voxels.run(1).rho_MOT4_SH, FC_Voxels.run(1).pval_MOT4_SH] = getVoxelCorrelations_half(MOT4Run1,AAL_img,AAL_ROI(idx_ROI(iROI)).ID,'SH');
    [FC_Voxels.run(2).rho_MOT4_SH, FC_Voxels.run(2).pval_MOT4_SH] = getVoxelCorrelations_half(MOT4Run2,AAL_img,AAL_ROI(idx_ROI(iROI)).ID,'SH');
 
    disp('Resting State - SH');
    [FC_Voxels.run(1).rho_RestingState_SH, FC_Voxels.run(1).pval_RestingState_SH] = getVoxelCorrelations_half(RestingStateRun1,AAL_img,AAL_ROI(idx_ROI(iROI)).ID,'SH');
    [FC_Voxels.run(2).rho_RestingState_SH, FC_Voxels.run(2).pval_RestingState_SH] = getVoxelCorrelations_half(RestingStateRun2,AAL_img,AAL_ROI(idx_ROI(iROI)).ID,'SH');
    
    save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Half-',area_label,'.mat'),'subject','file','area_label','FC_Voxels','-v7.3');

    clear FC_Voxels
    
end

end

function doRegionByRegion(settings,idx_ROI)

subject = settings.folders.subject;

%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

lowhigh_load_all_data_FSL;

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(label_ROI{iROI});
    
    disp('High Attention');
    [FC_Voxels.run(1).rho_MOT4, FC_Voxels.run(1).pval_MOT4] = getVoxelCorrelations(MOT4Run1,AAL_img,AAL_ROI(idx_ROI(iROI)).ID);
    [FC_Voxels.run(2).rho_MOT4, FC_Voxels.run(2).pval_MOT4] = getVoxelCorrelations(MOT4Run2,AAL_img,AAL_ROI(idx_ROI(iROI)).ID);

%     disp('Low Attention');
%     [FC_Voxels.run(1).rho_MOT2, FC_Voxels.run(1).pval_MOT2] = getCorrelations(MOT2Run1,AAL_img,AAL_ROI(idx_ROI(iROI)).ID);
%     [FC_Voxels.run(2).rho_MOT2, FC_Voxels.run(2).pval_MOT2] = getCorrelations(MOT2Run2,AAL_img,AAL_ROI(idx_ROI(iROI)).ID);
    
    disp('Resting State');
    [FC_Voxels.run(1).rho_RestingState, FC_Voxels.run(1).pval_RestingState] = getVoxelCorrelations(RestingStateRun1,AAL_img,AAL_ROI(idx_ROI(iROI)).ID);
    [FC_Voxels.run(2).rho_RestingState, FC_Voxels.run(2).pval_RestingState] = getVoxelCorrelations(RestingStateRun2,AAL_img,AAL_ROI(idx_ROI(iROI)).ID);
    
    save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area_label,'.mat'),'subject','file','area_label','FC_Voxels');

    clear FC_Voxels
    
end

end

function computeClustersRegionByRegion(settings,idx_ROI)

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(label_ROI{iROI});
    
    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area_label,'.mat'));

    [FC_Kmeans.run(1).high_IdxClusters, FC_Kmeans.run(1).high_Tidx, FC_Kmeans.run(1).high_Nclusters] = ClusterWithKmeans( FC_Voxels.run(1).rho_MOT4, FC_Voxels.run(1).pval_MOT4 );
    [FC_Kmeans.run(2).high_IdxClusters, FC_Kmeans.run(2).high_Tidx, FC_Kmeans.run(2).high_Nclusters] = ClusterWithKmeans( FC_Voxels.run(2).rho_MOT4, FC_Voxels.run(2).pval_MOT4 );

    [FC_Kmeans.run(1).rest_IdxClusters, FC_Kmeans.run(1).rest_Tidx, FC_Kmeans.run(1).rest_Nclusters] = ClusterWithKmeans( FC_Voxels.run(1).rho_RestingState, FC_Voxels.run(1).pval_RestingState );
    [FC_Kmeans.run(2).rest_IdxClusters, FC_Kmeans.run(1).rest_Tidx, FC_Kmeans.run(2).rest_Nclusters] = ClusterWithKmeans( FC_Voxels.run(2).rho_RestingState, FC_Voxels.run(2).pval_RestingState );
    
    save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area_label,'-Kmeans','.mat'),'area_label','FC_Kmeans');

    clear FC_Kmeans
    
end

end

function computeClustersRegionByRegion_half(settings,idx_ROI)

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
%load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(label_ROI{iROI});
    
    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Half-',area_label,'.mat'));

    [FC_Kmeans.run(1).high_IdxClusters_FH, FC_Kmeans.run(1).high_Tidx_FH, FC_Kmeans.run(1).high_Nclusters_FH] = ClusterWithKmeans( FC_Voxels.run(1).rho_MOT4_FH, FC_Voxels.run(1).pval_MOT4_FH );
    [FC_Kmeans.run(2).high_IdxClusters_FH, FC_Kmeans.run(2).high_Tidx_FH, FC_Kmeans.run(2).high_Nclusters_FH] = ClusterWithKmeans( FC_Voxels.run(2).rho_MOT4_FH, FC_Voxels.run(2).pval_MOT4_FH );

    [FC_Kmeans.run(1).rest_IdxClusters_FH, FC_Kmeans.run(1).rest_Tidx_FH, FC_Kmeans.run(1).rest_Nclusters_FH] = ClusterWithKmeans( FC_Voxels.run(1).rho_RestingState_FH, FC_Voxels.run(1).pval_RestingState_FH );
    [FC_Kmeans.run(2).rest_IdxClusters_FH, FC_Kmeans.run(1).rest_Tidx_FH, FC_Kmeans.run(2).rest_Nclusters_FH] = ClusterWithKmeans( FC_Voxels.run(2).rho_RestingState_FH, FC_Voxels.run(2).pval_RestingState_FH );
    
    [FC_Kmeans.run(1).high_IdxClusters_SH, FC_Kmeans.run(1).high_Tidx_SH, FC_Kmeans.run(1).high_Nclusters_SH] = ClusterWithKmeans( FC_Voxels.run(1).rho_MOT4_SH, FC_Voxels.run(1).pval_MOT4_SH );
    [FC_Kmeans.run(2).high_IdxClusters_SH, FC_Kmeans.run(2).high_Tidx_SH, FC_Kmeans.run(2).high_Nclusters_SH] = ClusterWithKmeans( FC_Voxels.run(2).rho_MOT4_SH, FC_Voxels.run(2).pval_MOT4_SH );

    [FC_Kmeans.run(1).rest_IdxClusters_SH, FC_Kmeans.run(1).rest_Tidx_SH, FC_Kmeans.run(1).rest_Nclusters_SH] = ClusterWithKmeans( FC_Voxels.run(1).rho_RestingState_SH, FC_Voxels.run(1).pval_RestingState_SH );
    [FC_Kmeans.run(2).rest_IdxClusters_SH, FC_Kmeans.run(1).rest_Tidx_SH, FC_Kmeans.run(2).rest_Nclusters_SH] = ClusterWithKmeans( FC_Voxels.run(2).rho_RestingState_SH, FC_Voxels.run(2).pval_RestingState_SH );
    
    save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Half-',area_label,'-Kmeans','.mat'),'area_label','FC_Kmeans');

    clear FC_Kmeans
    
end

end

function [IdxClusters, Tidx, Ncluster] = ClusterWithKmeans( rho_input, pval_input )

tic

nVoxels = size(rho_input,1);

%% keep only significant correlations
    
    pcriterion = 0.01;

    rho_input( pval_input > pcriterion ) = 0;
    
%% find NaNs

    kk = find(~any(rho_input,2));
    
    rho_input(kk,:) = [];
    rho_input(:,kk) = [];

%% compute the clusters

    num_Clust = 10;

    [Tidx,C,sumd,D] = kmeans(rho_input, num_Clust , 'distance', 'correlation', 'display', 'final','replicate',20);

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

function plotFCClusteredRegionByRegion(settings,idx_ROI)

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(label_ROI{iROI});
    
    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area_label,'.mat'));

    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area_label,'-Kmeans','.mat'));

    Ncomponent = size(FC_Voxels.run(1).rho_MOT4,1);
    cluster_assignment = FC_Kmeans.run(1).high_IdxClusters;
    Ncluster = FC_Kmeans.run(1).high_Nclusters;
    
    plotFCClustered(FC_Voxels.run(1).rho_MOT4,cluster_assignment, Ncluster, Ncomponent,strcat(area_label,'-','High','-','Run','-',int2str(1),'-FC-Clustered'));
    
    Ncomponent = size(FC_Voxels.run(2).rho_MOT4,1);
    cluster_assignment = FC_Kmeans.run(2).high_IdxClusters;
    Ncluster = FC_Kmeans.run(2).high_Nclusters;
    
    plotFCClustered(FC_Voxels.run(2).rho_MOT4,cluster_assignment, Ncluster, Ncomponent,strcat(area_label,'-','High','-','Run','-',int2str(2),'-FC-Clustered'));
    
    Ncomponent = size(FC_Voxels.run(1).rho_RestingState,1);
    cluster_assignment = FC_Kmeans.run(1).rest_IdxClusters;
    Ncluster = FC_Kmeans.run(1).rest_Nclusters;
    
    plotFCClustered(FC_Voxels.run(1).rho_RestingState,cluster_assignment, Ncluster, Ncomponent,strcat(area_label,'-','Rest','-','Run','-',int2str(1),'-FC-Clustered'));
    
    Ncomponent = size(FC_Voxels.run(2).rho_RestingState,1);
    cluster_assignment = FC_Kmeans.run(2).rest_IdxClusters;
    Ncluster = FC_Kmeans.run(2).rest_Nclusters;
    
    plotFCClustered(FC_Voxels.run(2).rho_RestingState,cluster_assignment, Ncluster, Ncomponent,strcat(area_label,'-','Rest','-','Run','-',int2str(2),'-FC-Clustered'));
    
    clear FC_Kmeans
    clear FC_Voxels
    
end

end

function plotFCClusteredRegionByRegion_half(settings,idx_ROI)

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
%load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(label_ROI{iROI});
    
    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Half-',area_label,'.mat'));

    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Half-',area_label,'-Kmeans','.mat'));

    
    %% FH
    
    Ncomponent = size(FC_Voxels.run(1).rho_MOT4_FH,1);
    cluster_assignment = FC_Kmeans.run(1).high_IdxClusters_FH;
    Ncluster = FC_Kmeans.run(1).high_Nclusters_FH;
    
    plotFCClustered(settings,FC_Voxels.run(1).rho_MOT4_FH,cluster_assignment, Ncluster, Ncomponent,strcat(area_label,'-','High','-','Run','-',int2str(1),'-FC-Clustered-FH'));
    
    Ncomponent = size(FC_Voxels.run(2).rho_MOT4_FH,1);
    cluster_assignment = FC_Kmeans.run(2).high_IdxClusters_FH;
    Ncluster = FC_Kmeans.run(2).high_Nclusters_FH;
    
    plotFCClustered(settings,FC_Voxels.run(2).rho_MOT4_FH,cluster_assignment, Ncluster, Ncomponent,strcat(area_label,'-','High','-','Run','-',int2str(2),'-FC-Clustered-FH'));
    
    Ncomponent = size(FC_Voxels.run(1).rho_RestingState_FH,1);
    cluster_assignment = FC_Kmeans.run(1).rest_IdxClusters_FH;
    Ncluster = FC_Kmeans.run(1).rest_Nclusters_FH;
    
    plotFCClustered(settings,FC_Voxels.run(1).rho_RestingState_FH,cluster_assignment, Ncluster, Ncomponent,strcat(area_label,'-','Rest','-','Run','-',int2str(1),'-FC-Clustered-FH'));
    
    Ncomponent = size(FC_Voxels.run(2).rho_RestingState_FH,1);
    cluster_assignment = FC_Kmeans.run(2).rest_IdxClusters_FH;
    Ncluster = FC_Kmeans.run(2).rest_Nclusters_FH;
    
    plotFCClustered(settings,FC_Voxels.run(2).rho_RestingState_FH,cluster_assignment, Ncluster, Ncomponent,strcat(area_label,'-','Rest','-','Run','-',int2str(2),'-FC-Clustered-FH'));
    
    
    %%% SH
    
    Ncomponent = size(FC_Voxels.run(1).rho_MOT4_SH,1);
    cluster_assignment = FC_Kmeans.run(1).high_IdxClusters_SH;
    Ncluster = FC_Kmeans.run(1).high_Nclusters_SH;
    
    plotFCClustered(settings,FC_Voxels.run(1).rho_MOT4_SH,cluster_assignment, Ncluster, Ncomponent,strcat(area_label,'-','High','-','Run','-',int2str(1),'-FC-Clustered-SH'));
    
    Ncomponent = size(FC_Voxels.run(2).rho_MOT4_SH,1);
    cluster_assignment = FC_Kmeans.run(2).high_IdxClusters_SH;
    Ncluster = FC_Kmeans.run(2).high_Nclusters_SH;
    
    plotFCClustered(settings,FC_Voxels.run(2).rho_MOT4_SH,cluster_assignment, Ncluster, Ncomponent,strcat(area_label,'-','High','-','Run','-',int2str(2),'-FC-Clustered-SH'));
    
    Ncomponent = size(FC_Voxels.run(1).rho_RestingState_SH,1);
    cluster_assignment = FC_Kmeans.run(1).rest_IdxClusters_SH;
    Ncluster = FC_Kmeans.run(1).rest_Nclusters_SH;
    
    plotFCClustered(settings,FC_Voxels.run(1).rho_RestingState_SH,cluster_assignment, Ncluster, Ncomponent,strcat(area_label,'-','Rest','-','Run','-',int2str(1),'-FC-Clustered-SH'));
    
    Ncomponent = size(FC_Voxels.run(2).rho_RestingState_SH,1);
    cluster_assignment = FC_Kmeans.run(2).rest_IdxClusters_SH;
    Ncluster = FC_Kmeans.run(2).rest_Nclusters_SH;
    
    plotFCClustered(settings,FC_Voxels.run(2).rho_RestingState_SH,cluster_assignment, Ncluster, Ncomponent,strcat(area_label,'-','Rest','-','Run','-',int2str(2),'-FC-Clustered-SH'));
    
    
    
    clear FC_Kmeans
    clear FC_Voxels
    
end

end

function plotFCClusteredRegionByRegion_half_SameOrder(settings,idx_ROI)

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
%load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(label_ROI{iROI});
    
    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area_label,'.mat'));

    high_run1 = FC_Voxels.run(1).rho_MOT4;
    high_run2 = FC_Voxels.run(2).rho_MOT4;
    
    rest_run1 = FC_Voxels.run(1).rho_RestingState;
    rest_run2 = FC_Voxels.run(2).rho_RestingState;
    
    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area_label,'-Kmeans','.mat'));

    high_run1_idx_clusters = FC_Kmeans.run(1).high_IdxClusters;
    high_run1_Nclusters = FC_Kmeans.run(1).high_Nclusters;
    
    high_run2_idx_clusters = FC_Kmeans.run(2).high_IdxClusters;
    high_run2_Nclusters = FC_Kmeans.run(2).high_Nclusters;
    
    rest_run1_idx_clusters = FC_Kmeans.run(1).rest_IdxClusters;
    rest_run1_Nclusters = FC_Kmeans.run(1).rest_Nclusters;
    
    rest_run2_idx_clusters = FC_Kmeans.run(2).rest_IdxClusters;
    rest_run2_Nclusters = FC_Kmeans.run(2).rest_Nclusters;
    
    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Half-',area_label,'.mat'));

    high_run1_FH = FC_Voxels.run(1).rho_MOT4_FH;
    high_run2_FH = FC_Voxels.run(2).rho_MOT4_FH;
    
    rest_run1_FH = FC_Voxels.run(1).rho_RestingState_FH;
    rest_run2_FH = FC_Voxels.run(2).rho_RestingState_FH;
    
    high_run1_SH = FC_Voxels.run(1).rho_MOT4_SH;
    high_run2_SH = FC_Voxels.run(2).rho_MOT4_SH;
    
    rest_run1_SH = FC_Voxels.run(1).rho_RestingState_SH;
    rest_run2_SH = FC_Voxels.run(2).rho_RestingState_SH;
    
    %load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Half-',area_label,'-Kmeans','.mat'));

    Ncomponent = size(high_run1,1);
    cluster_assignment = high_run1_idx_clusters;
    Ncluster = high_run1_Nclusters;
    
    plotFCClusteredFullHalfSameOrder(settings,high_run1,high_run1_FH,high_run1_SH,cluster_assignment, Ncluster, Ncomponent, strcat(area_label,'-','High','-','Run','-',int2str(1),'-FC-Clustered'));
    
    Ncomponent = size(high_run2,1);
    cluster_assignment = high_run2_idx_clusters;
    Ncluster = high_run2_Nclusters;
    
    plotFCClusteredFullHalfSameOrder(settings,high_run2,high_run2_FH,high_run2_SH,cluster_assignment, Ncluster, Ncomponent, strcat(area_label,'-','High','-','Run','-',int2str(2),'-FC-Clustered'));
    
    Ncomponent = size(rest_run1,1);
    cluster_assignment = rest_run1_idx_clusters;
    Ncluster = rest_run1_Nclusters;
    
    plotFCClusteredFullHalfSameOrder(settings,rest_run1,rest_run1_FH,rest_run1_SH,cluster_assignment, Ncluster, Ncomponent, strcat(area_label,'-','Rest','-','Run','-',int2str(1),'-FC-Clustered'));
    
    Ncomponent = size(rest_run2,1);
    cluster_assignment = rest_run2_idx_clusters;
    Ncluster = rest_run2_Nclusters;
    
    plotFCClusteredFullHalfSameOrder(settings,rest_run2,rest_run2_FH,rest_run2_SH,cluster_assignment, Ncluster, Ncomponent, strcat(area_label,'-','Rest','-','Run','-',int2str(2),'-FC-Clustered'));
    
    clear FC_Kmeans
    clear FC_Voxels
    
end

end

function plotFCClustered(settings,rho_mat,cluster_assignment, Ncluster, Ncomponent, plot_label)

    source_order = 1:Ncomponent;
    [target_order, cluster_idx_sorted] = TargetOrder( cluster_assignment, Ncluster );

    [rho_sorted, ~] = SortMatrix( rho_mat, source_order, target_order, Ncomponent );
    
    rlimit = 0.5;

    rmin = -rlimit + rlimit/32;
    rmax =  rlimit;

    Tick = [1 Ncomponent];
    TickLabel = { '1' ; num2str(Ncomponent) };

    f = figure;

    clrmp = colormap('jet');
    clrmp(32,:) = [1 1 1];

    hold on;
    caxis([rmin rmax]);
    h = pcolor( 0.5:Ncomponent-0.5, 0.5:Ncomponent-0.5, rho_sorted );
    set( h, 'EdgeColor', 'none');
    colormap( clrmp);
    hold off;
    axis 'square'; 
    axis([0 Ncomponent 0 Ncomponent]);
    xlabel(plot_label);
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

    h = colorbar;
    set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1]);

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'.jpg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'.pdf'));
    
end

function plotFCClusteredFullHalfSameOrder(settings,rho_mat_full,rho_mat_fh,rho_mat_sh,cluster_assignment, Ncluster, Ncomponent, plot_label)

    source_order = 1:Ncomponent;
    [target_order, cluster_idx_sorted] = TargetOrder( cluster_assignment, Ncluster );

    [rho_sorted_full, ~] = SortMatrix( rho_mat_full, source_order, target_order, Ncomponent );
    [rho_sorted_fh, ~] = SortMatrix( rho_mat_fh, source_order, target_order, Ncomponent );
    [rho_sorted_sh, ~] = SortMatrix( rho_mat_sh, source_order, target_order, Ncomponent );
    
    rlimit = 0.5;

    rmin = -rlimit + rlimit/32;
    rmax =  rlimit;

    Tick = [1 Ncomponent];
    TickLabel = { '1' ; num2str(Ncomponent) };

    f = figure;

    clrmp = colormap('jet');
    clrmp(32,:) = [1 1 1];

    hold on;
    caxis([rmin rmax]);
    h = pcolor( 0.5:Ncomponent-0.5, 0.5:Ncomponent-0.5, rho_sorted_full );
    set( h, 'EdgeColor', 'none');
    colormap( clrmp);
    hold off;
    axis 'square'; 
    axis([0 Ncomponent 0 Ncomponent]);
    xlabel(plot_label);
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );
    
    title('Full');

    h = colorbar;
    set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1]);

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'-Full-SameOrder','.jpg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'-Full-SameOrder','.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'-Full-SameOrder','.pdf'));
    
    f = figure;

    clrmp = colormap('jet');
    clrmp(32,:) = [1 1 1];

    hold on;
    caxis([rmin rmax]);
    h = pcolor( 0.5:Ncomponent-0.5, 0.5:Ncomponent-0.5, rho_sorted_fh );
    set( h, 'EdgeColor', 'none');
    colormap( clrmp);
    hold off;
    axis 'square'; 
    axis([0 Ncomponent 0 Ncomponent]);
    xlabel(plot_label);
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );
    
    title('FirstHalf');

    h = colorbar;
    set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1]);

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'-FH-SameOrder','.jpg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'-FH-SameOrder','.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'-FH-SameOrder','.pdf'));
    
    f = figure;

    clrmp = colormap('jet');
    clrmp(32,:) = [1 1 1];

    hold on;
    caxis([rmin rmax]);
    h = pcolor( 0.5:Ncomponent-0.5, 0.5:Ncomponent-0.5, rho_sorted_sh );
    set( h, 'EdgeColor', 'none');
    colormap( clrmp);
    hold off;
    axis 'square'; 
    axis([0 Ncomponent 0 Ncomponent]);
    xlabel(plot_label);
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );
    
    title('SecondHalf');

    h = colorbar;
    set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1]);

    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'-SH-SameOrder','.jpg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'-SH-SameOrder','.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'-SH-SameOrder','.pdf'));
    
end

function computeFCMeanSQClustersAllROIs(settings,idx_ROI)

nROI = length(idx_ROI);
nClusters = 10;
nTR = 300;

%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

lowhigh_load_all_data_FSL;

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

all_ROI_means_MOT4Run1 = zeros(nROI*nClusters,nTR);
all_ROI_means_MOT4Run2 = zeros(nROI*nClusters,nTR);
all_ROI_means_RestingStateRun1 = zeros(nROI*nClusters,nTR);
all_ROI_means_RestingStateRun2 = zeros(nROI*nClusters,nTR);

iROICluster = 0;

for iROI=1:nROI
    
    idx_AAL = AAL_ROI(idx_ROI(iROI)).ID;
   
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(label_ROI{iROI});
    
    idx_voxels = find(AAL_img==idx_AAL);
    
    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area_label,'-Kmeans','.mat'));

    for iCluster=1:nClusters
        
        iROICluster = iROICluster + 1;
        
        idx_voxels_on_cluster = find(FC_Kmeans.run(1).high_IdxClusters==iCluster);
        selected_voxels = idx_voxels(idx_voxels_on_cluster);
        nVoxels = length(selected_voxels);
        run = zeros(nVoxels,nTR);
        for iVoxel=1:nVoxels
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),selected_voxels(iVoxel));
            run(iVoxel,:) = MOT4Run1(idxx,idxy,idxz,:);
        end
        
        FCMean_n.run(1).ROI(iROI).Cluster(iCluster).high = nVoxels;
        all_ROI_means_MOT4Run1(iROICluster,:) = mean(run,1);
        
        clear run
        clear selected_voxels
        
        idx_voxels_on_cluster = find(FC_Kmeans.run(2).high_IdxClusters==iCluster);
        selected_voxels = idx_voxels(idx_voxels_on_cluster);
        nVoxels = length(selected_voxels);
        run = zeros(nVoxels,nTR);
        for iVoxel=1:nVoxels
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),selected_voxels(iVoxel));
            run(iVoxel,:) = MOT4Run2(idxx,idxy,idxz,:);
        end
        
        FCMean_n.run(2).ROI(iROI).Cluster(iCluster).high = nVoxels;
        all_ROI_means_MOT4Run2(iROICluster,:) = mean(run,1);
        
        clear run
        clear selected_voxels
        
        idx_voxels_on_cluster = find(FC_Kmeans.run(1).rest_IdxClusters==iCluster);
        selected_voxels = idx_voxels(idx_voxels_on_cluster);
        nVoxels = length(selected_voxels);
        run = zeros(nVoxels,nTR);
        for iVoxel=1:nVoxels
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),selected_voxels(iVoxel));
            run(iVoxel,:) = RestingStateRun1(idxx,idxy,idxz,:);
        end
        
        FCMean_n.run(1).ROI(iROI).Cluster(iCluster).rest = nVoxels;
        all_ROI_means_RestingStateRun1(iROICluster,:) = mean(run,1);
        
        clear run
        clear selected_voxels
        
        idx_voxels_on_cluster = find(FC_Kmeans.run(2).rest_IdxClusters==iCluster);
        selected_voxels = idx_voxels(idx_voxels_on_cluster);
        nVoxels = length(selected_voxels);
        run = zeros(nVoxels,nTR);
        for iVoxel=1:nVoxels
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),selected_voxels(iVoxel));
            run(iVoxel,:) = RestingStateRun2(idxx,idxy,idxz,:);
        end
        
        FCMean_n.run(2).ROI(iROI).Cluster(iCluster).rest = nVoxels;
        all_ROI_means_RestingStateRun2(iROICluster,:) = mean(run,1);
        
        clear run
        clear selected_voxels
        
    end
    
    clear FC_Kmeans
    
end

all_ROI_means_MOT4Run1 = all_ROI_means_MOT4Run1';
all_ROI_means_MOT4Run2 = all_ROI_means_MOT4Run2';

all_ROI_means_RestingStateRun1 = all_ROI_means_RestingStateRun1';
all_ROI_means_RestingStateRun2 = all_ROI_means_RestingStateRun2';

disp('doing Pearson');

[FCMean_rp.run(1).rho_high,FCMean_rp.run(1).pval_high] = corr(all_ROI_means_MOT4Run1);
[FCMean_rp.run(2).rho_high,FCMean_rp.run(2).pval_high] = corr(all_ROI_means_MOT4Run2);

[FCMean_rp.run(1).rho_rest,FCMean_rp.run(1).pval_rest] = corr(all_ROI_means_RestingStateRun1);
[FCMean_rp.run(2).rho_rest,FCMean_rp.run(2).pval_rest] = corr(all_ROI_means_RestingStateRun2);

disp('doing Sum of Square');

sq_MOT4Run1 = zeros(nROI,nROI);
sq_MOT4Run2 = zeros(nROI,nROI);
sq_RestingStateRun1 = zeros(nROI,nROI);
sq_RestingStateRun2 = zeros(nROI,nROI);

for iROI=1:nROI
    
   for iiROI=1:nROI
       
      sq_clusters_mot4run1 = zeros(nClusters,nClusters);
      sq_clusters_mot4run2 = zeros(nClusters,nClusters);
      sq_clusters_restingstaterun1 = zeros(nClusters,nClusters);
      sq_clusters_restingstaterun2 = zeros(nClusters,nClusters);
      for iCluster=1:nClusters
          
          for iiCluster=1:nClusters
              
                sq_clusters_mot4run1(iCluster,iiCluster) = FCMean_rp.run(1).rho_high(((iROI-1)*nClusters + iCluster),((iiROI-1)*nClusters + iiCluster)) * FCMean_n.run(1).ROI(iROI).Cluster(iCluster).high * FCMean_n.run(1).ROI(iiROI).Cluster(iiCluster).high;
                sq_clusters_mot4run2(iCluster,iiCluster) = FCMean_rp.run(2).rho_high(((iROI-1)*nClusters + iCluster),((iiROI-1)*nClusters + iiCluster)) * FCMean_n.run(2).ROI(iROI).Cluster(iCluster).high * FCMean_n.run(2).ROI(iiROI).Cluster(iiCluster).high;
 
                sq_clusters_restingstaterun1(iCluster,iiCluster) = FCMean_rp.run(1).rho_rest(((iROI-1)*nClusters + iCluster),((iiROI-1)*nClusters + iiCluster)) * FCMean_n.run(1).ROI(iROI).Cluster(iCluster).rest * FCMean_n.run(1).ROI(iiROI).Cluster(iiCluster).rest;
                sq_clusters_restingstaterun2(iCluster,iiCluster) = FCMean_rp.run(2).rho_rest(((iROI-1)*nClusters + iCluster),((iiROI-1)*nClusters + iiCluster)) * FCMean_n.run(2).ROI(iROI).Cluster(iCluster).rest * FCMean_n.run(2).ROI(iiROI).Cluster(iiCluster).rest;
 
          end
          
      end
    
      sq_MOT4Run1(iROI,iiROI) = sumsqr(sq_clusters_mot4run1);
      sq_MOT4Run2(iROI,iiROI) = sumsqr(sq_clusters_mot4run2);
      sq_RestingStateRun1(iROI,iiROI) = sumsqr(sq_clusters_restingstaterun1);
      sq_RestingStateRun2(iROI,iiROI) = sumsqr(sq_clusters_restingstaterun2);

   end
    
end

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Mean-Clusters-SQ','.mat'),'all_ROI_means_MOT4Run1','all_ROI_means_MOT4Run2','all_ROI_means_RestingStateRun1','all_ROI_means_RestingStateRun2','FCMean_n','FCMean_rp','sq_MOT4Run1','sq_MOT4Run2','sq_RestingStateRun1','sq_RestingStateRun2');

end

function computeFCMeanSQClustersAllROIs_half(settings,idx_ROI)

nROI = length(idx_ROI);
nClusters = 10;
nTR = 300;

%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

lowhigh_load_all_data_FSL;

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
%load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

all_ROI_means_MOT4Run1_FH = zeros(nROI*nClusters,nTR/2);
all_ROI_means_MOT4Run2_FH = zeros(nROI*nClusters,nTR/2);
all_ROI_means_RestingStateRun1_FH = zeros(nROI*nClusters,nTR/2);
all_ROI_means_RestingStateRun2_FH = zeros(nROI*nClusters,nTR/2);

all_ROI_means_MOT4Run1_SH = zeros(nROI*nClusters,nTR/2);
all_ROI_means_MOT4Run2_SH = zeros(nROI*nClusters,nTR/2);
all_ROI_means_RestingStateRun1_SH = zeros(nROI*nClusters,nTR/2);
all_ROI_means_RestingStateRun2_SH = zeros(nROI*nClusters,nTR/2);

iROICluster = 0;

for iROI=1:nROI
    
    idx_AAL = AAL_ROI(idx_ROI(iROI)).ID;
   
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(label_ROI{iROI});
    
    idx_voxels = find(AAL_img==idx_AAL);
    
    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Half-',area_label,'-Kmeans','.mat'));

    for iCluster=1:nClusters
        
        iROICluster = iROICluster + 1;
        
        idx_voxels_on_cluster = find(FC_Kmeans.run(1).high_IdxClusters_FH==iCluster);
        selected_voxels = idx_voxels(idx_voxels_on_cluster);
        nVoxels = length(selected_voxels);
        run = zeros(nVoxels,nTR/2);
        for iVoxel=1:nVoxels
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),selected_voxels(iVoxel));
            run(iVoxel,:) = MOT4Run1(idxx,idxy,idxz,1:(nTR/2));
        end
        
        FCMean_n.run(1).ROI(iROI).Cluster(iCluster).high_FH = nVoxels;
        all_ROI_means_MOT4Run1_FH(iROICluster,:) = mean(run,1);
        
        clear run
        clear selected_voxels
        
        idx_voxels_on_cluster = find(FC_Kmeans.run(2).high_IdxClusters_FH==iCluster);
        selected_voxels = idx_voxels(idx_voxels_on_cluster);
        nVoxels = length(selected_voxels);
        run = zeros(nVoxels,nTR/2);
        for iVoxel=1:nVoxels
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),selected_voxels(iVoxel));
            run(iVoxel,:) = MOT4Run2(idxx,idxy,idxz,1:(nTR/2));
        end
        
        FCMean_n.run(2).ROI(iROI).Cluster(iCluster).high_FH = nVoxels;
        all_ROI_means_MOT4Run2_FH(iROICluster,:) = mean(run,1);
        
        clear run
        clear selected_voxels
        
        idx_voxels_on_cluster = find(FC_Kmeans.run(1).rest_IdxClusters_FH==iCluster);
        selected_voxels = idx_voxels(idx_voxels_on_cluster);
        nVoxels = length(selected_voxels);
        run = zeros(nVoxels,nTR/2);
        for iVoxel=1:nVoxels
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),selected_voxels(iVoxel));
            run(iVoxel,:) = RestingStateRun1(idxx,idxy,idxz,1:(nTR/2));
        end
        
        FCMean_n.run(1).ROI(iROI).Cluster(iCluster).rest_FH = nVoxels;
        all_ROI_means_RestingStateRun1_FH(iROICluster,:) = mean(run,1);
        
        clear run
        clear selected_voxels
        
        idx_voxels_on_cluster = find(FC_Kmeans.run(2).rest_IdxClusters_FH==iCluster);
        selected_voxels = idx_voxels(idx_voxels_on_cluster);
        nVoxels = length(selected_voxels);
        run = zeros(nVoxels,nTR/2);
        for iVoxel=1:nVoxels
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),selected_voxels(iVoxel));
            run(iVoxel,:) = RestingStateRun2(idxx,idxy,idxz,1:(nTR/2));
        end
        
        FCMean_n.run(2).ROI(iROI).Cluster(iCluster).rest_FH = nVoxels;
        all_ROI_means_RestingStateRun2_FH(iROICluster,:) = mean(run,1);
        
        clear run
        clear selected_voxels
        
    end
    
    for iCluster=1:nClusters
        
        iROICluster = iROICluster + 1;
        
        idx_voxels_on_cluster = find(FC_Kmeans.run(1).high_IdxClusters_SH==iCluster);
        selected_voxels = idx_voxels(idx_voxels_on_cluster);
        nVoxels = length(selected_voxels);
        run = zeros(nVoxels,nTR/2);
        for iVoxel=1:nVoxels
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),selected_voxels(iVoxel));
            run(iVoxel,:) = MOT4Run1(idxx,idxy,idxz,((nTR/2) + 1):nTR);
        end
        
        FCMean_n.run(1).ROI(iROI).Cluster(iCluster).high_SH = nVoxels;
        all_ROI_means_MOT4Run1_SH(iROICluster,:) = mean(run,1);
        
        clear run
        clear selected_voxels
        
        idx_voxels_on_cluster = find(FC_Kmeans.run(2).high_IdxClusters_SH==iCluster);
        selected_voxels = idx_voxels(idx_voxels_on_cluster);
        nVoxels = length(selected_voxels);
        run = zeros(nVoxels,nTR/2);
        for iVoxel=1:nVoxels
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),selected_voxels(iVoxel));
            run(iVoxel,:) = MOT4Run2(idxx,idxy,idxz,((nTR/2) + 1):nTR);
        end
        
        FCMean_n.run(2).ROI(iROI).Cluster(iCluster).high_SH = nVoxels;
        all_ROI_means_MOT4Run2_SH(iROICluster,:) = mean(run,1);
        
        clear run
        clear selected_voxels
        
        idx_voxels_on_cluster = find(FC_Kmeans.run(1).rest_IdxClusters_SH==iCluster);
        selected_voxels = idx_voxels(idx_voxels_on_cluster);
        nVoxels = length(selected_voxels);
        run = zeros(nVoxels,nTR/2);
        for iVoxel=1:nVoxels
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),selected_voxels(iVoxel));
            run(iVoxel,:) = RestingStateRun1(idxx,idxy,idxz,((nTR/2) + 1):nTR);
        end
        
        FCMean_n.run(1).ROI(iROI).Cluster(iCluster).rest_SH = nVoxels;
        all_ROI_means_RestingStateRun1_SH(iROICluster,:) = mean(run,1);
        
        clear run
        clear selected_voxels
        
        idx_voxels_on_cluster = find(FC_Kmeans.run(2).rest_IdxClusters_SH==iCluster);
        selected_voxels = idx_voxels(idx_voxels_on_cluster);
        nVoxels = length(selected_voxels);
        run = zeros(nVoxels,nTR/2);
        for iVoxel=1:nVoxels
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),selected_voxels(iVoxel));
            run(iVoxel,:) = RestingStateRun2(idxx,idxy,idxz,((nTR/2) + 1):nTR);
        end
        
        FCMean_n.run(2).ROI(iROI).Cluster(iCluster).rest_SH = nVoxels;
        all_ROI_means_RestingStateRun2_SH(iROICluster,:) = mean(run,1);
        
        clear run
        clear selected_voxels
        
    end
    
    clear FC_Kmeans
    
end

all_ROI_means_MOT4Run1_FH = all_ROI_means_MOT4Run1_FH';
all_ROI_means_MOT4Run2_FH = all_ROI_means_MOT4Run2_FH';

all_ROI_means_RestingStateRun1_FH = all_ROI_means_RestingStateRun1_FH';
all_ROI_means_RestingStateRun2_FH = all_ROI_means_RestingStateRun2_FH';

all_ROI_means_MOT4Run1_SH = all_ROI_means_MOT4Run1_SH';
all_ROI_means_MOT4Run2_SH = all_ROI_means_MOT4Run2_SH';

all_ROI_means_RestingStateRun1_SH = all_ROI_means_RestingStateRun1_SH';
all_ROI_means_RestingStateRun2_SH = all_ROI_means_RestingStateRun2_SH';

disp('doing Pearson');

[FCMean_rp.run(1).rho_high_FH,FCMean_rp.run(1).pval_high_FH] = corr(all_ROI_means_MOT4Run1_FH);
[FCMean_rp.run(2).rho_high_FH,FCMean_rp.run(2).pval_high_FH] = corr(all_ROI_means_MOT4Run2_FH);

[FCMean_rp.run(1).rho_rest_FH,FCMean_rp.run(1).pval_rest_FH] = corr(all_ROI_means_RestingStateRun1_FH);
[FCMean_rp.run(2).rho_rest_FH,FCMean_rp.run(2).pval_rest_FH] = corr(all_ROI_means_RestingStateRun2_FH);

[FCMean_rp.run(1).rho_high_SH,FCMean_rp.run(1).pval_high_SH] = corr(all_ROI_means_MOT4Run1_SH);
[FCMean_rp.run(2).rho_high_SH,FCMean_rp.run(2).pval_high_SH] = corr(all_ROI_means_MOT4Run2_SH);

[FCMean_rp.run(1).rho_rest_SH,FCMean_rp.run(1).pval_rest_SH] = corr(all_ROI_means_RestingStateRun1_SH);
[FCMean_rp.run(2).rho_rest_SH,FCMean_rp.run(2).pval_rest_SH] = corr(all_ROI_means_RestingStateRun2_SH);

disp('doing Sum of Square');

sq_MOT4Run1_FH = zeros(nROI,nROI);
sq_MOT4Run2_FH = zeros(nROI,nROI);
sq_RestingStateRun1_FH = zeros(nROI,nROI);
sq_RestingStateRun2_FH = zeros(nROI,nROI);

sq_MOT4Run1_SH = zeros(nROI,nROI);
sq_MOT4Run2_SH = zeros(nROI,nROI);
sq_RestingStateRun1_SH = zeros(nROI,nROI);
sq_RestingStateRun2_SH = zeros(nROI,nROI);

for iROI=1:nROI
    
   for iiROI=1:nROI
       
      sq_clusters_mot4run1_FH = zeros(nClusters,nClusters);
      sq_clusters_mot4run2_FH = zeros(nClusters,nClusters);
      sq_clusters_restingstaterun1_FH = zeros(nClusters,nClusters);
      sq_clusters_restingstaterun2_FH = zeros(nClusters,nClusters);
      
      sq_clusters_mot4run1_SH = zeros(nClusters,nClusters);
      sq_clusters_mot4run2_SH = zeros(nClusters,nClusters);
      sq_clusters_restingstaterun1_SH = zeros(nClusters,nClusters);
      sq_clusters_restingstaterun2_SH = zeros(nClusters,nClusters);
      
      for iCluster=1:nClusters
          
          for iiCluster=1:nClusters
              
                sq_clusters_mot4run1_FH(iCluster,iiCluster) = FCMean_rp.run(1).rho_high_FH(((iROI-1)*nClusters + iCluster),((iiROI-1)*nClusters + iiCluster)) * FCMean_n.run(1).ROI(iROI).Cluster(iCluster).high_FH * FCMean_n.run(1).ROI(iiROI).Cluster(iiCluster).high_FH;
                sq_clusters_mot4run2_FH(iCluster,iiCluster) = FCMean_rp.run(2).rho_high_FH(((iROI-1)*nClusters + iCluster),((iiROI-1)*nClusters + iiCluster)) * FCMean_n.run(2).ROI(iROI).Cluster(iCluster).high_FH * FCMean_n.run(2).ROI(iiROI).Cluster(iiCluster).high_FH;
 
                sq_clusters_restingstaterun1_FH(iCluster,iiCluster) = FCMean_rp.run(1).rho_rest_FH(((iROI-1)*nClusters + iCluster),((iiROI-1)*nClusters + iiCluster)) * FCMean_n.run(1).ROI(iROI).Cluster(iCluster).rest_FH * FCMean_n.run(1).ROI(iiROI).Cluster(iiCluster).rest_FH;
                sq_clusters_restingstaterun2_FH(iCluster,iiCluster) = FCMean_rp.run(2).rho_rest_FH(((iROI-1)*nClusters + iCluster),((iiROI-1)*nClusters + iiCluster)) * FCMean_n.run(2).ROI(iROI).Cluster(iCluster).rest_FH * FCMean_n.run(2).ROI(iiROI).Cluster(iiCluster).rest_FH;
 
                sq_clusters_mot4run1_SH(iCluster,iiCluster) = FCMean_rp.run(1).rho_high_SH(((iROI-1)*nClusters + iCluster),((iiROI-1)*nClusters + iiCluster)) * FCMean_n.run(1).ROI(iROI).Cluster(iCluster).high_SH * FCMean_n.run(1).ROI(iiROI).Cluster(iiCluster).high_SH;
                sq_clusters_mot4run2_SH(iCluster,iiCluster) = FCMean_rp.run(2).rho_high_SH(((iROI-1)*nClusters + iCluster),((iiROI-1)*nClusters + iiCluster)) * FCMean_n.run(2).ROI(iROI).Cluster(iCluster).high_SH * FCMean_n.run(2).ROI(iiROI).Cluster(iiCluster).high_SH;
 
                sq_clusters_restingstaterun1_SH(iCluster,iiCluster) = FCMean_rp.run(1).rho_rest_SH(((iROI-1)*nClusters + iCluster),((iiROI-1)*nClusters + iiCluster)) * FCMean_n.run(1).ROI(iROI).Cluster(iCluster).rest_SH * FCMean_n.run(1).ROI(iiROI).Cluster(iiCluster).rest_SH;
                sq_clusters_restingstaterun2_SH(iCluster,iiCluster) = FCMean_rp.run(2).rho_rest_SH(((iROI-1)*nClusters + iCluster),((iiROI-1)*nClusters + iiCluster)) * FCMean_n.run(2).ROI(iROI).Cluster(iCluster).rest_SH * FCMean_n.run(2).ROI(iiROI).Cluster(iiCluster).rest_SH;
 
          end
          
      end
    
      sq_MOT4Run1_FH(iROI,iiROI) = sumsqr(sq_clusters_mot4run1_FH);
      sq_MOT4Run2_FH(iROI,iiROI) = sumsqr(sq_clusters_mot4run2_FH);
      sq_RestingStateRun1_FH(iROI,iiROI) = sumsqr(sq_clusters_restingstaterun1_FH);
      sq_RestingStateRun2_FH(iROI,iiROI) = sumsqr(sq_clusters_restingstaterun2_FH);
      
      sq_MOT4Run1_SH(iROI,iiROI) = sumsqr(sq_clusters_mot4run1_SH);
      sq_MOT4Run2_SH(iROI,iiROI) = sumsqr(sq_clusters_mot4run2_SH);
      sq_RestingStateRun1_SH(iROI,iiROI) = sumsqr(sq_clusters_restingstaterun1_SH);
      sq_RestingStateRun2_SH(iROI,iiROI) = sumsqr(sq_clusters_restingstaterun2_SH);

   end
    
end

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Mean-Clusters-SQ-Half-',area_label,'.mat'),...
    'all_ROI_means_MOT4Run1_FH','all_ROI_means_MOT4Run2_FH','all_ROI_means_RestingStateRun1_FH','all_ROI_means_RestingStateRun2_FH','sq_MOT4Run1_FH','sq_MOT4Run2_FH','sq_RestingStateRun1_FH','sq_RestingStateRun2_FH',...
    'all_ROI_means_MOT4Run1_SH','all_ROI_means_MOT4Run2_SH','all_ROI_means_RestingStateRun1_SH','all_ROI_means_RestingStateRun2_SH','sq_MOT4Run1_SH','sq_MOT4Run2_SH','sq_RestingStateRun1_SH','sq_RestingStateRun2_SH','FCMean_n','FCMean_rp');
end

function plotFCMeanSQClustersAllROIs(settings,idx_ROI)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Mean-Clusters-SQ','.mat'));

mean_val = mean([mean(sq_MOT4Run1(:)),mean(sq_MOT4Run2(:)),mean(sq_RestingStateRun1(:)),mean(sq_RestingStateRun2(:))]);

std_val = mean([std(sq_MOT4Run1(:)),std(sq_MOT4Run2(:)),std(sq_RestingStateRun1(:)),std(sq_RestingStateRun2(:))]);

%CLIM = [0 (mean_val + 2*std_val)];
CLIM = [0 (mean_val + std_val)];

plotSQ(settings,sq_MOT4Run1,'High-Run1',CLIM);
plotSQ(settings,sq_MOT4Run2,'High-Run2',CLIM);

disp(strcat('Pearson - High Attention - Runs 1 & 2:',num2str(corr(sq_MOT4Run1(:),sq_MOT4Run2(:)))));

plotSQ(settings,sq_RestingStateRun1,'RestingState-Run1',CLIM);
plotSQ(settings,sq_RestingStateRun2,'RestingState-Run2',CLIM);

disp(strcat('Pearson - Resting State - Runs 1 & 2:',num2str(corr(sq_RestingStateRun1(:),sq_RestingStateRun2(:)))));

end

function plotFCMeanSQClustersAllROIs_half(settings,idx_ROI)

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
%load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

iROI = 1;

idx_AAL = AAL_ROI(idx_ROI(iROI)).ID;
   
label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
area_label = strrep(label_ROI{iROI},'_','-');

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Mean-Clusters-SQ-Half-',area_label,'.mat'));

mean_val = mean([mean(sq_MOT4Run1(:)),mean(sq_MOT4Run2(:)),mean(sq_RestingStateRun1(:)),mean(sq_RestingStateRun2(:))]);

std_val = mean([std(sq_MOT4Run1(:)),std(sq_MOT4Run2(:)),std(sq_RestingStateRun1(:)),std(sq_RestingStateRun2(:))]);

%CLIM = [0 (mean_val + 2*std_val)];
CLIM = [0 (mean_val + std_val)];

plotSQ(settings,sq_MOT4Run1,'High-Run1',CLIM);
plotSQ(settings,sq_MOT4Run2,'High-Run2',CLIM);

disp(strcat('Pearson - High Attention - Runs 1 & 2:',num2str(corr(sq_MOT4Run1(:),sq_MOT4Run2(:)))));

plotSQ(settings,sq_RestingStateRun1,'RestingState-Run1',CLIM);
plotSQ(settings,sq_RestingStateRun2,'RestingState-Run2',CLIM);

disp(strcat('Pearson - Resting State - Runs 1 & 2:',num2str(corr(sq_RestingStateRun1(:),sq_RestingStateRun2(:)))));

end

function plotSQ(settings,rho,label,CLIM)

F_n = 22;
O_n = 12;
P_n = 18;

f = figure;

%CLIM = [min_aal max_aal];

imagesc(rho,CLIM);

title(label);

clrmp = colormap('jet');
clrmp(32,:) = [1 1 1];
colormap(clrmp);
h = colorbar;
set(h,'YLim',round(CLIM), 'YTick',round(CLIM), 'PlotBoxAspectRatio', [1 20 1]);

Tick = [F_n/2 (F_n+P_n/2) (F_n+P_n+O_n/2)];
TickLabel = { 'F' ; 'P' ; 'O' };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

hold on;

line_position = F_n;
plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = F_n + P_n;
plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Mean-Clusters-SQ','-',label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Mean-Clusters-SQ','-',label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Mean-Clusters-SQ','-',label,'.pdf'));


end

function computePearsonBtwROI(settings)

subject = settings.folders.subject;

%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

lowhigh_load_all_data_FSL;

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

% area1_ID = 7; %% FRONTAL MID L
% area2_ID = 52; %% OCCIPITAL MID R

area1_ID = 13; %% FRONTAL INF TRI L
area2_ID = 66; %% ANGULAR R

label_ROI_area1 = AAL_ROI(area1_ID).Nom_L;
label_ROI_area2 = AAL_ROI(area2_ID).Nom_L;

idx_AAL_area1 = AAL_ROI(area1_ID).ID;
idx_AAL_area2 = AAL_ROI(area2_ID).ID;

area1_label = strrep(label_ROI_area1,'_','-');
area2_label = strrep(label_ROI_area2,'_','-');

disp(strcat(label_ROI_area1,'-',label_ROI_area2));

disp('High Attention');
[FC_Voxels.run(1).rho_MOT4, FC_Voxels.run(1).pval_MOT4] = getVoxelCorrelations2Areas(MOT4Run1,AAL_img,idx_AAL_area1,idx_AAL_area2);
[FC_Voxels.run(2).rho_MOT4, FC_Voxels.run(2).pval_MOT4] = getVoxelCorrelations2Areas(MOT4Run2,AAL_img,idx_AAL_area1,idx_AAL_area2);

%     disp('Low Attention');
%     [FC_Voxels.run(1).rho_MOT2, FC_Voxels.run(1).pval_MOT2] = getCorrelations(MOT2Run1,AAL_img,AAL_ROI(idx_ROI(iROI)).ID);
%     [FC_Voxels.run(2).rho_MOT2, FC_Voxels.run(2).pval_MOT2] = getCorrelations(MOT2Run2,AAL_img,AAL_ROI(idx_ROI(iROI)).ID);

disp('Resting State');
[FC_Voxels.run(1).rho_RestingState, FC_Voxels.run(1).pval_RestingState] = getVoxelCorrelations2Areas(RestingStateRun1,AAL_img,idx_AAL_area1,idx_AAL_area2);
[FC_Voxels.run(2).rho_RestingState, FC_Voxels.run(2).pval_RestingState] = getVoxelCorrelations2Areas(RestingStateRun2,AAL_img,idx_AAL_area1,idx_AAL_area2);

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area1_label,'-',area2_label,'.mat'),'subject','file','area1_label','area2_label','FC_Voxels');

clear FC_Voxels



end

function plotFCClusterTwoRegions_all_rhos(settings,area1_ID,area2_ID)

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

label_ROI = AAL_ROI(area1_ID).Nom_L;
area1_label = strrep(label_ROI,'_','-');
disp(label_ROI);

label_ROI = AAL_ROI(area2_ID).Nom_L;
area2_label = strrep(label_ROI,'_','-');
disp(label_ROI);

area1.FC = load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area1_label,'.mat'));

area1.Kmeans = load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area1_label,'-Kmeans','.mat'));

area2.FC = load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area2_label,'.mat'));

area2.Kmeans = load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area2_label,'-Kmeans','.mat'));

area1_area2.FC = load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area1_label,'-',area2_label,'.mat'));

area1_rho = area1.FC.FC_Voxels.run(1).rho_MOT4;
area1_pval = area1.FC.FC_Voxels.run(1).pval_MOT4;
area1_IdxClusters = area1.Kmeans.FC_Kmeans.run(1).high_IdxClusters;
area1_Nclusters = area1.Kmeans.FC_Kmeans.run(1).high_Nclusters;
area1_label_run = strcat(area1_label,'-','High-Run-1');
area2_rho = area2.FC.FC_Voxels.run(1).rho_MOT4;
area2_pval = area2.FC.FC_Voxels.run(1).pval_MOT4;
area2_IdxClusters = area2.Kmeans.FC_Kmeans.run(1).high_IdxClusters;
area2_Nclusters = area2.Kmeans.FC_Kmeans.run(1).high_Nclusters;
area2_label_run = strcat(area2_label,'-','High-Run-1');
area1_area1_rho = area1_area2.FC.FC_Voxels.run(1).rho_MOT4;
area1_area2_pval = area1_area2.FC.FC_Voxels.run(1).rho_MOT4;
plotFCClusterTwoRegions_one_rho(settings,area1_rho,area1_pval,area1_IdxClusters,area1_Nclusters,area1_label_run,area2_rho,area2_pval,area2_IdxClusters,area2_Nclusters,area2_label_run,area1_area1_rho,area1_area2_pval);

area1_rho = area1.FC.FC_Voxels.run(2).rho_MOT4;
area1_pval = area1.FC.FC_Voxels.run(2).pval_MOT4;
area1_IdxClusters = area1.Kmeans.FC_Kmeans.run(2).high_IdxClusters;
area1_Nclusters = area1.Kmeans.FC_Kmeans.run(2).high_Nclusters;
area1_label_run = strcat(area1_label,'-','High-Run-2');
area2_rho = area2.FC.FC_Voxels.run(2).rho_MOT4;
area2_pval = area2.FC.FC_Voxels.run(2).pval_MOT4;
area2_IdxClusters = area2.Kmeans.FC_Kmeans.run(2).high_IdxClusters;
area2_Nclusters = area2.Kmeans.FC_Kmeans.run(2).high_Nclusters;
area2_label_run = strcat(area2_label,'-','High-Run-2');
area1_area1_rho = area1_area2.FC.FC_Voxels.run(2).rho_MOT4;
area1_area2_pval = area1_area2.FC.FC_Voxels.run(2).rho_MOT4;
plotFCClusterTwoRegions_one_rho(settings,area1_rho,area1_pval,area1_IdxClusters,area1_Nclusters,area1_label_run,area2_rho,area2_pval,area2_IdxClusters,area2_Nclusters,area2_label_run,area1_area1_rho,area1_area2_pval);

area1_rho = area1.FC.FC_Voxels.run(1).rho_RestingState;
area1_pval = area1.FC.FC_Voxels.run(1).pval_RestingState;
area1_IdxClusters = area1.Kmeans.FC_Kmeans.run(1).rest_IdxClusters;
area1_Nclusters = area1.Kmeans.FC_Kmeans.run(1).rest_Nclusters;
area1_label_run = strcat(area1_label,'-','Rest-Run-1');
area2_rho = area2.FC.FC_Voxels.run(1).rho_RestingState;
area2_pval = area2.FC.FC_Voxels.run(1).pval_RestingState;
area2_IdxClusters = area2.Kmeans.FC_Kmeans.run(1).rest_IdxClusters;
area2_Nclusters = area2.Kmeans.FC_Kmeans.run(1).rest_Nclusters;
area2_label_run = strcat(area2_label,'-','Rest-Run-1');
area1_area1_rho = area1_area2.FC.FC_Voxels.run(1).rho_RestingState;
area1_area2_pval = area1_area2.FC.FC_Voxels.run(1).rho_RestingState;
plotFCClusterTwoRegions_one_rho(settings,area1_rho,area1_pval,area1_IdxClusters,area1_Nclusters,area1_label_run,area2_rho,area2_pval,area2_IdxClusters,area2_Nclusters,area2_label_run,area1_area1_rho,area1_area2_pval);

area1_rho = area1.FC.FC_Voxels.run(2).rho_RestingState;
area1_pval = area1.FC.FC_Voxels.run(2).pval_RestingState;
area1_IdxClusters = area1.Kmeans.FC_Kmeans.run(2).rest_IdxClusters;
area1_Nclusters = area1.Kmeans.FC_Kmeans.run(2).rest_Nclusters;
area1_label_run = strcat(area1_label,'-','Rest-Run-2');
area2_rho = area2.FC.FC_Voxels.run(2).rho_RestingState;
area2_pval = area2.FC.FC_Voxels.run(2).pval_RestingState;
area2_IdxClusters = area2.Kmeans.FC_Kmeans.run(2).rest_IdxClusters;
area2_Nclusters = area2.Kmeans.FC_Kmeans.run(2).rest_Nclusters;
area2_label_run = strcat(area2_label,'-','Rest-Run-2');
area1_area1_rho = area1_area2.FC.FC_Voxels.run(2).rho_RestingState;
area1_area2_pval = area1_area2.FC.FC_Voxels.run(2).rho_RestingState;
plotFCClusterTwoRegions_one_rho(settings,area1_rho,area1_pval,area1_IdxClusters,area1_Nclusters,area1_label_run,area2_rho,area2_pval,area2_IdxClusters,area2_Nclusters,area2_label_run,area1_area1_rho,area1_area2_pval);


end

function plotFCClusterTwoRegions_one_rho(settings,area1_rho,area1_pval,area1_IdxClusters,area1_Nclusters,area1_label,area2_rho,area2_pval,area2_IdxClusters,area2_Nclusters,area2_label,area1_area2_rho,area1_area2_pval)

% prepare correlation plots
fs = 14;
pcriterion = 0.001;
rholimit   = 0.5;
rhomin = -rholimit + rholimit/32;
rhomax =  rholimit;
CLIM = [rhomin rhomax];
rhorange = linspace( rhomin, rhomax, 64 );
clrmp = colormap('jet');
clrmp(32,:) = [1 1 1];

kk = find( area1_pval > pcriterion );
area1_rho(kk) = zeros(size(kk));

kk = find( area2_pval > pcriterion );
area2_rho(kk) = zeros(size(kk));

kk = find( area1_area2_pval > pcriterion );
area1_area2_rho(kk) = zeros(size(kk));


area1_Ncomponent = size(area1_rho,1);
area1_cluster_assignment = area1_IdxClusters;

area1_source_order = 1:area1_Ncomponent;
[area1_target_order, area1_cluster_idx_sorted] = TargetOrder( area1_cluster_assignment, area1_Nclusters );

[area1_rho_sorted, area1_source_order] = SortMatrix( area1_rho, area1_source_order, area1_target_order, area1_Ncomponent );


area2_Ncomponent = size(area2_rho,1);
area2_cluster_assignment = area2_IdxClusters;

area2_source_order = 1:area2_Ncomponent;
[area2_target_order, area2_cluster_idx_sorted] = TargetOrder( area2_cluster_assignment, area2_Nclusters );

[area2_rho_sorted, area2_source_order] = SortMatrix( area2_rho, area2_source_order, area2_target_order, area2_Ncomponent );


X = size(area1_rho,1);
Y = size(area2_rho,1);
A12_rho = zeros(X,Y);
for i=1:X
    for j=1:Y
        A12_rho(i,j) = area1_area2_rho(area1_source_order(i),area2_source_order(j));
    end
end

X = size(area2_rho,1);
Y = size(area1_rho,1);
A21_rho = zeros(X,Y);
for i=1:X
    for j=1:Y
        A21_rho(i,j) = area1_area2_rho(area1_source_order(j),area2_source_order(i));
    end
end

f = figure;
imagesc(area1_rho_sorted,CLIM);
disp(strcat('SQ-',area1_label,':',num2str(sumsqr(area1_rho_sorted))));
colormap(clrmp);
caxis(CLIM);
title(area1_label, 'FontSize', fs );
h = colorbar;
set(h,'YLim',CLIM, 'YTick',CLIM, 'PlotBoxAspectRatio', [1 20 1]);
print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters','-',area1_label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters','-',area1_label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters','-',area1_label,'.pdf'));

f = figure;
area1_rho_sorted_positive = area1_rho_sorted;
area1_rho_sorted_positive(find(area1_rho_sorted<0)) = 0;
imagesc(area1_rho_sorted_positive,CLIM);
disp(strcat('SQ-positive-',area1_label,':',num2str(sumsqr(area1_rho_sorted_positive))));
colormap(clrmp);
caxis(CLIM);
title(strcat(area1_label,'-','positive'), 'FontSize', fs );
h = colorbar;
set(h,'YLim',CLIM, 'YTick',CLIM, 'PlotBoxAspectRatio', [1 20 1]);
print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Positive','-',area1_label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Positive','-',area1_label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Positive','-',area1_label,'.pdf'));

f = figure;
area1_rho_sorted_negative = area1_rho_sorted;
area1_rho_sorted_negative(find(area1_rho_sorted>0)) = 0;
imagesc(area1_rho_sorted_negative,CLIM);
disp(strcat('SQ-negative-',area1_label,':',num2str(sumsqr(area1_rho_sorted_negative))));
colormap(clrmp);
caxis(CLIM);
title(strcat(area1_label,'-','negative'), 'FontSize', fs );
h = colorbar;
set(h,'YLim',CLIM, 'YTick',CLIM, 'PlotBoxAspectRatio', [1 20 1]);
print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Negative','-',area1_label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Negative','-',area1_label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Negative','-',area1_label,'.pdf'));

f = figure;
imagesc(area2_rho_sorted,CLIM);
disp(strcat('SQ-',area2_label,':',num2str(sumsqr(area2_rho_sorted))));
colormap(clrmp);
caxis(CLIM);
title(area2_label, 'FontSize', fs );
h = colorbar;
set(h,'YLim',CLIM, 'YTick',CLIM, 'PlotBoxAspectRatio', [1 20 1]);
print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters','-',area2_label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters','-',area2_label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters','-',area2_label,'.pdf'));

f = figure;
area2_rho_sorted_positive = area2_rho_sorted;
area2_rho_sorted_positive(find(area2_rho_sorted<0)) = 0;
imagesc(area2_rho_sorted_positive,CLIM);
disp(strcat('SQ-positive-',area2_label,':',num2str(sumsqr(area2_rho_sorted_positive))));
colormap(clrmp);
caxis(CLIM);
title(strcat(area2_label,'-','positive'), 'FontSize', fs );
h = colorbar;
set(h,'YLim',CLIM, 'YTick',CLIM, 'PlotBoxAspectRatio', [1 20 1]);
print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Positive','-',area2_label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Positive','-',area2_label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Positive','-',area2_label,'.pdf'));

f = figure;
area2_rho_sorted_negative = area2_rho_sorted;
area2_rho_sorted_negative(find(area2_rho_sorted>0)) = 0;
imagesc(area2_rho_sorted_negative,CLIM);
disp(strcat('SQ-negative-',area2_label,':',num2str(sumsqr(area2_rho_sorted_negative))));
colormap(clrmp);
caxis(CLIM);
title(strcat(area2_label,'-','negative'), 'FontSize', fs );
h = colorbar;
set(h,'YLim',CLIM, 'YTick',CLIM, 'PlotBoxAspectRatio', [1 20 1]);
print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Negative','-',area2_label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Negative','-',area2_label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Negative','-',area2_label,'.pdf'));

f = figure;
imagesc(A12_rho,CLIM);
disp(strcat('SQ-','X-',area1_label,'-Y-',area2_label,':',num2str(sumsqr(A12_rho))));
colormap(clrmp);
caxis(CLIM);
title(strcat(area1_label,'-',area2_label), 'FontSize', fs );
h = colorbar;
set(h,'YLim',CLIM, 'YTick',CLIM, 'PlotBoxAspectRatio', [1 20 1]);
print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters','-X-',area1_label,'-Y-',area2_label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters','-X-',area1_label,'-Y-',area2_label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters','-X-',area1_label,'-Y-',area2_label,'.pdf'));

f = figure;
A12_rho_positive = A12_rho;
A12_rho_positive(find(A12_rho<0)) = 0;
imagesc(A12_rho_positive,CLIM);
disp(strcat('SQ-positive-',area1_label,'-',area2_label,':',num2str(sumsqr(A12_rho_positive))));
colormap(clrmp);
caxis(CLIM);
title(strcat(area1_label,'-',area2_label,'-','positive'), 'FontSize', fs );
h = colorbar;
set(h,'YLim',CLIM, 'YTick',CLIM, 'PlotBoxAspectRatio', [1 20 1]);
print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Positive','-X-',area1_label,'-Y-',area2_label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Positive','-X-',area1_label,'-Y-',area2_label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Positive','-X-',area1_label,'-Y-',area2_label,'.pdf'));

f = figure;
A12_rho_negative = A12_rho;
A12_rho_negative(find(A12_rho<0)) = 0;
imagesc(A12_rho_negative,CLIM);
disp(strcat('SQ-positive-',area1_label,'-',area2_label,':',num2str(sumsqr(A12_rho_negative))));
colormap(clrmp);
caxis(CLIM);
title(strcat(area1_label,'-',area2_label,'-','positive'), 'FontSize', fs );
h = colorbar;
set(h,'YLim',CLIM, 'YTick',CLIM, 'PlotBoxAspectRatio', [1 20 1]);
print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Negative','-X-',area1_label,'-Y-',area2_label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Negative','-X-',area1_label,'-Y-',area2_label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Negative','-X-',area1_label,'-Y-',area2_label,'.pdf'));

f = figure;
imagesc(A21_rho,CLIM);
disp(strcat('SQ-','X-',area2_label,'-Y-',area1_label,':',num2str(sumsqr(A21_rho))));
colormap(clrmp);
caxis(CLIM);
title(strcat(area2_label,'-',area1_label), 'FontSize', fs );
h = colorbar;
set(h,'YLim',CLIM, 'YTick',CLIM, 'PlotBoxAspectRatio', [1 20 1]);
print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters','-X-',area2_label,'-Y-',area1_label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters','-X-',area2_label,'-Y-',area1_label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters','-X-',area2_label,'-Y-',area1_label,'.pdf'));

f = figure;
A21_rho_positive = A21_rho;
A21_rho_positive(find(A21_rho<0)) = 0;
imagesc(A21_rho_positive,CLIM);
disp(strcat('SQ-positive-',area2_label,'-',area1_label,':',num2str(sumsqr(A21_rho_positive))));
colormap(clrmp);
caxis(CLIM);
title(strcat(area2_label,'-',area1_label,'-','positive'), 'FontSize', fs );
h = colorbar;
set(h,'YLim',CLIM, 'YTick',CLIM, 'PlotBoxAspectRatio', [1 20 1]);
print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Positive','-X-',area2_label,'-Y-',area1_label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Positive','-X-',area2_label,'-Y-',area1_label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Positive','-X-',area2_label,'-Y-',area1_label,'.pdf'));

f = figure;
A21_rho_negative = A12_rho;
A21_rho_negative(find(A21_rho<0)) = 0;
imagesc(A21_rho_negative,CLIM);
disp(strcat('SQ-positive-',area2_label,'-',area1_label,':',num2str(sumsqr(A21_rho_negative))));
colormap(clrmp);
caxis(CLIM);
title(strcat(area2_label,'-',area1_label,'-','positive'), 'FontSize', fs );
h = colorbar;
set(h,'YLim',CLIM, 'YTick',CLIM, 'PlotBoxAspectRatio', [1 20 1]);
print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Negative','-X-',area2_label,'-Y-',area1_label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Negative','-X-',area2_label,'-Y-',area1_label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Clusters-Negative','-X-',area2_label,'-Y-',area1_label,'.pdf'));

 

close all

end

function plotFCBtwConditionFractionDistance(settings,area1_ID,area2_ID)

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

nTR = 300;
pcriterion = 0.001;

AAL_ROI = load_roi.ROI;

label_ROI = AAL_ROI(area1_ID).Nom_L;
area1_label = strrep(label_ROI,'_','-');
disp(label_ROI);

label_ROI = AAL_ROI(area2_ID).Nom_L;
area2_label = strrep(label_ROI,'_','-');
disp(label_ROI);

area1.FC = load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area1_label,'.mat'));

area2.FC = load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area2_label,'.mat'));

area1_area2.FC = load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area1_label,'-',area2_label,'.mat'));

%% AREA 1

all_distance_area1 = getAllDistance(AAL_img,AAL_ROI,area1_ID,area1_ID);
all_distance_area1(all_distance_area1==0) = [];

rho1 = area1.FC.FC_Voxels.run(1).rho_MOT4;
rho2 = area1.FC.FC_Voxels.run(1).rho_RestingState;
n1 = nTR;
n2 = nTR;
pval = testBtwRhos(rho1,rho2,n1,n2);

label = strcat(area1_label,'-','Run-1');
significative = double( pval < pcriterion );
significative(find(eye(size(significative)))) = 0;
fraction = getFractionPerDistance(AAL_img,AAL_ROI,area1_ID,area1_ID,significative);
plotFractionPerDistance(settings,fraction,label,all_distance_area1);

rho1 = area1.FC.FC_Voxels.run(2).rho_MOT4;
rho2 = area1.FC.FC_Voxels.run(2).rho_RestingState;
n1 = nTR;
n2 = nTR;
pval = testBtwRhos(rho1,rho2,n1,n2);

label = strcat(area1_label,'-','Run-2');
significative = double( pval < pcriterion );
significative(find(eye(size(significative)))) = 0;
fraction = getFractionPerDistance(AAL_img,AAL_ROI,area1_ID,area1_ID,significative);
plotFractionPerDistance(settings,fraction,label,all_distance_area1);

%% AREA 2

all_distance_area2 = getAllDistance(AAL_img,AAL_ROI,area2_ID,area2_ID);
all_distance_area2(all_distance_area2==0) = [];

rho1 = area2.FC.FC_Voxels.run(1).rho_MOT4;
rho2 = area2.FC.FC_Voxels.run(1).rho_RestingState;
n1 = nTR;
n2 = nTR;
pval = testBtwRhos(rho1,rho2,n1,n2);

label = strcat(area2_label,'-','Run-1');
significative = double( pval < pcriterion );
significative(find(eye(size(significative)))) = 0;
fraction = getFractionPerDistance(AAL_img,AAL_ROI,area2_ID,area2_ID,significative);
plotFractionPerDistance(settings,fraction,label,all_distance_area2);

rho1 = area2.FC.FC_Voxels.run(2).rho_MOT4;
rho2 = area2.FC.FC_Voxels.run(2).rho_RestingState;
n1 = nTR;
n2 = nTR;
pval = testBtwRhos(rho1,rho2,n1,n2);

label = strcat(area2_label,'-','Run-2');
significative = double( pval < pcriterion );
significative(find(eye(size(significative)))) = 0;
fraction = getFractionPerDistance(AAL_img,AAL_ROI,area2_ID,area2_ID,significative);
plotFractionPerDistance(settings,fraction,label,all_distance_area2);

%% AREA 1 & 2

all_distance_area1_area2 = getAllDistance(AAL_img,AAL_ROI,area1_ID,area2_ID);
all_distance_area1_area2(all_distance_area1_area2==0) = [];

rho1 = area1_area2.FC.FC_Voxels.run(1).rho_MOT4;
rho2 = area1_area2.FC.FC_Voxels.run(1).rho_RestingState;
n1 = nTR;
n2 = nTR;
pval = testBtwRhos(rho1,rho2,n1,n2);

label = strcat(area1_label,'-',area2_label,'-','Run-1');
significative = double( pval < pcriterion );
significative(find(eye(size(significative)))) = 0;
fraction = getFractionPerDistance(AAL_img,AAL_ROI,area1_ID,area2_ID,significative);
plotFractionPerDistance(settings,fraction,label,all_distance_area1_area2);

rho1 = area1_area2.FC.FC_Voxels.run(2).rho_MOT4;
rho2 = area1_area2.FC.FC_Voxels.run(2).rho_RestingState;
n1 = nTR;
n2 = nTR;
pval = testBtwRhos(rho1,rho2,n1,n2);

label = strcat(area1_label,'-',area2_label,'-','Run-2');
significative = double( pval < pcriterion );
significative(find(eye(size(significative)))) = 0;
fraction = getFractionPerDistance(AAL_img,AAL_ROI,area1_ID,area2_ID,significative);
plotFractionPerDistance(settings,fraction,label,all_distance_area1_area2);


end

function plotFCDistributionAndDifference(settings,area1_ID,area2_ID)

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

label_ROI = AAL_ROI(area1_ID).Nom_L;
area1_label = strrep(label_ROI,'_','-');
disp(label_ROI);

label_ROI = AAL_ROI(area2_ID).Nom_L;
area2_label = strrep(label_ROI,'_','-');
disp(label_ROI);

area1.FC = load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area1_label,'.mat'));

area2.FC = load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area2_label,'.mat'));

area1_area2.FC = load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area1_label,'-',area2_label,'.mat'));

label = strcat(area1_label,'-','Run-1');
rho_MOT4 = area1.FC.FC_Voxels.run(1).rho_MOT4;
rho_RestingState = area1.FC.FC_Voxels.run(1).rho_RestingState;

plotFCDistributionAndDifferenceOnePlot(settings,rho_MOT4,rho_RestingState,label);

label = strcat(area1_label,'-','Run-2');
rho_MOT4 = area1.FC.FC_Voxels.run(2).rho_MOT4;
rho_RestingState = area1.FC.FC_Voxels.run(2).rho_RestingState;

plotFCDistributionAndDifferenceOnePlot(settings,rho_MOT4,rho_RestingState,label);

label = strcat(area2_label,'-','Run-1');
rho_MOT4 = area2.FC.FC_Voxels.run(1).rho_MOT4;
rho_RestingState = area2.FC.FC_Voxels.run(1).rho_RestingState;

plotFCDistributionAndDifferenceOnePlot(settings,rho_MOT4,rho_RestingState,label);

label = strcat(area2_label,'-','Run-2');
rho_MOT4 = area2.FC.FC_Voxels.run(2).rho_MOT4;
rho_RestingState = area2.FC.FC_Voxels.run(2).rho_RestingState;

plotFCDistributionAndDifferenceOnePlot(settings,rho_MOT4,rho_RestingState,label);

label = strcat(area1_label,'-',area2_label,'-','Run-1');
rho_MOT4 = area1_area2.FC.FC_Voxels.run(1).rho_MOT4;
rho_RestingState = area1_area2.FC.FC_Voxels.run(1).rho_RestingState;

plotFCDistributionAndDifferenceOnePlot(settings,rho_MOT4,rho_RestingState,label);

label = strcat(area1_label,'-',area2_label,'-','Run-2');
rho_MOT4 = area1_area2.FC.FC_Voxels.run(2).rho_MOT4;
rho_RestingState = area1_area2.FC.FC_Voxels.run(2).rho_RestingState;

plotFCDistributionAndDifferenceOnePlot(settings,rho_MOT4,rho_RestingState,label);

end

function plotFCDistributionAndDifferenceOnePlot(settings,rho_MOT4,rho_RestingState,label)

nTR = 300;
pcriterion = 0.001;
points = 1000;

max_y = 0.02;

max_val = 0.5 * log( (1 + 0.99)./(1 - 0.99) );
min_val = 0.5 * log( (1 + (-0.99))./(1 - (-0.99)) );
xi = linspace(min_val,max_val,points);

rho_MOT4(find(eye(size(rho_MOT4)))) = NaN;
rho_RestingState(find(eye(size(rho_RestingState)))) = NaN;
rho_MOT4_fisher = 0.5 * log( (1 + rho_MOT4)./(1 - rho_MOT4) );
rho_RestingState_fisher = 0.5 * log( (1 + rho_RestingState)./(1 - rho_RestingState) );

max_val_MOT4 = max(rho_MOT4_fisher(:));
min_val_MOT4 = min(rho_MOT4_fisher(:));
xi_MOT4 = linspace(min_val_MOT4,max_val_MOT4,points);

max_val_RestingState = max(rho_RestingState_fisher(:));
min_val_RestingState = min(rho_RestingState_fisher(:));
xi_RestingState = linspace(min_val_RestingState,max_val_RestingState,points);

max_val_MOT4_RestingState = max(rho_MOT4_fisher(:) - rho_RestingState_fisher(:));
min_val_MOT4_RestingState = min(rho_MOT4_fisher(:) - rho_RestingState_fisher(:));
xi_MOT4_RestingState = linspace(min_val_MOT4_RestingState,max_val_MOT4_RestingState,points);

% min_val = min([min_val_MOT4, min_val_RestingState, min_val_MOT4_RestingState]);
% max_val = max([max_val_MOT4, max_val_RestingState, max_val_MOT4_RestingState]);

% k_MOT4 = ksdensity(rho_MOT4_fisher(:),xi_MOT4);
% k_RestingState = ksdensity(rho_RestingState_fisher(:),xi_RestingState);
% k = ksdensity(rho_MOT4_fisher(:)-rho_RestingState_fisher(:),xi_MOT4_RestingState);

interval = xi(2) - xi(1);

k_MOT4 = ksdensity(rho_MOT4_fisher(:),xi);
k_RestingState = ksdensity(rho_RestingState_fisher(:),xi);
k = ksdensity(rho_MOT4_fisher(:)-rho_RestingState_fisher(:),xi);

k_MOT4 = k_MOT4.*interval;
k_RestingState = k_RestingState.*interval;
k = k.*interval;

P1 = k_MOT4 + eps;
P2 = k_RestingState + eps;
KL = sum(P1.*log2(P1./P2));

% k_MOT4 = k_MOT4./max(k_MOT4);
% k_RestingState = k_RestingState./max(k_RestingState);
% k = k./max(k);

% Y1 = k_MOT4;
% Y2 = k_RestingState;
% cost_name = 'KL_kNN_k';
% mult = 1;
% co = D_initialization(cost_name,mult);
% D = D_estimation(Y1,Y2,co);

std_k = std(k);

idx = find(k_MOT4 == 0);
k_MOT4(idx) = [];
xi_MOT4 = xi;
xi_MOT4(idx) = [];

idx = find(k_RestingState == 0);
k_RestingState(idx) = [];
xi_RestingState = xi;
xi_RestingState(idx) = [];

idx = find(k == 0);
k(idx) = [];
xi_MOT4_RestingState = xi;
xi_MOT4_RestingState(idx) = [];

z_fisher = ( xi_MOT4_RestingState ) ./ ( (1/(nTR-3)) + (1/(nTR-3)) );

z = abs(norminv(pcriterion,0,1));

idx_positive = find(z_fisher>z);
idx_negative = find(z_fisher<-z);

percent_significant = 1;
n_idx_positive_significant = ceil((length(idx_positive)/100)*(percent_significant/2));
n_idx_negative_significant = ceil((length(idx_negative)/100)*(percent_significant/2));

xx = xi(end-300);
yy1 = (max_y/2);
yy2 = (max_y/2) + 0.1*max_y;

f = figure;
plot(xi_MOT4,k_MOT4,'b:','LineWidth',2);
hold on
plot(xi_RestingState,k_RestingState,'g:','LineWidth',2);
plot(xi_MOT4_RestingState,k,'r:','LineWidth',2);
plot(xi_MOT4_RestingState(idx_positive(end-n_idx_positive_significant+1:end)),k(idx_positive(end-n_idx_positive_significant+1:end)),'ro','MarkerFaceColor','r'); 
plot(xi_MOT4_RestingState(idx_negative(1:n_idx_negative_significant)),k(idx_negative(1:n_idx_negative_significant)),'ro','MarkerFaceColor','r'); 
xlim([min_val max_val]);
ylim([-0.1*max_y (max_y+0.1*max_y)]);
xlabel('Pearson - Fisher Z Transformed');
ylabel('Probability');
legend({'High Attention' 'Resting State' 'Attention-Rest'});
text(xx,yy1,strcat('KL(Att||Rest):',num2str(KL)));
text(xx,yy2,strcat('std:',num2str(std_k)));
title(label);
N = 4;
set(gca,'YTIck',linspace(0,max_y,N));

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-Distr.-High-Rest','-',label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-Distr.-High-Rest','-',label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-Distr.-High-Rest','-',label,'.pdf'));

end

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
    
    significant = Adecrease ~= 0;
    AdecreaseFraction = sum(significant,1)./Ncomponent;

    significant = Aincrease ~= 0;
    AincreaseFraction = sum(significant,1)./Ncomponent;
    
    significant = Rdecrease ~= 0;
    RdecreaseFraction = sum(significant,1)./Ncomponent;

    significant = Rincrease ~= 0;
    RincreaseFraction = sum(significant,1)./Ncomponent;
    
end

function plotFCDifferencePositiveNegativeFraction(settings,area_ID)

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

label_ROI = AAL_ROI(area_ID).Nom_L;
area_label = strrep(label_ROI,'_','-');
disp(label_ROI);

points = 1000;

fs = 12;

shouldIDoCluster = 1;

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area_label,'.mat'));

for irun=1:2
    
    rho_att = FC_Voxels.run(irun).rho_MOT4;
    pval_att = FC_Voxels.run(irun).pval_MOT4;
    
    rho_rest = FC_Voxels.run(irun).rho_RestingState;
    pval_rest = FC_Voxels.run(irun).pval_RestingState; 
    
    [AincreaseFraction, RincreaseFraction, AdecreaseFraction, RdecreaseFraction, Achanges, Rchanges, Aincrease, Rincrease, Adecrease, Rdecrease, Dfisher, Afisher, Rfisher, Asorted, Rsorted, target_order, A_source_order, R_source_order, source_order, Ncomponent, cluster_assignment, Ncluster] = getFCDifferencePositiveNegativeFraction(rho_att,pval_att,rho_rest,pval_rest,shouldIDoCluster);

    significant = Rchanges ~= 0;
    R_distance = getFractionPerDistanceSorted(AAL_img,AAL_ROI,area_ID,area_ID,significant,R_source_order);
    
    significant = Achanges ~= 0;
    A_distance = getFractionPerDistanceSorted(AAL_img,AAL_ROI,area_ID,area_ID,significant,A_source_order);

    all_distance = getAllDistance(AAL_img,AAL_ROI,area_ID,area_ID);
    
    rlimit = 0.5;

    rmin = -rlimit + rlimit/32;
    rmax =  rlimit;

    Tick = [1 Ncomponent];
    TickLabel = { '1' ; num2str(Ncomponent) };
    
    Ri = zeros( Ncomponent+1, Ncomponent+1 );
    Ri(1:Ncomponent, 1:Ncomponent) = Rincrease;

    Ai = zeros( Ncomponent+1, Ncomponent+1 );
    Ai(1:Ncomponent, 1:Ncomponent) = Aincrease;
    
    Rd = zeros( Ncomponent+1, Ncomponent+1 );
    Rd(1:Ncomponent, 1:Ncomponent) = Rdecrease;

    Ad = zeros( Ncomponent+1, Ncomponent+1 );
    Ad(1:Ncomponent, 1:Ncomponent) = Adecrease;
    
    Rc = zeros( Ncomponent+1, Ncomponent+1 );
    Rc(1:Ncomponent, 1:Ncomponent) = Rchanges;

    Ac = zeros( Ncomponent+1, Ncomponent+1 );
    Ac(1:Ncomponent, 1:Ncomponent) = Achanges;
    
    xi = linspace(min(all_distance),max(all_distance),points);
    interval = xi(2) - xi(1);
    
    [A_k,xi] = ksdensity(A_distance,xi);
    [R_k,xi] = ksdensity(R_distance,xi);
    [k,xi] = ksdensity(all_distance,xi);
    
    A_k = A_k.*interval;
    R_k = R_k.*interval;
    k = k.*interval;
    
    f = figure;
    clrmp = colormap('jet');
    clrmp(32,:) = [1 1 1];
    hold on;
    caxis([rmin rmax]);
    h = pcolor(0.5:Ncomponent+0.5, 0.5:Ncomponent+0.5, Ai );
    set( h, 'EdgeColor', 'none');
    colormap(clrmp);
    hold off;
    axis 'square'; 
    axis([0 Ncomponent+1 0 Ncomponent+1]);
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel, 'FontSize', fs );
    title(strcat(area_label,'-Attention-Increase','-Run-',int2str(irun)), 'FontSize', fs );
    h = colorbar;
    set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1], 'FontSize', fs);
    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Ai-',area_label,'-R',int2str(irun),'.jpg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Ai-',area_label,'-R',int2str(irun),'.eps'));
    print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Ai-',area_label,'-R',int2str(irun),'.pdf'));

    f = figure;
    clrmp = colormap('jet');
    clrmp(32,:) = [1 1 1];
    hold on;
    caxis([rmin rmax]);
    h = pcolor(0.5:Ncomponent+0.5, 0.5:Ncomponent+0.5, Ad );
    set( h, 'EdgeColor', 'none');
    colormap(clrmp);
    hold off;
    axis 'square'; 
    axis([0 Ncomponent+1 0 Ncomponent+1]);
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel, 'FontSize', fs );
    title(strcat(area_label,'-Attention-Decrease','-Run-',int2str(irun)), 'FontSize', fs );
    h = colorbar;
    set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1], 'FontSize', fs);
    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Ad-',area_label,'-R',int2str(irun),'.jpg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Ad-',area_label,'-R',int2str(irun),'.eps'));
    print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Ad-',area_label,'-R',int2str(irun),'.pdf'));

    f = figure;
    clrmp = colormap('jet');
    clrmp(32,:) = [1 1 1];
    hold on;
    caxis([rmin rmax]);
    h = pcolor(0.5:Ncomponent+0.5, 0.5:Ncomponent+0.5, Ri );
    set( h, 'EdgeColor', 'none');
    colormap(clrmp);
    hold off;
    axis 'square'; 
    axis([0 Ncomponent+1 0 Ncomponent+1]);
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel, 'FontSize', fs );
    title(strcat(area_label,'-RestingState-Increase','-Run-',int2str(irun)), 'FontSize', fs );
    h = colorbar;
    set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1], 'FontSize', fs);
    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Ri-',area_label,'-R',int2str(irun),'.jpg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Ri-',area_label,'-R',int2str(irun),'.eps'));
    print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Ri-',area_label,'-R',int2str(irun),'.pdf'));

    f = figure;
    clrmp = colormap('jet');
    clrmp(32,:) = [1 1 1];
    hold on;
    caxis([rmin rmax]);
    h = pcolor(0.5:Ncomponent+0.5, 0.5:Ncomponent+0.5, Rd );
    set( h, 'EdgeColor', 'none');
    colormap(clrmp);
    hold off;
    axis 'square'; 
    axis([0 Ncomponent+1 0 Ncomponent+1]);
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel, 'FontSize', fs );
    title(strcat(area_label,'-RestingState-Decrease','-Run-',int2str(irun)), 'FontSize', fs );
    h = colorbar;
    set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1], 'FontSize', fs);
    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Rd-',area_label,'-R',int2str(irun),'.jpg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Rd-',area_label,'-R',int2str(irun),'.eps'));
    print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Rd-',area_label,'-R',int2str(irun),'.pdf'));

    f = figure;
    clrmp = colormap('jet');
    clrmp(32,:) = [1 1 1];
    hold on;
    caxis([rmin rmax]);
    h = pcolor(0.5:Ncomponent+0.5, 0.5:Ncomponent+0.5, Rc );
    set( h, 'EdgeColor', 'none');
    colormap(clrmp);
    hold off;
    axis 'square'; 
    axis([0 Ncomponent+1 0 Ncomponent+1]);
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel, 'FontSize', fs );
    title(strcat(area_label,'-RestingState-Changes','-Run-',int2str(irun)), 'FontSize', fs );
    h = colorbar;
    set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1], 'FontSize', fs);
    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Rc-',area_label,'-R',int2str(irun),'.jpg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Rc-',area_label,'-R',int2str(irun),'.eps'));
    print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Rc-',area_label,'-R',int2str(irun),'.pdf'));

    f = figure;
    clrmp = colormap('jet');
    clrmp(32,:) = [1 1 1];
    hold on;
    caxis([rmin rmax]);
    h = pcolor(0.5:Ncomponent+0.5, 0.5:Ncomponent+0.5, Ac );
    set( h, 'EdgeColor', 'none');
    colormap(clrmp);
    hold off;
    axis 'square'; 
    axis([0 Ncomponent+1 0 Ncomponent+1]);
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel, 'FontSize', fs );
    title(strcat(area_label,'-Attention-Changes','-Run-',int2str(irun)), 'FontSize', fs );
    h = colorbar;
    set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1], 'FontSize', fs);
    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Ac-',area_label,'-R',int2str(irun),'.jpg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Ac-',area_label,'-R',int2str(irun),'.eps'));
    print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Diff-','Ac-',area_label,'-R',int2str(irun),'.pdf'));

    f = figure;
    plot(AincreaseFraction);
    title(strcat(area_label,'-Attention-Fraction-Increase','-Run-',int2str(irun)));
    xlabel('Voxels');
    ylabel('Fraction');
    xlim([0 Ncomponent]);
    ylim([0 1]);
    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Ai-',area_label,'-R',int2str(irun),'.jpg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Ai-',area_label,'-R',int2str(irun),'.eps'));
    print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Ai-',area_label,'-R',int2str(irun),'.pdf'));

    f = figure;
    plot(AdecreaseFraction);
    title(strcat(area_label,'-Attention-Fraction-Decrease','-Run-',int2str(irun)));
    xlabel('Voxels');
    ylabel('Fraction');
    xlim([0 Ncomponent]);
    ylim([0 1]);
    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Ad-',area_label,'-R',int2str(irun),'.jpg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Ad-',area_label,'-R',int2str(irun),'.eps'));
    print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Ad-',area_label,'-R',int2str(irun),'.pdf'));

    f = figure;
    plot(RincreaseFraction);
    title(strcat(area_label,'-Resting-Fraction-Increase','-Run-',int2str(irun)));
    xlabel('Voxels');
    ylabel('Fraction');
    xlim([0 Ncomponent]);
    ylim([0 1]);
    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Ri-',area_label,'-R',int2str(irun),'.jpg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Ri-',area_label,'-R',int2str(irun),'.eps'));
    print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Ri-',area_label,'-R',int2str(irun),'.pdf'));

    f = figure;
    plot(RdecreaseFraction);
    title(strcat(area_label,'-Resting-Fraction-Decrease','-Run-',int2str(irun)));
    xlabel('Voxels');
    ylabel('Fraction');
    xlim([0 Ncomponent]);
    ylim([0 1]);
    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Rd-',area_label,'-R',int2str(irun),'.jpg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Rd-',area_label,'-R',int2str(irun),'.eps'));
    print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Rd-',area_label,'-R',int2str(irun),'.pdf'));

    f = figure; 
    plot(xi,A_k,'b');
    hold on
    plot(xi,k,'r');
    title(strcat(area_label,'-Attention-Distance','-Run-',int2str(irun)));
    xlabel('Distance (mm)');
    ylabel('Fraction');
    xlim([0 max(xi)]);
    ylim([0 max([max(A_k) max(k)])]);
    legend({'Changed Pairs' 'All Pairs'});
    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Adis-',area_label,'-R',int2str(irun),'.jpg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Adis-',area_label,'-R',int2str(irun),'.eps'));
    print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Adis-',area_label,'-R',int2str(irun),'.pdf'));

    f = figure;
    plot(xi,R_k,'b');
    hold on
    plot(xi,k,'r');
    title(strcat(area_label,'-Resting-Distance','-Run-',int2str(irun)));
    xlabel('Distance (mm)');
    ylabel('Fraction');
    xlim([0 max(xi)]);
    ylim([0 max([max(R_k) max(k)])]);
    legend({'Changed Pairs' 'All Pairs'});
    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Rdis-',area_label,'-R',int2str(irun),'.jpg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Rdis-',area_label,'-R',int2str(irun),'.eps'));
    print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fra-','Rdis-',area_label,'-R',int2str(irun),'.pdf'));

    close all
    
end

end
%% funtions for sorting clusters of rows/columns by size

function computeAllFractionsFCDifferencePositiveNegative(settings,idx_areas)

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
%load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

shouldIDoCluster = 0;

for iROI=1:length(idx_areas)
    
    tic 
    
    label_ROI = AAL_ROI(idx_areas(iROI)).Nom_L;
    area_label = strrep(label_ROI,'_','-');
    disp(label_ROI);

    load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-',area_label,'.mat'));

    for irun=1:2

        rho_att = FC_Voxels.run(irun).rho_MOT4;
        pval_att = FC_Voxels.run(irun).pval_MOT4;

        rho_rest = FC_Voxels.run(irun).rho_RestingState;
        pval_rest = FC_Voxels.run(irun).pval_RestingState; 

        [AincreaseFraction, RincreaseFraction, AdecreaseFraction, RdecreaseFraction, Achanges, Rchanges, Aincrease, Rincrease, Adecrease, Rdecrease, Dfisher, Afisher, Rfisher, Asorted, Rsorted, target_order, A_source_order, R_source_order, source_order, Ncomponent, cluster_assignment, Ncluster] = getFCDifferencePositiveNegativeFraction(rho_att,pval_att,rho_rest,pval_rest,shouldIDoCluster);

        FC_Fraction.run(irun).AincreaseFraction = AincreaseFraction;
        FC_Fraction.run(irun).RincreaseFraction = RincreaseFraction;
        
        FC_Fraction.run(irun).AdecreaseFraction = AdecreaseFraction;
        FC_Fraction.run(irun).RdecreaseFraction = RdecreaseFraction;
        
    end
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fraction-',area_label,'.mat'),'source_order','area_label','target_order','A_source_order','R_source_order','FC_Fraction');
    
    clear FC_Voxels
    clear FC_Fraction
    
    toc
    
end

end

function generate3DimgFractions(settings)

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26 27 28];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [17 18 19 20 57 58 59 60 61 62 63 64 65 66 67 68 69 70];
idx_areas = [idx_frontal,idx_occipital,idx_parietal];
%idx_areas = 1:90;

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
%load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

%nROI = 90;

for irun=1:2
    
    AincreaseFraction3D = zeros(size(AAL_img));
    AdecreaseFraction3D = zeros(size(AAL_img));

    RincreaseFraction3D = zeros(size(AAL_img));
    RdecreaseFraction3D = zeros(size(AAL_img));
    
    for iROI=1:length(idx_areas)
    
        label_ROI = AAL_ROI(idx_areas(iROI)).Nom_L;
        area_label = strrep(label_ROI,'_','-');
        disp(label_ROI);

        idx_AAL = AAL_ROI(idx_areas(iROI)).ID;
        idx_voxels = find(AAL_img == idx_AAL);

        load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fraction-',area_label,'.mat'));

        for iVoxel=1:length(idx_voxels)

            [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

            AincreaseFraction3D(idxx,idxy,idxz) = FC_Fraction.run(irun).AincreaseFraction(iVoxel);
            RincreaseFraction3D(idxx,idxy,idxz) = FC_Fraction.run(irun).RincreaseFraction(iVoxel);
            
            AdecreaseFraction3D(idxx,idxy,idxz) = FC_Fraction.run(irun).AdecreaseFraction(iVoxel);
            RdecreaseFraction3D(idxx,idxy,idxz) = FC_Fraction.run(irun).RdecreaseFraction(iVoxel);

        end
    
    end
    
   
    scl_slope = 1;
    scl_inter = 0; 
    dim = [size(AAL_img,1),size(AAL_img,2),size(AAL_img,3)];
    dtype = 'FLOAT32';
    offset = 0;
    nifti_file.mat = load_aal.mat;
    nifti_file.mat_intent = load_aal.mat_intent;
    nifti_file.mat0 = load_aal.mat0;
    nifti_file.mat0_intent = load_aal.mat0_intent;

    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fraction-Ai-R',int2str(irun),'.nii');
    descrip = 'AincreaseFraction3D';
    input_data = AincreaseFraction3D; 
    lowhigh_save_image;
    
    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fraction-Ad-R',int2str(irun),'.nii');
    descrip = 'AdecreaseFraction3D';
    input_data = AdecreaseFraction3D; 
    lowhigh_save_image;
    
    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fraction-Ri-R',int2str(irun),'.nii');
    descrip = 'RincreaseFraction3D';
    input_data = RincreaseFraction3D; 
    lowhigh_save_image;
    
    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-per-ROI-AAL-Fraction-Rd-R',int2str(irun),'.nii');
    descrip = 'RdecreaseFraction3D';
    input_data = RdecreaseFraction3D; 
    lowhigh_save_image;
    
end

end

function [target_order, cluster_idx_sorted] = TargetOrder( cluster_assignment, ncluster )

% cluster_assignment assigns nc rows/columns to ncluster clusters

% target_order assigns nc matrix locations to nc rows/columns

% the first locations are filled with rows/columns from the largest cluster

% the next locations are filled with rows/columns from the next smaller
% cluster

% and so on

% cluster_idx_sorted lists the clusters in descending order of size ...

nc = length( cluster_assignment );

kk = [];
kk = find(isnan(cluster_assignment));

if ~isempty(kk)
    
    nNaN = length(kk);
    
    ncluster = ncluster + 1;
    cluster_assignment(kk) = ncluster;
    
end


for kc = 1 : ncluster   % loop over clusters
    
    cluster_size(kc) = length( find( cluster_assignment == kc ) );
    
end

[cluster_size_sorted, cluster_idx_sorted] = sort( cluster_size, 'descend' );   % c_size_sorted = c_size( c_idx_sorted )

cluster_size_sorted;

target_order  = [];  

for kc = 1 : ncluster % loop over sorted clusters
    'target';
    size(target_order);
    'find';
    size(find( cluster_idx_sorted(kc) == cluster_assignment ));
    target_order = [target_order; find( cluster_idx_sorted(kc) == cluster_assignment )];  % append indices of matrix rows/columns in current cluster
       
end

target_order = reshape( target_order, 1, nc );

end

%% funtions for sorting correlation matrices
function [m_sorted, source_order] = SortMatrix( m_original, source_order, target_order, nc )

% typically the rows/columns of a matrix are numbered by location, ie., 1st, 2nd, 3rd, ...
% however, when we want to rearrange a matrix, we need to distinguish between
% row/column number and location in the matrix

% source_order           lists row/column numbers in the source order (1st, 2nd, 3rd, ... location)
% target_order            lists row/column numbers in the target order (1st, 2nd, 3rd, ... location)

m_sorted = m_original;

ix_must_swap = find( source_order ~= target_order );    % find positions where component numbers differ in source and target

ix = ix_must_swap(1);     % first location to be rectified

Nswap = length( ix_must_swap );

for ic = 1 : Nswap                  % loop number of necessary swaps
    
    i_from   = ix;                  % location of component to be moved away
    c_source = source_order(ix);    % number of component to be moved away
 
    c_target = target_order(ix);                 % number of component to replace it
    i_to     = find( c_target == source_order );   % location of other component to replace it


    
    if i_from ~= i_to               % locations are different are different
        
        [m_sorted, source_order] = MSwap( m_sorted, source_order, i_from, i_to ); % swap components
                                            
                                            % now location i_from has component c_to, which is correct
                                            % but location i_to has component c_from, which is incorrect
                                            
        ix = i_to;                          % next location to be rectified                                
        
              
    elseif ic < Nswap
        ix_must_swap = find( source_order ~= target_order ); 
        ix = ix_must_swap(1);
    end
        
end

end

%% remove nan rows/columns 
function [rho_output,pval_output] = removeNaNs(rho_input,pval_input)

[Ncomponent, ~] = size( rho_input );

source_order = 1 : Ncomponent;

kk = find( isnan( rho_input(1,:) ) );

['remove ' num2str(length(kk)) ' NaN rows/columns']

for i = length(kk) : -1 : 1
    
    idelete = kk(i);
    
    [temp_rr, ~]  = MSwap( rho_input, source_order, idelete, Ncomponent );
    [temp_pr, ~]  = MSwap( pval_input, source_order, idelete, Ncomponent );
    
    clear rho_input pval_input;
    
    Ncomponent = Ncomponent - 1;
    source_order = 1 : Ncomponent;
    
    rho_input  = temp_rr( 1:Ncomponent, 1: Ncomponent);
    pval_input = temp_pr( 1:Ncomponent, 1: Ncomponent);
    
end

rho_ouput = rho_input;
pval_output = pval_input;

end

% swap locations i and k, plus update order accordingly
function [m_out, order_out]  = MSwap( m_in, order_in, i, k )

temp = m_in;

temp(i,:) = m_in(k,:);

temp(k,:) = m_in(i,:);

m_out = temp;

m_out(:,i) = temp(:,k);

m_out(:,k) = temp(:,i);

order_out = order_in;

order_out(i) = order_in(k);

order_out(k) = order_in(i);

end

function pval = testBtwRhos(rho1,rho2,n1,n2)

rho1_fisher = 0.5 * log( (1 + rho1)./(1 - rho1) );

rho2_fisher = 0.5 * log( (1 + rho2)./(1 - rho2) );

z_test = ( rho1_fisher - rho2_fisher ) ./ sqrt( (1/(n1 - 3)) + (1/(n2 - 3)) );

pval = 1-normcdf(z_test,0,1);

end

function distance = getFractionPerDistance(AAL_img,AAL_ROI,area1_ID,area2_ID,significative)

voxel_size = 2; %% mm

area1_AAL_ID = AAL_ROI(area1_ID).ID;
area2_AAL_ID = AAL_ROI(area2_ID).ID;

idx_voxels_area1 = find(AAL_img==area1_AAL_ID);
idx_voxels_area2 = find(AAL_img==area2_AAL_ID);

idx_significative = find(significative);

for ipair=1:length(idx_significative)
    
    [idx,idy] = ind2sub(size(significative),idx_significative(ipair));
    
    [idxx1,idxy1,idxz1] = ind2sub(size(AAL_img),idx_voxels_area1(idx));
    
    [idxx2,idxy2,idxz2] = ind2sub(size(AAL_img),idx_voxels_area2(idy));
    
    distance(ipair) = norm([idxx1,idxy1,idxz1] - [idxx2,idxy2,idxz2]) * voxel_size;
    
end

end

function distance = getFractionPerDistanceSorted(AAL_img,AAL_ROI,area1_ID,area2_ID,significative,source_order)

voxel_size = 2; %% mm

area1_AAL_ID = AAL_ROI(area1_ID).ID;
area2_AAL_ID = AAL_ROI(area2_ID).ID;

idx_voxels_area1 = find(AAL_img==area1_AAL_ID);
idx_voxels_area2 = find(AAL_img==area2_AAL_ID);

idx_significative = find(significative);

for ipair=1:length(idx_significative)
    
    [idx,idy] = ind2sub(size(significative),idx_significative(ipair));
    
    [idxx1,idxy1,idxz1] = ind2sub(size(AAL_img),idx_voxels_area1(source_order(idx)));
    
    [idxx2,idxy2,idxz2] = ind2sub(size(AAL_img),idx_voxels_area2(source_order(idy)));
    
    distance(ipair) = norm([idxx1,idxy1,idxz1] - [idxx2,idxy2,idxz2]) * voxel_size;
    
end

end

function all_distance = getAllDistance(AAL_img,AAL_ROI,area1_ID,area2_ID)

voxel_size = 2; %% mm

area1_AAL_ID = AAL_ROI(area1_ID).ID;
area2_AAL_ID = AAL_ROI(area2_ID).ID;

idx_voxels_area1 = find(AAL_img==area1_AAL_ID);
idx_voxels_area2 = find(AAL_img==area2_AAL_ID);

ipair = 0;
for ivoxel=1:length(idx_voxels_area1)
    
    [idxx1,idxy1,idxz1] = ind2sub(size(AAL_img),idx_voxels_area1(ivoxel));
    
    for iivoxel=1:length(idx_voxels_area2)
        
        ipair = ipair + 1;
  
        [idxx2,idxy2,idxz2] = ind2sub(size(AAL_img),idx_voxels_area2(iivoxel));
    
        all_distance(ipair) = norm([idxx1,idxy1,idxz1] - [idxx2,idxy2,idxz2]) * voxel_size;
        
    end
    
end

end

function plotFractionPerDistance(settings,fraction,label,all_distance)

f = figure;

[k,xi] = ksdensity(fraction);
plot(xi,k,'b');
hold on
[k,xi] = ksdensity(all_distance);
plot(xi,k,'r');

title(label);
xlabel('Distance (mm)');
ylabel('Probability Density');
legend({'Changed Pairs' 'All Pairs'});

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Distance-Density-',label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Distance-Density-',label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Voxels-per-ROI-AAL-Distance-Density-',label,'.pdf'));


end

% function [Tidx, Ncluster] = ClusterWithKmeans( Rvalue, Ncomponent )
% 
% num_Clust = 10;
% 
% [Tidx,C,sumd,D] = kmeans(Rvalue, num_Clust , 'distance', 'correlation', 'display', 'final','replicate',20);
% 
% Ncluster = max(Tidx);
% 
% end
