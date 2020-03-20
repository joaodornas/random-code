function real_methods_paper_v2

% ROI_spatial = getSPATIALClusters;
% ROI_MD758 = getMD758Clusters;
% 
% file_step_wavelet = 'filtered_func_data_mcf_unwarp2standard-clean-voxel';
% RestingState = getMagdeburgDataRS(file_step_wavelet);
% getMeanClusters(RestingState,ROI_spatial,'SPATIAL-WAVELET');
% getMeanClusters(RestingState,ROI_MD758,'MD758-WAVELET');
% getCorrClusters('SPATIAL-WAVELET');
% getCorrClusters('MD758-WAVELET');
% plotMeanFCClusters('SPATIAL-WAVELET');
% plotMeanFCClusters('MD758-WAVELET');
% plotVarianceClusters('SPATIAL-WAVELET');
% plotVarianceClusters('MD758-WAVELET');

% file_step_residual = 'filtered_func_data_mcf_unwarp2standard-clean-voxel-res';
% RestingState = getMagdeburgDataRS(file_step_residual);
% getMeanClusters(RestingState,ROI_spatial,'SPATIAL-RESIDUAL');
% getMeanClusters(RestingState,ROI_MD758,'MD758-RESIDUAL');
% getCorrClusters('SPATIAL-RESIDUAL');
% getCorrClusters('MD758-RESIDUAL');
% plotMeanFCClusters('SPATIAL-RESIDUAL');
% plotMeanFCClusters('MD758-RESIDUAL');
% plotVarianceClusters('SPATIAL-RESIDUAL');
% plotVarianceClusters('MD758-RESIDUAL');
% 
% file_step_filtered = 'filtered_func_data_mcf_unwarp2standard';
% RestingState = getMagdeburgDataRS(file_step_filtered);
% getMeanClusters(RestingState,ROI_spatial,'SPATIAL-FILTERED');
% getMeanClusters(RestingState,ROI_MD758,'MD758-FILTERED');
% getCorrClusters('SPATIAL-FILTERED');
% getCorrClusters('MD758-FILTERED');
% plotMeanFCClusters('SPATIAL-FILTERED');
% plotMeanFCClusters('MD758-FILTERED');
% plotVarianceClusters('SPATIAL-FILTERED');
% plotVarianceClusters('MD758-FILTERED');

% ROI_spatial = getSPATIALClusters;
% ROI_MD758 = getMD758Clusters;

% label = 'SP758';
% rho = getSPATIALClustersCorr;
% ROI = getSPATIALClusters;
% getCrossCorrClusterVoxelLevelInsideAAL(ROI,rho,label);
% 
% label = 'MD758';
% rho = getMD758ClustersCorr;
% ROI = getMD758Clusters;
% getCrossCorrClusterVoxelLevelInsideAAL(ROI,rho,label);

% plotDistributionCrossCorrClusterVoxelLevelInsideAAL; %%% It does not make sense because there is no correspondence btw. clusters of different parcellation

% plotDistributionCrossCorrrWholeBrain;

% parcellation_label = 'SP758';
% ROI = getSPATIALClusters;
% getCorrVoxelsInsideAClusterV2(ROI,parcellation_label);
% 
% parcellation_label = 'MD758';
% ROI = getMD758Clusters;
% getCorrVoxelsInsideAClusterV2(ROI,parcellation_label);

% label = 'MD758';
% rho = getMD758ClustersCorr;
% ROI = getMD758Clusters;
% getCrossCorrClusterVoxelLevelInsideAALAnyCluster(ROI,rho,label);
% 
% label = 'SP758';
% rho = getSPATIALClustersCorr;
% ROI = getSPATIALClusters;
% getCrossCorrClusterVoxelLevelInsideAALAnyCluster(ROI,rho,label);

% plotDistributionCorrVoxelsInsideAClusterWholeBrain;

% plotDistributionCrossCorrrWholeBrainAnyCluster;

% label = 'SP758';
% rho = getSPATIALClustersCorr;
% ROI = getSPATIALClusters;
% plotDistributionCorrClusterDiffCluSameAAL(ROI,rho,label);
% 
% label = 'MD758';
% rho = getMD758ClustersCorr;
% ROI = getMD758Clusters;
% plotDistributionCorrClusterDiffCluSameAAL(ROI,rho,label);

% label = 'SP758';
% rho = getSPATIALClustersCorr;
% ROI = getSPATIALClusters;
% plotDistributionCorrClusterDiffCluSameAALOnlyBiggest(ROI,rho,label);
% 
% label = 'MD758';
% rho = getMD758ClustersCorr;
% ROI = getMD758Clusters;
% plotDistributionCorrClusterDiffCluSameAALOnlyBiggest(ROI,rho,label);
 
% label = 'SP758';
% rho = getSPATIALClustersCorr;
% ROI = getSPATIALClusters;
% plotDistributionCorrClusterDiffCluDiffAAL(ROI,rho,label);
% 
% label = 'MD758';
% rho = getMD758ClustersCorr;
% ROI = getMD758Clusters;
% plotDistributionCorrClusterDiffCluDiffAAL(ROI,rho,label);

% label = 'SP758';
% [rho_rest, rho_passive, rho_track] = getSPATIALClustersCorrAttPassRest;
% ROI_SP = getSPATIALClusters;
% 
% label = 'SP758-Rest';
% plotDistributionCorrClusterDiffCluSameAAL(ROI_SP,rho_rest,label);
% plotDistributionCorrClusterDiffCluSameAALOnlyBiggest(ROI_SP,rho_rest,label);
% plotDistributionCorrClusterDiffCluDiffAAL(ROI_SP,rho_rest,label);
% 
% label = 'SP758-Passive';
% plotDistributionCorrClusterDiffCluSameAAL(ROI_SP,rho_passive,label);
% plotDistributionCorrClusterDiffCluSameAALOnlyBiggest(ROI_SP,rho_passive,label);
% plotDistributionCorrClusterDiffCluDiffAAL(ROI_SP,rho_passive,label);
% 
% label = 'SP758-Track';
% plotDistributionCorrClusterDiffCluSameAAL(ROI_SP,rho_track,label);
% plotDistributionCorrClusterDiffCluSameAALOnlyBiggest(ROI_SP,rho_track,label);
% plotDistributionCorrClusterDiffCluDiffAAL(ROI_SP,rho_track,label);
% 
% label = 'MD758';
% [rho_rest, rho_passive, rho_track] = getMD758ClustersCorrAttPassRest;
% ROI_MD = getMD758Clusters;
% 
% label = 'MD758-Rest';
% plotDistributionCorrClusterDiffCluSameAAL(ROI_MD,rho_rest,label);
% plotDistributionCorrClusterDiffCluSameAALOnlyBiggest(ROI_MD,rho_rest,label);
% plotDistributionCorrClusterDiffCluDiffAAL(ROI_MD,rho_rest,label);
% 
% label = 'MD758-Passive';
% plotDistributionCorrClusterDiffCluSameAAL(ROI_MD,rho_passive,label);
% plotDistributionCorrClusterDiffCluSameAALOnlyBiggest(ROI_MD,rho_passive,label);
% plotDistributionCorrClusterDiffCluDiffAAL(ROI_MD,rho_passive,label);
% 
% label = 'MD758-Track';
% plotDistributionCorrClusterDiffCluSameAAL(ROI_MD,rho_track,label);
% plotDistributionCorrClusterDiffCluSameAALOnlyBiggest(ROI_MD,rho_track,label);
% plotDistributionCorrClusterDiffCluDiffAAL(ROI_MD,rho_track,label);

% parcellation_label = 'MD758';
% [rho_rest, rho_passive, rho_track] = getMD758ClustersCorrAttPassRest;
% ROI_MD = getMD758Clusters;
% 
% condition_label = 'RestingState';
% getCrossCorrClusterVoxelLevelInsideAALAnyCluster(ROI_MD,rho_rest,parcellation_label,condition_label);
% getCrossCorrClusterVoxelLevelInsideAAL(ROI_MD,rho_rest,parcellation_label,condition_label);
% getCorrVoxelsInsideAClusterV2(ROI_MD,parcellation_label,condition_label);
% 
% condition_label = 'Passive';
% getCrossCorrClusterVoxelLevelInsideAALAnyCluster(ROI_MD,rho_passive,parcellation_label,condition_label);
% getCrossCorrClusterVoxelLevelInsideAAL(ROI_MD,parcellation_label,condition_label);
% getCrossCorrClusterVoxelLevelInsideAAL(ROI_MD,rho_passive,parcellation_label,condition_label);
% getCorrVoxelsInsideAClusterV2(ROI_MD,parcellation_label,condition_label);
% 
% condition_label = 'Track';
% getCrossCorrClusterVoxelLevelInsideAALAnyCluster(ROI_MD,rho_track,parcellation_label,condition_label);
% getCrossCorrClusterVoxelLevelInsideAAL(ROI_MD,parcellation_label,condition_label);
% getCrossCorrClusterVoxelLevelInsideAAL(ROI_MD,rho_track,parcellation_label,condition_label);
% getCorrVoxelsInsideAClusterV2(ROI_MD,parcellation_label,condition_label);

% parcellation_label = 'SP758';
% [rho_rest, rho_passive, rho_track] = getSPATIALClustersCorrAttPassRest;
% ROI_SP = getSPATIALClusters;
% 
% condition_label = 'RestingState';
% getCrossCorrClusterVoxelLevelInsideAALAnyCluster(ROI_SP,rho_rest,parcellation_label,condition_label);
% getCrossCorrClusterVoxelLevelInsideAAL(ROI_SP,rho_rest,parcellation_label,condition_label);
% getCorrVoxelsInsideAClusterV2(ROI_SP,parcellation_label,condition_label);
% 
% condition_label = 'Passive';
% getCrossCorrClusterVoxelLevelInsideAALAnyCluster(ROI_SP,rho_passive,parcellation_label,condition_label);
% getCrossCorrClusterVoxelLevelInsideAAL(ROI_SP,rho_passive,parcellation_label,condition_label);
% getCorrVoxelsInsideAClusterV2(ROI_SP,parcellation_label,condition_label);
% 
% condition_label = 'Track';
% getCrossCorrClusterVoxelLevelInsideAALAnyCluster(ROI_SP,rho_track,parcellation_label,condition_label);
% getCrossCorrClusterVoxelLevelInsideAAL(ROI_SP,rho_track,parcellation_label,condition_label);
% getCorrVoxelsInsideAClusterV2(ROI_SP,parcellation_label,condition_label);

% condition_label = 'RestingState';
% plotDistributionCorrVoxelsInsideAClusterWholeBrain(condition_label); DONE
% plotDistributionCrossCorrrWholeBrainAnyCluster(condition_label); DONE
% plotDistributionCrossCorrClusterVoxelLevelInsideAAL(condition_label);
% DONE

% condition_label = 'Passive';

% plotDistributionCorrVoxelsInsideAClusterWholeBrain(condition_label); DONE
% plotDistributionCrossCorrrWholeBrainAnyCluster(condition_label); DONE
% plotDistributionCrossCorrClusterVoxelLevelInsideAAL(condition_label);

% condition_label = 'Track';
% plotDistributionCorrVoxelsInsideAClusterWholeBrain(condition_label); DONE
% plotDistributionCrossCorrrWholeBrainAnyCluster(condition_label); DONE
% plotDistributionCrossCorrClusterVoxelLevelInsideAAL(condition_label);

% MyContourCorrelation_v5_redo;

% [rho_rest, rho_passive, rho_track] = getMD758ClustersCorrAttPassRest;
% plotContourMD758Conditions(rho_rest,'Rest');
% plotContourMD758Conditions(rho_passive,'Passive');
% plotContourMD758Conditions(rho_track,'Track');

% parcellation_label = 'MD758';
% condition_label = 'RestingState';
% ROI_clusters = getMD758Clusters;
% getMeasuresOfClusterValidityMD758(condition_label,parcellation_label,ROI_clusters);
% plotContourMeasuresOfClusterValidity(parcellation_label,condition_label);

% parcellation_label = 'SP758';
% condition_label = 'RestingState';
% ROI_clusters = getSPATIALClusters;
% getMeasuresOfClusterValiditySP758(condition_label,parcellation_label,ROI_clusters);
% plotContourMeasuresOfClusterValidity(parcellation_label,condition_label);

nSUBJ = 8;
for SUBJ_ID=2:nSUBJ
    
    parcellation_label = 'MD758';
    condition_label = 'RestingState';
    ROI_clusters = getMD758Clusters;
    getMeasuresOfClusterValidityPerSUBJMD758(condition_label,parcellation_label,ROI_clusters,SUBJ_ID);

end

for SUBJ_ID=1:nSUBJ
    
    parcellation_label = 'SP758';
    condition_label = 'RestingState';
    ROI_clusters = getSPATIALClusters;
    getMeasuresOfClusterValidityPerSUBJSP758(condition_label,parcellation_label,ROI_clusters,SUBJ_ID);
    
end

end

%%% COMPARE MD758 with SPATIAL758

function RestingState = getMagdeburgDataRS(file)

nRuns = 4;

all_settings = getAllSettings;
nSubjects = length(all_settings);

iiRun = 0;
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    % file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
    
    kind = 'RestingState';
    for irun=1:nRuns;
        iiRun = iiRun + 1;
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
end

end

function ROI = getSPATIALClusters

% load_aal = nifti('ROI_MNI_V4.nii');
% load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
% AAL_img = load_aal.dat(:,:,:);
% load_roi = load('ROI_MNI_V4_List.mat');
% AAL_ROI = load_roi.ROI;
% 
% nROI = 90;
% 
% for iROI=1:nROI
%     
%     ROI_label = strrep(AAL_ROI(iROI).Nom_L,'_','-');
%     
%     load(strcat('Clusters-Voxels-Per-3D-Coordinate-v2-',ROI_label,'.mat'));
%     
%     ROI(iROI).clusters = clusters;
%     
%     clear clusters
%     
% end

load('Z:\_DATA\Parcellation\Spatial-Cluster\v2-758-clusters\Anatomical-Spatial-Clusters-v2a.mat');

end

function rho = getSPATIALClustersCorr

load('Z:\_DATA\Parcellation\Spatial-Cluster\v2-758-clusters\corr-clusters\Anatomical-Spatial-Clusters-v2-Corr.mat');

nTotalRuns = 32;
nTotalClusters = 758;

rho = zeros(nTotalRuns,nTotalClusters,nTotalClusters);

for iRun=1:nTotalRuns
    
    rho(iRun,:,:) = Spatial_Clusters(iRun).rest_corr(:,:);
    
end

end

function [rho_rest, rho_passive, rho_track] = getSPATIALClustersCorrAttPassRest

load('Z:\_DATA\Parcellation\Spatial-Cluster\v2-758-clusters\corr-clusters\Anatomical-Spatial-Clusters-v2-Corr.mat');

nTotalRuns = 32;
nTotalClusters = 758;

rho_rest = zeros(nTotalRuns,nTotalClusters,nTotalClusters);
rho_passive = zeros(nTotalRuns,nTotalClusters,nTotalClusters);
rho_track = zeros(nTotalRuns,nTotalClusters,nTotalClusters);

for iRun=1:nTotalRuns
    
    rho_rest(iRun,:,:) = Spatial_Clusters(iRun).rest_corr(:,:);
    rho_passive(iRun,:,:) = Spatial_Clusters(iRun).passive_corr(:,:);
    rho_track(iRun,:,:) = Spatial_Clusters(iRun).track_corr(:,:);
    
end

end

function ROI = getMD758Clusters

load('Z:\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat');

end

function rho = getMD758ClustersCorr

load('Z:\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-corrected.mat');

nTotalRuns = 32;
nTotalClusters = 758;

rho = zeros(nTotalRuns,nTotalClusters,nTotalClusters);

for iRun=1:nTotalRuns
    
    rho(iRun,:,:) = FC_clusters(iRun).rest_rho(:,:);
    
end

end

function [rho_rest, rho_passive, rho_track] = getMD758ClustersCorrAttPassRest

load('Z:\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-corrected.mat');

nTotalRuns = 32;
nTotalClusters = 758;

rho_rest = zeros(nTotalRuns,nTotalClusters,nTotalClusters);
rho_passive = zeros(nTotalRuns,nTotalClusters,nTotalClusters);
rho_track = zeros(nTotalRuns,nTotalClusters,nTotalClusters);

for iRun=1:nTotalRuns
    
    rho_rest(iRun,:,:) = FC_clusters(iRun).rest_rho(:,:);
    rho_passive(iRun,:,:) = FC_clusters(iRun).passive_rho(:,:);
    rho_track(iRun,:,:) = FC_clusters(iRun).track_rho(:,:);
    
end

end

function ROI_info = getROIInfoOnClusters

load('Z:\_DATA\Parcellation\758-Cluster\ROI_info.mat');

end

function getMeanClusters(run,ROI,label)

nROI = 90;
nTotalRuns = length(run);
nTR = 150;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for iRun=1:nTotalRuns

    disp(strcat('Run:',int2str(iRun)));
    
    iiCluster = 0;

    for iROI=1:nROI
        
        nClusters = length(ROI(iROI).clusters);
        
        for iCluster=1:nClusters
            
            iiCluster = iiCluster + 1;

            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);

            area_RestingState = zeros(nVoxels,nTR);

            %izr = 0;

            %zeros_RestingState = [];

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

                area_RestingState(iVoxel,:) = run(iRun).run(idxx,idxy,idxz,:);

                %if sum(area_RestingState(iVoxel,:)) == 0

                    %izr = izr + 1;
                    %zeros_RestingState(izr) = iVoxel;

                %end

            end

            %area_RestingState(zeros_RestingState,:) = [];

            mn_RestingState = nanmean(area_RestingState,1);

            mean_area_RestingState_voxel(iiCluster,:) = mn_RestingState(:);

        end
        
    end

    mean_run(iRun).rest = mean_area_RestingState_voxel;

    clear mean_area_RestingState_voxel

end 

save(strcat('Clusters-ROI-mean-run-',label,'.mat'),'mean_run');

end

function getCorrClusters(label)

load(strcat('Clusters-ROI-mean-run-',label,'.mat'));

nTotalRuns = length(mean_run);

for iRun=1:nTotalRuns
    
   run = mean_run(iRun).rest;
   run = run';
   rho(iRun).rest_corr = corr(run);
    
end

save(strcat('Clusters-ROI-mean-run-corr-',label,'.mat'),'rho');

end

function all_runs = getCorrVoxelLevel(iROI,condition_label)

ROI = getMD758Clusters;

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\_ATLAS\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

nROI = 90;
nSubjects = 8;
nRuns = 4;
   
ROI_label = strrep(AAL_ROI(iROI).Nom_L,'_','-');

nVoxels = ROI(iROI).nVoxels;

all_runs = zeros(nRuns*nSubjects,nVoxels,nVoxels);

iiRun = 0;
for iSubject=1:nSubjects
    
    disp(strcat('iSubject:',int2str(iSubject)));

    load(strcat('LHR-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-',condition_label,'-',ROI_label,'.mat'));

    for iRun=1:nRuns

        iiRun = iiRun + 1;

        eval(strcat('all_runs(iiRun,:,:) = FC_Voxels.run(iRun).rho_',condition_label,'(:,:);'));

    end

    clear FC_Voxels

end

end

function plotMeanFCClusters(label)

load(strcat('Clusters-ROI-mean-run-corr-',label,'.mat'));

nROI = 758;
nTotalRuns = length(rho);
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                pairs_rest_rho(iRun,iPair) = single_rho;
            
            end
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z,1);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ mean(pairs_rest_z,1);



mean_corr = zeros(nROI,nROI);
s_corr = zeros(nROI,nROI);

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            %if my_cv_aal(iPair) < cv_criterion & my_cv_aal(iPair) > 0 & p_pairs(iPair) < p_criterion
                
                mean_corr(iROI,iiROI) =  m_pairs_z(iPair);
                s_corr(iROI,iiROI) =  s_pairs_z(iPair);
                
            %end
            
        end
        
    end
    
end

plot_label = strcat('FC-Clusters-',label);
plotFCgeneral(mean_corr,s_corr,plot_label);

end

function plotVarianceClusters(label)

load(strcat('Clusters-ROI-mean-run-corr-',label,'.mat'));

nROI = 758;
nTotalRuns = length(rho);

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            end
            
        end
        
    end
    
end

my_var_aal = var(pairs_rest);
my_cv_aal = std(pairs_rest) ./ abs(mean(pairs_rest,1));
my_mean_aal = mean(pairs_rest,1);
my_std_aal = std(pairs_rest);

vector_one = my_std_aal;
vector_two = my_mean_aal;
vector_three = my_cv_aal;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = strcat('Clusters-',label);

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

fraction = FractionSignificantConsistent( vector_two, vector_one, z_threshold, cv_threshold )

end

function plotContourMD758Conditions(rho,label)

nROI = 758;
nTotalRuns = length(rho);

rho_z = (1/2).*log((1 + rho)./(1 - rho));

my_mean = mean(rho_z,1);
my_std = std(rho_z,0,1);

title_label = strcat('Clusters-',label);

MyContourCorrelation_v5(my_mean(:),my_std(:),title_label,1);

end

%%% Best Pairings
function getCrossCorrClusterVoxelLevelInsideAAL(ROI,rho,parcellation_label,condition_label)

nROI = 90;
nTotalClusters = 758;
nTotalRuns = 32;

%%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ROI_info = getROIInfoOnClusters;

rho_z = (1/2).*log((1 + rho)./(1 - rho));
m_rho_z = squeeze(mean(rho_z,1));

iiCluster = 0;
% iiCluster = 127;

for iROI=1:nROI
    % for iROI=9:9
    
    idx_voxels_this_ROI = find(AAL_img==AAL_ROI(iROI).ID);
    
    try
    
       ROI_label = ROI(iROI).label; 

       disp(ROI_label); 

       all_runs = getCorrVoxelLevel(iROI,condition_label);

       nVoxels = ROI(iROI).nVoxels;
       nClusters = ROI(iROI).nClusters;

       for iCluster=1:nClusters

           disp(strcat('iCluster:',int2str(iCluster)));
           
           idx_voxels_this_cluster = ROI(iROI).clusters(iCluster).idx_voxels;

           iiCluster = iiCluster + 1;

           vector = squeeze(m_rho_z(iiCluster,:));
           vector(isinf(vector)) = 0;
           vector(isnan(vector)) = 0;

           if iROI==1

                vector((nClusters+1):end) = 0;

           else

               vector(1:(ROI_info{iROI,4}-1)) = 0;
               vector((ROI_info{iROI,5}+1):end) = 0;

           end
           
           vector(iiCluster) = -2;

           [s,I] = sort(vector,'descend');

           idx_g_max_cross_corr_cluster = I(1);

           vector(vector==0) = [];

           [s,I] = sort(vector,'descend');

           idx_l_max_cross_corr_cluster = I(1);
           
           idx_voxels_cross_cluster = ROI(iROI).clusters(idx_l_max_cross_corr_cluster).idx_voxels;

           all_corrs = all_runs(1:nTotalRuns,ismember(idx_voxels_this_ROI,idx_voxels_this_cluster),ismember(idx_voxels_this_ROI,idx_voxels_cross_cluster));
           % all_corrs = all_corrs(:); %%% DO NOT VECTORIZE

           save(strcat('Cross-Corr-Voxel-BasedOnCluster','-',parcellation_label,'-',condition_label,'-',ROI(iROI).label,'-','cluster','-',int2str(iCluster),'-',int2str(iiCluster),'.mat'),'all_corrs','iCluster','iiCluster','idx_g_max_cross_corr_cluster','idx_l_max_cross_corr_cluster');

           clear all_corrs

       end

       clear final_voxel_amount all_corrs end_idx_voxel_cross_cluster end_idx_voxel_this_cluster start_idx_voxel_cross_cluster start_idx_voxel_this_cluster idx_l_max_cross_corr_cluster idx_g_max_cross_corr_cluster vector Cluster_info nClusters nVoxels all_runs ROI_label
   
    catch ME
        
        msg = getReport(ME);
        
        disp(msg);
        
        save('tmp_whole_memory.mat');
        
        return
        
    end
   
end

end

%%% All Pairings
function getCrossCorrClusterVoxelLevelInsideAALAnyCluster(ROI,rho,parcellation_label,condition_label)

nROI = 90;
nTotalClusters = 758;
nTotalRuns = 32;

%%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\_ATLAS\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ROI_info = getROIInfoOnClusters;

rho_z = (1/2).*log((1 + rho)./(1 - rho));
m_rho_z = squeeze(mean(rho_z,1));

iiCluster = 0;
% iiCluster = 127;

for iROI=1:nROI
    % for iROI=9:9
    
    idx_voxels_this_ROI = find(AAL_img==AAL_ROI(iROI).ID);
    
    try
    
       ROI_label = ROI(iROI).label; 

       disp(ROI_label); 

       all_runs = getCorrVoxelLevel(iROI,condition_label);

       nVoxels = ROI(iROI).nVoxels;
       nClusters = ROI(iROI).nClusters;

       for iCluster=1:nClusters

           disp(strcat('iCluster:',int2str(iCluster)));
           
           idx_voxels_this_cluster = ROI(iROI).clusters(iCluster).idx_voxels;

           iiCluster = iiCluster + 1;
           
           idx_g_cross_corr_cluster = ROI_info{iROI,4}:ROI_info{iROI,5};
           
           idx_g_cross_corr_cluster(idx_g_cross_corr_cluster==iiCluster) = 0;
           
           idx_l_cross_corr_cluster = find(idx_g_cross_corr_cluster~=0);
           
           idx_voxels_cross_cluster = [];
           
           for iiiCluster=idx_l_cross_corr_cluster
               
                tmp_idx_voxels = ROI(iROI).clusters(iiiCluster).idx_voxels(:);
           
                idx_voxels_cross_cluster = [idx_voxels_cross_cluster(:); tmp_idx_voxels(:)];
           
           end
           
           all_corrs = all_runs(1:nTotalRuns,ismember(idx_voxels_this_ROI,idx_voxels_this_cluster),ismember(idx_voxels_this_ROI,idx_voxels_cross_cluster));
           % all_corrs = all_corrs(:); %%% DO NOT VECTORIZE

           save(strcat('Cross-Corr-Voxel-BasedOnClusterAnyCluster','-',parcellation_label,'-',condition_label,'-',ROI(iROI).label,'-','cluster','-',int2str(iCluster),'-',int2str(iiCluster),'.mat'),'all_corrs','iCluster','iiCluster','idx_g_cross_corr_cluster','idx_l_cross_corr_cluster');

           clear all_corrs

       end

       clear final_voxel_amount all_corrs end_idx_voxel_cross_cluster end_idx_voxel_this_cluster start_idx_voxel_cross_cluster start_idx_voxel_this_cluster idx_l_max_cross_corr_cluster idx_g_max_cross_corr_cluster vector Cluster_info nClusters nVoxels all_runs ROI_label
   
    catch ME
        
        msg = getReport(ME);
        
        disp(msg);
        
        save('tmp_whole_memory.mat');
        
        return
        
    end
   
end

end

function all_over_all = computeDistributionCrossCorrClusterVoxelLevelInsideAAL(parcellation_label)

%%% Execute the function inside the folder with the
%%% Cross-Corr-Voxel-BasedOneCluster FILES !!!

files = dir(strcat('Cross-Corr-Voxel-BasedOnCluster','-',parcellation_label,'-','*'));

nFiles = length(files);

all_over_all = [];

for iFile=1:nFiles
    
    name = files(iFile).name;
    
    if ~strcmp('.',name) && ~strcmp('..',name)
       
        load(name);
        
        all_over_all = [all_over_all; all_corrs];
        
        clear all_corrs
        
    end
    
end

end

function plotDistributionCrossCorrClusterVoxelLevelInsideAALPerCluster

nROI = 90;
nTotalClusters = 758;
nTotalRuns = 32;

ROI_info = getROIInfoOnClusters;

iiCluster = 0;
for iROI=1:nROI
    
    nClusters = ROI_info{iROI,3};
    
    for iCluster=1:nClusters
        
        iiCluster = iiCluster + 1;
        
        parcellation_label = 'MD758';
        load(strcat('Cross-Corr-Voxel-BasedOnCluster','-',parcellation_label,'-',condition_label,'-',ROI_info{iROI,2},'-','cluster','-',int2str(iCluster),'-',int2str(iiCluster),'.mat'));

        all_over_all_MD = all_corrs;
        clear all_corrs
        
        parcellation_label = 'SP758';
        load(strcat('Cross-Corr-Voxel-BasedOnCluster','-',parcellation_label,'-',condition_label,'-',ROI_info{iROI,2},'-','cluster','-',int2str(iCluster),'-',int2str(iiCluster),'.mat'));

        all_over_all_SP = all_corrs;
        clear all_corrs
        
        all_over_all_MD_z = (1/2)*log((1 + all_over_all_MD)./(1 - all_over_all_MD));
        all_over_all_SP_z = (1/2)*log((1 + all_over_all_SP)./(1 - all_over_all_SP));
        
        m_all_over_all_MD = squeeze(nanmean(all_over_all_MD,1));
        m_all_over_all_SP = squeeze(nanmean(all_over_all_SP,1));
        
        s_all_over_all_MD = squeeze(nanstd(all_over_all_MD,0,1));
        s_all_over_all_SP = squeeze(nanstd(all_over_all_SP,0,1));
        
        m_all_over_all_MD = m_all_over_all_MD(:);
        m_all_over_all_SP = m_all_over_all_SP(:);
        
        s_all_over_all_MD = s_all_over_all_MD(:);
        s_all_over_all_SP = s_all_over_all_SP(:);

        MyContourCorrelation_v5(m_all_over_all_MD,s_all_over_all_MD,strcat('Cross-Cluster-Corr-MD-',int2str(iiCluster)),1);
        MyContourCorrelation_v5(m_all_over_all_SP,s_all_over_all_SP,strcat('Cross-Cluster-Corr-SP-',int2str(iiCluster)),1);
        
        plotDistributionDensity(m_all_over_all_MD,m_all_over_all_SP,'MD758','SP758',strcat('Cross-Cluster-Corr-MD-SP-',condition_label,'-',int2str(iiCluster)),'Correlations','Probability');

        close all
        
    end
    
end

end

function plotDistributionCrossCorrrWholeBrain

nROI = 90;
nTotalClusters = 758;
nTotalRuns = 32;

ROI_info = getROIInfoOnClusters;

all_over_all_MD = [];
all_over_all_SP = [];

iiCluster = 0;
for iROI=1:nROI
    
    nClusters = ROI_info{iROI,3};
    
    for iCluster=1:nClusters
        
        iiCluster = iiCluster + 1;
        
        parcellation_label = 'MD758';
        load(strcat('Cross-Corr-Voxel-BasedOnCluster','-',parcellation_label,'-',ROI_info{iROI,2},'-','cluster','-',int2str(iCluster),'-',int2str(iiCluster),'.mat'));
        
        [l,c,d] = size(all_corrs);
        all_corrs = reshape(all_corrs,[l,c*d]);
        
        all_over_all_MD = [all_over_all_MD,all_corrs];
        clear all_corrs
        
        parcellation_label = 'SP758';
        load(strcat('Cross-Corr-Voxel-BasedOnCluster','-',parcellation_label,'-',ROI_info{iROI,2},'-','cluster','-',int2str(iCluster),'-',int2str(iiCluster),'.mat'));

        [l,c,d] = size(all_corrs);
        all_corrs = reshape(all_corrs,[l,c*d]);
        
        all_over_all_SP = [all_over_all_SP,all_corrs];
        clear all_corrs
        
    end
    
end
        
all_over_all_MD_z = (1/2)*log((1 + all_over_all_MD)./(1 - all_over_all_MD));
all_over_all_SP_z = (1/2)*log((1 + all_over_all_SP)./(1 - all_over_all_SP));

m_all_over_all_MD = squeeze(nanmean(all_over_all_MD_z,1));
m_all_over_all_SP = squeeze(nanmean(all_over_all_SP_z,1));

s_all_over_all_MD = squeeze(nanstd(all_over_all_MD_z,0,1));
s_all_over_all_SP = squeeze(nanstd(all_over_all_SP_z,0,1));

MyContourCorrelation_v5(m_all_over_all_MD,s_all_over_all_MD,strcat('Cross-Cluster-Corr-MD-',int2str(iiCluster)),1);
MyContourCorrelation_v5(m_all_over_all_SP,s_all_over_all_SP,strcat('Cross-Cluster-Corr-SP-',int2str(iiCluster)),1);

plotDistributionDensity(m_all_over_all_MD,m_all_over_all_SP,'MD758','SP758',strcat('Cross-Cluster-Corr-MD-SP-',int2str(iiCluster)),'Correlations','Probability');

close all

end

function getCorrVoxelsInsideAClusterV2(ROI,parcellation_label,condition_label)

nROI = 90;
iiClusters = 0;
nTR = 150;
MNI_img = [91 109 91];
nTotalRuns = 32;
nRuns = 4;

all_settings = getAllSettings;
nSubjects = length(all_settings);

iiRun = 0;
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
        
    kind = condition_label;
    for irun=1:nRuns;
        iiRun = iiRun + 1;
        [RUN(iiRun).run, RUN(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
end

for iROI=1:nROI
    
    nClusters = length(ROI(iROI).clusters);
    
    for iCluster=1:nClusters
        
        iiClusters = iiClusters + 1;
            
            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
            
            nVoxels = length(idx_voxels);
            
            voxels_corr = zeros(nTotalRuns,nVoxels,nVoxels);
            
            for iRun=1:nTotalRuns
            
                voxels_ts = zeros(nTR,nVoxels);
            
                for iVoxel=1:nVoxels
               
                    [idxx,idxy,idxz] = ind2sub(MNI_img,idx_voxels(iVoxel));
                
                    voxels_ts(:,iVoxel) = squeeze(RUN(iRun).run(idxx,idxy,idxz,:));
                    
                end
                
                voxels_corr(iRun,:,:) = corr(voxels_ts);
                
            end
            
            save(strcat('Voxels-Inside-Cluster-Corr-',parcellation_label,'-',condition_label,'-',int2str(iiClusters),'.mat'),'voxels_corr');
        
    end

end

end

function plotDistributionCorrVoxelsInsideAClusterWholeBrain(condition_label)

nROI = 90;
nTotalClusters = 758;
nTotalRuns = 32;

ROI_info = getROIInfoOnClusters;

all_over_all_MD = [];
all_over_all_SP = [];

iiCluster = 0;
for iROI=1:nROI
    
    nClusters = ROI_info{iROI,3};
    
    for iCluster=1:nClusters
        
        iiCluster = iiCluster + 1;
        
        parcellation_label = 'MD758';
        load(strcat('Voxels-Inside-Cluster-Corr','-',parcellation_label,'-',condition_label,'-',int2str(iiCluster),'.mat'));
        
        [l,c,d] = size(voxels_corr);
        all_corrs = reshape(voxels_corr,[l,c*d]);
        
        all_over_all_MD = [all_over_all_MD,all_corrs];
        clear all_corrs
        
        parcellation_label = 'SP758';
        load(strcat('Voxels-Inside-Cluster-Corr','-',parcellation_label,'-',condition_label,'-',int2str(iiCluster),'.mat'));

        [l,c,d] = size(voxels_corr);
        all_corrs = reshape(voxels_corr,[l,c*d]);
        
        all_over_all_SP = [all_over_all_SP,all_corrs];
        clear all_corrs
        
    end
    
end
        
all_over_all_MD_z = (1/2)*log((1 + all_over_all_MD)./(1 - all_over_all_MD));
all_over_all_SP_z = (1/2)*log((1 + all_over_all_SP)./(1 - all_over_all_SP));

m_all_over_all_MD = squeeze(nanmean(all_over_all_MD_z,1));
m_all_over_all_SP = squeeze(nanmean(all_over_all_SP_z,1));

s_all_over_all_MD = squeeze(nanstd(all_over_all_MD_z,0,1));
s_all_over_all_SP = squeeze(nanstd(all_over_all_SP_z,0,1));

MyContourCorrelation_v5(m_all_over_all_MD,s_all_over_all_MD,strcat('Voxel-Inside-Cluster-Corr-MD-',condition_label),1);
MyContourCorrelation_v5(m_all_over_all_SP,s_all_over_all_SP,strcat('Voxel-Inside-Cluster-Corr-SP-',condition_label),1);

plotDistributionDensity(m_all_over_all_MD,m_all_over_all_SP,'MD758','SP758',strcat('Voxel-Inside-Cluster-Corr-MD-SP-',condition_label),'Correlations','Probability');

close all

end

function plotDistributionCrossCorrrWholeBrainAnyCluster(condition_label)

nROI = 90;
nTotalClusters = 758;
nTotalRuns = 32;

ROI_info = getROIInfoOnClusters;

all_over_all_MD = [];
all_over_all_SP = [];

iiCluster = 0;
for iROI=1:nROI
    
    nClusters = ROI_info{iROI,3};
    
    for iCluster=1:nClusters
        
        iiCluster = iiCluster + 1;
        
        parcellation_label = 'MD758';
        load(strcat('Cross-Corr-Voxel-BasedOnClusterAnyCluster','-',parcellation_label,'-',condition_label,'-',ROI_info{iROI,2},'-','cluster','-',int2str(iCluster),'-',int2str(iiCluster),'.mat'));
        
        if ~isempty(all_corrs)
            
            [l,c,d] = size(all_corrs);
            all_corrs = reshape(all_corrs,[l,c*d]);
        
            all_over_all_MD = [all_over_all_MD,all_corrs];
            clear all_corrs
            
        end
        
        parcellation_label = 'SP758';
        load(strcat('Cross-Corr-Voxel-BasedOnClusterAnyCluster','-',parcellation_label,'-',condition_label,'-',ROI_info{iROI,2},'-','cluster','-',int2str(iCluster),'-',int2str(iiCluster),'.mat'));

        if ~isempty(all_corrs)
            
            [l,c,d] = size(all_corrs);
            all_corrs = reshape(all_corrs,[l,c*d]);

            all_over_all_SP = [all_over_all_SP,all_corrs];
            clear all_corrs
        
        end
        
    end
    
end
        
all_over_all_MD_z = (1/2)*log((1 + all_over_all_MD)./(1 - all_over_all_MD));
all_over_all_SP_z = (1/2)*log((1 + all_over_all_SP)./(1 - all_over_all_SP));

m_all_over_all_MD = squeeze(nanmean(all_over_all_MD_z,1));
m_all_over_all_SP = squeeze(nanmean(all_over_all_SP_z,1));

s_all_over_all_MD = squeeze(nanstd(all_over_all_MD_z,0,1));
s_all_over_all_SP = squeeze(nanstd(all_over_all_SP_z,0,1));

MyContourCorrelation_v5(m_all_over_all_MD,s_all_over_all_MD,strcat('Cross-AnyCluster-Corr-MD-',condition_label),1);
MyContourCorrelation_v5(m_all_over_all_SP,s_all_over_all_SP,strcat('Cross-AnyCluster-Corr-SP-',condition_label),1);

plotDistributionDensity(m_all_over_all_MD,m_all_over_all_SP,'MD758','SP758',strcat('Cross-AnyCluster-Corr-MD-SP-',condition_label),'Correlations','Probability');

close all

end

function plotDistributionCrossCorrClusterVoxelLevelInsideAAL(condition_label)

nROI = 90;
nTotalClusters = 758;
nTotalRuns = 32;

ROI_info = getROIInfoOnClusters;

all_over_all_MD = [];
all_over_all_SP = [];

iiCluster = 0;
for iROI=1:nROI
    
    nClusters = ROI_info{iROI,3};
    
    for iCluster=1:nClusters
        
        iiCluster = iiCluster + 1;
        
        parcellation_label = 'MD758';
        load(strcat('Cross-Corr-Voxel-BasedOnCluster','-',parcellation_label,'-',condition_label,'-',ROI_info{iROI,2},'-','cluster','-',int2str(iCluster),'-',int2str(iiCluster),'.mat'));
        
        if ~isempty(all_corrs)
            
            [l,c,d] = size(all_corrs);
            all_corrs = reshape(all_corrs,[l,c*d]);
        
            all_over_all_MD = [all_over_all_MD,all_corrs];
            clear all_corrs
            
        end
        
        parcellation_label = 'SP758';
        load(strcat('Cross-Corr-Voxel-BasedOnCluster','-',parcellation_label,'-',condition_label,'-',ROI_info{iROI,2},'-','cluster','-',int2str(iCluster),'-',int2str(iiCluster),'.mat'));

        if ~isempty(all_corrs)
            
            [l,c,d] = size(all_corrs);
            all_corrs = reshape(all_corrs,[l,c*d]);

            all_over_all_SP = [all_over_all_SP,all_corrs];
            clear all_corrs
        
        end
        
    end
    
end
        
all_over_all_MD_z = (1/2)*log((1 + all_over_all_MD)./(1 - all_over_all_MD));
all_over_all_SP_z = (1/2)*log((1 + all_over_all_SP)./(1 - all_over_all_SP));

m_all_over_all_MD = squeeze(nanmean(all_over_all_MD_z,1));
m_all_over_all_SP = squeeze(nanmean(all_over_all_SP_z,1));

s_all_over_all_MD = squeeze(nanstd(all_over_all_MD_z,0,1));
s_all_over_all_SP = squeeze(nanstd(all_over_all_SP_z,0,1));

MyContourCorrelation_v5(m_all_over_all_MD,s_all_over_all_MD,strcat('Cross-Corr-VoxelLevelInsideAAL-MD-',condition_label),1);
MyContourCorrelation_v5(m_all_over_all_SP,s_all_over_all_SP,strcat('Cross-Corr-VoxelLevelInsideAAL-SP-',condition_label),1);

plotDistributionDensity(m_all_over_all_MD,m_all_over_all_SP,'MD758','SP758',strcat('Cross-Corr-VoxelLevelInsideAAL-MD-SP-',condition_label),'Correlations','Probability');

close all

end

function plotDistributionCorrClusterDiffCluSameAAL(ROI,rho,label)

nROI = 90;

ROI_info = getROIInfoOnClusters;

all_total_corr = [];

for iROI=1:nROI
    
    start_clu = ROI_info{iROI,4};
    end_clu = ROI_info{iROI,5};
    
    idx_global = start_clu:end_clu;
    
    nClusters = ROI(iROI).nClusters;
    
    for iCluster=1:nClusters
        
        idx_global_removed = idx_global;
        idx_global_removed(iCluster) = [];
    
        all_corr = rho(:,idx_global(iCluster),idx_global_removed(:));
        
        all_corr = squeeze(all_corr);
        
        all_total_corr = [all_total_corr, all_corr];
        
    end
    
end

all_over_all_z = (1/2)*log((1 + all_total_corr)./(1 - all_total_corr));

m_all_over_all = squeeze(nanmean(all_over_all_z,1));

s_all_over_all = squeeze(nanstd(all_over_all_z,0,1));

MyContourCorrelation_v5(m_all_over_all,s_all_over_all,strcat('Cross-DiffCluSameAAL-Corr-CluLevel-',label),1);


end

function plotDistributionCorrClusterDiffCluSameAALOnlyBiggest(ROI,rho,label)

nROI = 90;

ROI_info = getROIInfoOnClusters;

rho_z = (1/2).*log((1 + rho)./(1 - rho));
m_rho_z = squeeze(mean(rho_z,1));

all_total_corr = [];

for iROI=1:nROI
    
    start_clu = ROI_info{iROI,4};
    end_clu = ROI_info{iROI,5};
    
    idx_global = start_clu:end_clu;
    
    nClusters = ROI(iROI).nClusters;
    
    if nClusters > 1
        
        for iCluster=1:nClusters

            idx_global_removed = idx_global;
            idx_global_removed(iCluster) = [];

            m_all_corr = squeeze(m_rho_z(idx_global(iCluster),idx_global_removed(:)));

            [s,I] = sort(m_all_corr,'descend');
            idx_biggest = I(1);

            all_corr = rho(:,idx_global(iCluster),idx_global_removed(idx_biggest));

            all_corr = squeeze(all_corr);

            all_total_corr = [all_total_corr, all_corr];

        end
    
    end
    
end

all_over_all_z = (1/2)*log((1 + all_total_corr)./(1 - all_total_corr));

m_all_over_all = squeeze(nanmean(all_over_all_z,1));

s_all_over_all = squeeze(nanstd(all_over_all_z,0,1));

MyContourCorrelation_v5(m_all_over_all,s_all_over_all,strcat('Cross-DiffCluSameAAL-Corr-CluLevel-OnlyBiggest-',label),1);

end

function plotDistributionCorrClusterDiffCluDiffAAL(ROI,rho,label)

nROI = 90;

ROI_info = getROIInfoOnClusters;

all_total_corr = [];

for iROI=1:nROI
    
    start_clu = ROI_info{iROI,4};
    end_clu = ROI_info{iROI,5};
    
    idx_global = start_clu:end_clu;
    
    nClusters = ROI(iROI).nClusters;
    
    for iCluster=1:nClusters
        
        all_corr = rho(:,idx_global(iCluster),:);
        all_corr = squeeze(all_corr);
        
        all_corr(:,idx_global(:)) = [];
        
        all_total_corr = [all_total_corr, all_corr];
        
    end
    
end

all_over_all_z = (1/2)*log((1 + all_total_corr)./(1 - all_total_corr));

m_all_over_all = squeeze(nanmean(all_over_all_z,1));

s_all_over_all = squeeze(nanstd(all_over_all_z,0,1));

MyContourCorrelation_v5(m_all_over_all,s_all_over_all,strcat('Cross-DiffCluDiffAAL-Corr-CluLevel-',label),1);


end

function MyContourCorrelation_v5_redo

bintype = 'Quantile';

bins_labels = {'Cross-AnyCluster-Corr-MD-RestingState' ...
    'Cross-AnyCluster-Corr-MD-Passive' ...
    'Cross-AnyCluster-Corr-MD-Track' ... 
    'Cross-AnyCluster-Corr-SP-RestingState' ...
    'Cross-AnyCluster-Corr-SP-Passive' ...
    'Cross-AnyCluster-Corr-SP-Track' ...
    'Cross-Corr-VoxelLevelInsideAAL-MD-RestingState' ...
    'Cross-Corr-VoxelLevelInsideAAL-MD-Passive' ...
    'Cross-Corr-VoxelLevelInsideAAL-MD-Track' ... 
    'Cross-Corr-VoxelLevelInsideAAL-SP-RestingState' ...
    'Cross-Corr-VoxelLevelInsideAAL-SP-Passive' ...
    'Cross-Corr-VoxelLevelInsideAAL-SP-Track' ...
    'Voxel-Inside-Cluster-Corr-MD-RestingState' ...
    'Voxel-Inside-Cluster-Corr-MD-Passive' ...
    'Voxel-Inside-Cluster-Corr-MD-Track' ... 
    'Voxel-Inside-Cluster-Corr-SP-RestingState' ...
    'Voxel-Inside-Cluster-Corr-SP-Passive' ...
    'Voxel-Inside-Cluster-Corr-SP-Track' ...
    'Voxels-Diff-AAL-Rest' ...
    'Voxels-Diff-AAL-passive-MD758' ...
    'Voxels-Diff-AAL-track-MD758' ... 
    'Voxels-Diff-AAL-rest-SP758' ...
    'Voxels-Diff-AAL-passive-SP758' ...
    'Voxels-Diff-AAL-track-SP758'}; 

for iLabel=1:length(bins_labels)
    
    load(strcat('Contour-',bins_labels{iLabel},'-',bintype,'-bins.mat'),'x_mid','y_mid','D_xy');

    plotcontour(bins_labels{iLabel}, bintype, 'Mean-STD', x_mid, y_mid, D_xy);

    close all
    
end

end

function getMeasuresOfClusterValidityMD758(condition_label,parcellation_label,ROI)

disp(condition_label);

nROI = 90;
nSubjects = 8;
nRuns = 4;

CH = struct('voxels',[]);
SC = struct('voxels',[]);

iiCluster = 0;
for iROI=1:nROI
    
    iiVoxel_CH = 0;
    iiVoxel_SC = 0;
    
    disp(strcat('ROI:',int2str(iROI)));
    
    ROI_label = ROI(iROI).label;
    
    nVoxels_ROI = ROI(iROI).nVoxels;
    
    rho = zeros(nSubjects*nRuns,nVoxels_ROI,nVoxels_ROI);
    
    disp('loading FC_Voxels + FC_KMeans...');
    
    load(strcat('LHR-All-Subjects-FC-Voxels-AAL-ROI-corr-',condition_label,'-',ROI_label,'-KMeans'));
    
    IdxClusters = FC_KMeans.all_runs.IdxClusters;
    
    iiRun = 0;
    for iSubject=1:nSubjects
    
        load(strcat('LHR-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-',condition_label,'-',ROI_label,'.mat'));
    
        for iRun=1:nRuns
            
            iiRun = iiRun + 1;
            
            eval(strcat('rho(',int2str(iiRun),',:,:) = FC_Voxels.run(',int2str(iRun),').rho_',condition_label,'(:,:);'));
            
        end
        
        clear FC_Voxels
        
    end
    
    rho(rho==1) = NaN;
    
    rho_z = (1/2)*log((1 + rho)./(1 - rho));
    
    disp('getting Correlation Profile...');
    
    corrProfile = zeros(nVoxels_ROI,nSubjects*nRuns);
    
    rho_z(isnan(rho_z)) = 0;
    
    for iVoxel=1:nVoxels_ROI
        
        for iiRun=1:nRuns*nSubjects
            
            corrProfile(iVoxel,iiRun) = nanmean(squeeze(rho_z(iiRun,iVoxel,:)));
            
        end
        
    end
    
    corrProfile = corrProfile';
    
    disp('getting overall Similarity and Dissimilarity...');
    
    corr_coef = corr(corrProfile);
    
    corr_coef(find(corr_coef==1)) = 0.999;
    
    similarity_z = (1/2)*log((1 + corr_coef)./(1 - corr_coef));
    
    dissimilarity_z = 1 - similarity_z;
    
    similarity = corr_coef;
    
    dissimilarity = 1 - similarity;
    
    nClusters = ROI(iROI).nClusters;
    
    for iCluster=1:nClusters
        
        iiCluster = iiCluster + 1;
        
        disp(strcat('Cluster:',int2str(iiCluster)));
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
        
        nVoxels = length(idx_voxels);
     
        idx_voxels_AAL = find(IdxClusters==iCluster);
        idx_voxels_non_AAL = find(IdxClusters~=iCluster);
        
        CH_tmp = similarity_z(idx_voxels_AAL,idx_voxels_AAL);
        
        for iVoxel=1:nVoxels
            
            iiVoxel_CH = iiVoxel_CH + 1;
        
            CH(iROI).voxels(iiVoxel_CH) = mean(CH_tmp(iVoxel,:));
            
        end

        SC_same = dissimilarity(idx_voxels_AAL,idx_voxels_AAL);
        
        SC_same_mean = mean(SC_same,2);
        
        SC_adj = dissimilarity(idx_voxels_AAL,idx_voxels_non_AAL);
        
        SC_adj_mean = mean(SC_adj,2);
        
        for iVoxel=1:nVoxels
            
            iiVoxel_SC = iiVoxel_SC + 1;
          
            SC(iROI).voxels(iiVoxel_SC) = ( SC_adj_mean(iVoxel) - SC_same_mean(iVoxel) ) / max(SC_adj_mean(iVoxel),SC_same_mean(iVoxel));
            
        end
        
    end
    
    clear FC_KMeans
    
end

save(strcat('Similarity-',condition_label,'-',parcellation_label,'.mat'),'SC','CH');

end

function getMeasuresOfClusterValidityPerSUBJMD758(condition_label,parcellation_label,ROI,SUBJ_ID)

disp(condition_label);

nROI = 90;
nSubjects = 8;
nRuns = 4;
nSubjects = 1;

CH = struct('voxels',[]);
SC = struct('voxels',[]);

iiCluster = 0;
for iROI=1:nROI
    
    iiVoxel_CH = 0;
    iiVoxel_SC = 0;
    
    disp(strcat('ROI:',int2str(iROI)));
    
    ROI_label = ROI(iROI).label;
    
    nVoxels_ROI = ROI(iROI).nVoxels;
    
    rho = zeros(nRuns,nVoxels_ROI,nVoxels_ROI);
    
    disp('loading FC_Voxels + FC_KMeans...');
    
    load(strcat('LHR-All-Subjects-FC-Voxels-AAL-ROI-corr-',condition_label,'-',ROI_label,'-KMeans'));
    
    IdxClusters = FC_KMeans.all_runs.IdxClusters;
    
    iiRun = 0;
    for iSubject=SUBJ_ID
    
        load(strcat('LHR-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-',condition_label,'-',ROI_label,'.mat'));
    
        for iRun=1:nRuns
            
            iiRun = iiRun + 1;
            
            eval(strcat('rho(',int2str(iiRun),',:,:) = FC_Voxels.run(',int2str(iRun),').rho_',condition_label,'(:,:);'));
            
        end
        
        clear FC_Voxels
        
    end
    
    rho(rho==1) = NaN;
    
    rho_z = (1/2)*log((1 + rho)./(1 - rho));
    
    disp('getting Correlation Profile...');
    
    corrProfile = zeros(nVoxels_ROI,nSubjects*nRuns);
    
    rho_z(isnan(rho_z)) = 0;
    
    for iVoxel=1:nVoxels_ROI
        
        for iiRun=1:nRuns*nSubjects
            
            corrProfile(iVoxel,iiRun) = nanmean(squeeze(rho_z(iiRun,iVoxel,:)));
            
        end
        
    end
    
    corrProfile = corrProfile';
    
    disp('getting overall Similarity and Dissimilarity...');
    
    corr_coef = corr(corrProfile);
    
    corr_coef(find(corr_coef==1)) = 0.999;
    
    similarity_z = (1/2)*log((1 + corr_coef)./(1 - corr_coef));
    
    dissimilarity_z = 1 - similarity_z;
    
    similarity = corr_coef;
    
    dissimilarity = 1 - similarity;
    
    nClusters = ROI(iROI).nClusters;
    
    for iCluster=1:nClusters
        
        iiCluster = iiCluster + 1;
        
        disp(strcat('Cluster:',int2str(iiCluster)));
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
        
        nVoxels = length(idx_voxels);
     
        idx_voxels_AAL = find(IdxClusters==iCluster);
        idx_voxels_non_AAL = find(IdxClusters~=iCluster);
        
        CH_tmp = similarity_z(idx_voxels_AAL,idx_voxels_AAL);
        
        for iVoxel=1:nVoxels
            
            iiVoxel_CH = iiVoxel_CH + 1;
        
            CH(iROI).voxels(iiVoxel_CH) = mean(CH_tmp(iVoxel,:));
            
        end

        SC_same = dissimilarity(idx_voxels_AAL,idx_voxels_AAL);
        
        SC_same_mean = mean(SC_same,2);
        
        SC_adj = dissimilarity(idx_voxels_AAL,idx_voxels_non_AAL);
        
        SC_adj_mean = mean(SC_adj,2);
        
        for iVoxel=1:nVoxels
            
            iiVoxel_SC = iiVoxel_SC + 1;
          
            SC(iROI).voxels(iiVoxel_SC) = ( SC_adj_mean(iVoxel) - SC_same_mean(iVoxel) ) / max(SC_adj_mean(iVoxel),SC_same_mean(iVoxel));
            
        end
        
    end
    
    clear FC_KMeans
    
end

save(strcat('Similarity-',condition_label,'-',parcellation_label,'-',int2str(SUBJ_ID),'-.mat'),'SUBJ_ID','SC','CH');

end

function getMeasuresOfClusterValiditySP758(condition_label,parcellation_label,ROI)

disp(condition_label);

nROI = 90;
nSubjects = 8;
nRuns = 4;

CH = struct('voxels',[]);
SC = struct('voxels',[]);

iiCluster = 0;
for iROI=1:nROI
    
    iiVoxel_CH = 0;
    iiVoxel_SC = 0;
    
    disp(strcat('ROI:',int2str(iROI)));
    
    ROI_label = ROI(iROI).label;
    
    nVoxels_ROI = ROI(iROI).nVoxels;
    
    rho = zeros(nSubjects*nRuns,nVoxels_ROI,nVoxels_ROI);
    
    disp('loading FC_Voxels + FC_KMeans...');
    
    IdxClusters = getIdxClustersSP758(iROI);
    
    iiRun = 0;
    for iSubject=1:nSubjects
    
        load(strcat('LHR-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-',condition_label,'-',ROI_label,'.mat'));
    
        for iRun=1:nRuns
            
            iiRun = iiRun + 1;
            
            eval(strcat('rho(',int2str(iiRun),',:,:) = FC_Voxels.run(',int2str(iRun),').rho_',condition_label,'(:,:);'));
            
        end
        
        clear FC_Voxels
        
    end
    
    rho(rho==1) = NaN;
    
    rho_z = (1/2)*log((1 + rho)./(1 - rho));
    
    disp('getting Correlation Profile...');
    
    corrProfile = zeros(nVoxels_ROI,nSubjects*nRuns);
    
    rho_z(isnan(rho_z)) = 0;
    
    for iVoxel=1:nVoxels_ROI
        
        for iiRun=1:nRuns*nSubjects
            
            corrProfile(iVoxel,iiRun) = nanmean(squeeze(rho_z(iiRun,iVoxel,:)));
            
        end
        
    end
    
    corrProfile = corrProfile';
    
    disp('getting overall Similarity and Dissimilarity...');
    
    corr_coef = corr(corrProfile);
    
    corr_coef(find(corr_coef==1)) = 0.999;
    
    similarity_z = (1/2)*log((1 + corr_coef)./(1 - corr_coef));
    
    dissimilarity_z = 1 - similarity_z;
    
    similarity = corr_coef;
    
    dissimilarity = 1 - similarity;
    
    nClusters = ROI(iROI).nClusters;
    
    for iCluster=1:nClusters
        
        iiCluster = iiCluster + 1;
        
        disp(strcat('Cluster:',int2str(iiCluster)));
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
        
        nVoxels = length(idx_voxels);
     
        idx_voxels_AAL = find(IdxClusters==iCluster);
        idx_voxels_non_AAL = find(IdxClusters~=iCluster);
        
        CH_tmp = similarity_z(idx_voxels_AAL,idx_voxels_AAL);
        
        for iVoxel=1:nVoxels
            
            iiVoxel_CH = iiVoxel_CH + 1;
        
            CH(iROI).voxels(iiVoxel_CH) = mean(CH_tmp(iVoxel,:));
            
        end

        SC_same = dissimilarity(idx_voxels_AAL,idx_voxels_AAL);
        
        SC_same_mean = mean(SC_same,2);
        
        SC_adj = dissimilarity(idx_voxels_AAL,idx_voxels_non_AAL);
        
        SC_adj_mean = mean(SC_adj,2);
        
        for iVoxel=1:nVoxels
            
            iiVoxel_SC = iiVoxel_SC + 1;
          
            SC(iROI).voxels(iiVoxel_SC) = ( SC_adj_mean(iVoxel) - SC_same_mean(iVoxel) ) / max(SC_adj_mean(iVoxel),SC_same_mean(iVoxel));
            
        end
        
    end
    
    clear FC_KMeans
    
end

save(strcat('Similarity-',condition_label,'-',parcellation_label,'.mat'),'SC','CH');

end

function getMeasuresOfClusterValidityPerSUBJSP758(condition_label,parcellation_label,ROI,SUBJ_ID)

disp(condition_label);

nROI = 90;
nSubjects = 8;
nRuns = 4;
nSubjects = 1;

CH = struct('voxels',[]);
SC = struct('voxels',[]);

iiCluster = 0;
for iROI=1:nROI
    
    iiVoxel_CH = 0;
    iiVoxel_SC = 0;
    
    disp(strcat('ROI:',int2str(iROI)));
    
    ROI_label = ROI(iROI).label;
    
    nVoxels_ROI = ROI(iROI).nVoxels;
    
    rho = zeros(nRuns,nVoxels_ROI,nVoxels_ROI);
    
    disp('loading FC_Voxels + FC_KMeans...');
    
    IdxClusters = getIdxClustersSP758(iROI);
    
    iiRun = 0;
    for iSubject=SUBJ_ID
    
        load(strcat('LHR-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-',condition_label,'-',ROI_label,'.mat'));
    
        for iRun=1:nRuns
            
            iiRun = iiRun + 1;
            
            eval(strcat('rho(',int2str(iiRun),',:,:) = FC_Voxels.run(',int2str(iRun),').rho_',condition_label,'(:,:);'));
            
        end
        
        clear FC_Voxels
        
    end
    
    rho(rho==1) = NaN;
    
    rho_z = (1/2)*log((1 + rho)./(1 - rho));
    
    disp('getting Correlation Profile...');
    
    corrProfile = zeros(nVoxels_ROI,nSubjects*nRuns);
    
    rho_z(isnan(rho_z)) = 0;
    
    for iVoxel=1:nVoxels_ROI
        
        for iiRun=1:nRuns*nSubjects
            
            corrProfile(iVoxel,iiRun) = nanmean(squeeze(rho_z(iiRun,iVoxel,:)));
            
        end
        
    end
    
    corrProfile = corrProfile';
    
    disp('getting overall Similarity and Dissimilarity...');
    
    corr_coef = corr(corrProfile);
    
    corr_coef(find(corr_coef==1)) = 0.999;
    
    similarity_z = (1/2)*log((1 + corr_coef)./(1 - corr_coef));
    
    dissimilarity_z = 1 - similarity_z;
    
    similarity = corr_coef;
    
    dissimilarity = 1 - similarity;
    
    nClusters = ROI(iROI).nClusters;
    
    for iCluster=1:nClusters
        
        iiCluster = iiCluster + 1;
        
        disp(strcat('Cluster:',int2str(iiCluster)));
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
        
        nVoxels = length(idx_voxels);
     
        idx_voxels_AAL = find(IdxClusters==iCluster);
        idx_voxels_non_AAL = find(IdxClusters~=iCluster);
        
        CH_tmp = similarity_z(idx_voxels_AAL,idx_voxels_AAL);
        
        for iVoxel=1:nVoxels
            
            iiVoxel_CH = iiVoxel_CH + 1;
        
            CH(iROI).voxels(iiVoxel_CH) = mean(CH_tmp(iVoxel,:));
            
        end

        SC_same = dissimilarity(idx_voxels_AAL,idx_voxels_AAL);
        
        SC_same_mean = mean(SC_same,2);
        
        SC_adj = dissimilarity(idx_voxels_AAL,idx_voxels_non_AAL);
        
        SC_adj_mean = mean(SC_adj,2);
        
        for iVoxel=1:nVoxels
            
            iiVoxel_SC = iiVoxel_SC + 1;
          
            SC(iROI).voxels(iiVoxel_SC) = ( SC_adj_mean(iVoxel) - SC_same_mean(iVoxel) ) / max(SC_adj_mean(iVoxel),SC_same_mean(iVoxel));
            
        end
        
    end
    
    clear FC_KMeans
    
end

save(strcat('Similarity-',condition_label,'-',parcellation_label,'-',int2str(SUBJ_ID),'-.mat'),'SUBJ_ID','SC','CH');

end

function IdxClusters = getIdxClustersSP758(iROI)

load('Z:\_DATA\Parcellation\Spatial-Cluster\v2-758-clusters\Anatomical-Spatial-Clusters-v2a.mat');

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\_ATLAS\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ROI_ID = AAL_ROI(iROI).ID;

idx_voxels = find(AAL_img == ROI_ID);

IdxClusters = zeros(size(idx_voxels));

nClusters = length(ROI(iROI).clusters);

for iCluster=1:nClusters
   
    idx_voxels_cluster = ROI(iROI).clusters(iCluster).idx_voxels;
    
    tmp_member = ismember(idx_voxels,idx_voxels_cluster);
    
    IdxClusters(find(tmp_member)) = iCluster;
    
end


end

function plotContourMeasuresOfClusterValidity(parcellation_label,condition_label)

load(strcat('Similarity-',condition_label,'-',parcellation_label,'.mat'));

nROI = 90;

for iROI=1:nROI
    
    CH_mean(iROI) = mean(CH(iROI).voxels(:));
    CH_std(iROI) = std(CH(iROI).voxels(:));
    
    SC_mean(iROI) = mean(SC(iROI).voxels(:));
    SC_std(iROI) = std(SC(iROI).voxels(:));
    
end
    
% rho_z = (1/2).*log((1 + rho)./(1 - rho));

% title_label = strcat('CH-',condition_label,'-',parcellation_label);
% MyContourCorrelation_v5(CH_mean(:),CH_std(:),title_label,1);
% 
% title_label = strcat('SC-',condition_label,'-',parcellation_label);
% MyContourCorrelation_v5(SC_mean(:),SC_std(:),title_label,1);

title_label = strcat('CH-SC-',condition_label,'-',parcellation_label);
MyContourCorrelation_v5(CH_mean(:),SC_mean(:),title_label,1);

end






