function real_methods_paper_v1_functional

% getCorrClustersDifferentClusterSameDifferentAAL;
 
% label = 'TOTAL';
% nStartCluster = 1;
% nEndCluster = 758;
% plotMeanFCFunctionalClusters(nStartCluster,nEndCluster,label);

% totalInformationContent;

% plotVarianceDistributionClustersDifferentClusterEverybody;
% plotVarianceDistributionClustersDifferentClusterSameDifferentAAL;

% getCorrClustersDifferentClusterEverybodyAALHCP;

end

%%% FUNCTIONAL CLUSTER

function getCorrClustersDifferentClusterSameDifferentAAL

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat');

nClusters = 758;
nROI = 90;
nRuns = 32;

for iRun=1:nRuns
    
    iiCluster = 0;
    
    for iROI=1:nROI
    
        nClusters = length(ROI(iROI).clusters);
        
        for iCluster=1:nClusters
      
            iiCluster = iiCluster + 1;
            
            func_clusters(iRun).rest_run(iiCluster,:) = ROI(iROI).clusters(iCluster).rest_mean(iRun,:);
            func_clusters(iRun).passive_run(iiCluster,:) = ROI(iROI).clusters(iCluster).passive_mean(iRun,:);
            func_clusters(iRun).track_run(iiCluster,:) = ROI(iROI).clusters(iCluster).track_mean(iRun,:);
            
        end
        
    end
    
end

for iRun=1:nRuns
    
   func_clusters(iRun).rest_corr = corr(func_clusters(iRun).rest_run');
   func_clusters(iRun).passive_corr = corr(func_clusters(iRun).passive_run');
   func_clusters(iRun).track_corr = corr(func_clusters(iRun).track_run');
    
end
       
save('758-Functional-Clusters-Corr-corrected.mat','func_clusters');

end

function getCorrClustersPARCELSpetraDATAmd

load('Z:\Dropbox (Uni Magdeburg)\_DATA\PETRA-RITTER\all_subjects\FC-Voxels-AAL-Kmeans-Volume\PETRA-FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nTotalClusters = 758;
nROI = 90;
nTotalRuns = 32;
nRuns = 4;
nTR = 150;
MNI_size = [91 109 91];

all_settings = getAllSettings;

nSubjects = length(all_settings);

iiRun = 0;
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;

    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;

    kind = 'RestingState';
    for irun=1:nRuns
        iiRun = iiRun + 1;
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
end

for iRun=1:nTotalRuns
    
    disp(strcat('iRun:',int2str(iRun)));
    
    func_clusters(iRun).corr = zeros(nTotalClusters);
    
    this_run_clusters = zeros(nTR,nTotalClusters);
    
    iiCluster = 0;
    
    for iROI=1:nROI
    
        nClusters = length(ROI(iROI).clusters);
        
        for iCluster=1:nClusters
      
            iiCluster = iiCluster + 1;
            
            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
            
            nVoxels = length(idx_voxels);
            
            voxels_this_cluster = zeros(nTR,nVoxels);
            
            for iVoxel=1:nVoxels
                
                [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
                
                voxels_this_cluster(:,iVoxel) = RestingState(iRun).run(idxx,idxy,idxz,1:nTR);
                
            end
            
            m_voxels_this_cluster = mean(voxels_this_cluster,2);
            
            this_run_clusters(:,iiCluster) = m_voxels_this_cluster(:);
            
        end
        
    end
    
    func_clusters(iRun).corr = corr(this_run_clusters);
    
end
       
save('758-Functional-Clusters-Corr-PARCELS-PETRA-DATA-MD.mat','func_clusters');

end

function getCorrClustersDifferentClusterEverybodyAALPETRA

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nClusters = 758;
nROI = 90;
nRuns = 9*3;
MNI_size = [91 109 91];
nTR = 150;
nSegments = 3;

run = 1;
all_settings = getAllSettingsPetra;

iRun = 0;
for iSet=1:length(all_settings)
    
    settings = all_settings(iSet).settings;

    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;

    [residual, settings] = real_get_data_FSL_Petra(settings,run,file,get_at_this_preprocessed_step);
    
    for iSeg=1:nSegments
        
        iRun = iRun + 1;
        
        Subject(iRun).residual(:,:,:,:) = residual(:,:,:,(1+(iSeg-1)*nTR):(nTR+(iSeg-1)*nTR));
        
    end

end

for iRun=1:nRuns
               
    func_clusters(iRun).rest_run = zeros(nTR,nClusters);
    
    iiCluster = 0;
    for iROI=1:nROI

        nCluster = length(ROI(iROI).clusters);

        for iCluster=1:nCluster
            
            iiCluster = iiCluster + 1;

            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);
            
            this_cluster_voxels = zeros(nVoxels,nTR);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

                this_cluster_voxels(iVoxel,:) = squeeze(Subject(iRun).residual(idxx,idxy,idxz,:));

            end
            
            func_clusters(iRun).rest_run(:,iiCluster) = mean(this_cluster_voxels,1);

        end

    end

end


for iRun=1:nRuns
    
   func_clusters(iRun).rest_corr = corr(func_clusters(iRun).rest_run);
    
end
       
save('758-Functional-Clusters-Petra-Corr.mat','func_clusters','-v7.3');

end

function getCorrClustersDifferentClusterEverybodyAALPETRACustomParcellation

load('PETRA-FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nClusters = 758;
nROI = 90;
nRuns = 9*3;
MNI_size = [91 109 91];
nTR = 150;
nSegments = 3;

run = 1;
all_settings = getAllSettingsPetra;

iRun = 0;
for iSet=1:length(all_settings)
    
    settings = all_settings(iSet).settings;

    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;

    [residual, settings] = real_get_data_FSL_Petra(settings,run,file,get_at_this_preprocessed_step);
    
    for iSeg=1:nSegments
        
        iRun = iRun + 1;
        
        Subject(iRun).residual(:,:,:,:) = residual(:,:,:,(1+(iSeg-1)*nTR):(nTR+(iSeg-1)*nTR));
        
    end

end

for iRun=1:nRuns
               
    func_clusters(iRun).rest_run = zeros(nTR,nClusters);
    
    iiCluster = 0;
    for iROI=1:nROI

        nCluster = length(ROI(iROI).clusters);

        for iCluster=1:nCluster
            
            iiCluster = iiCluster + 1;

            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);
            
            this_cluster_voxels = zeros(nVoxels,nTR);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

                this_cluster_voxels(iVoxel,:) = squeeze(Subject(iRun).residual(idxx,idxy,idxz,:));

            end
            
            func_clusters(iRun).rest_run(:,iiCluster) = mean(this_cluster_voxels,1);

        end

    end

end


for iRun=1:nRuns
    
   func_clusters(iRun).rest_corr = corr(func_clusters(iRun).rest_run);
    
end
       
save('758-Functional-Clusters-Petra-Custom-Corr.mat','func_clusters','-v7.3');

end

function getCorrClustersDifferentClusterEverybodyAALHCP

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat');

nClusters = 758;
nROI = 90;
nTotalRuns = 32;
nRuns = 4;
MNI_size = [91 109 91];
nTR = 1200;

all_settings = getAllSettingsHCP;
nSubjects = length(all_settings);

iiRun = 0;
for iSet=1:5
    
    settings = all_settings(iSet).settings;
    
    disp(strcat('Subj:',int2str(iSet)));

    for iRun=1:nRuns
        
        iiRun = iiRun + 1;
        
        [all_runs(iiRun).run, settings] = real_get_data_HCP(settings,iRun);
    
    end
    
end

i_global_Run = 0;

for iRun=1:20
    
    i_global_Run = i_global_Run + 1;
               
    func_clusters(i_global_Run).rest_run = zeros(nTR,nClusters);
    
    iiCluster = 0;
    for iROI=1:nROI

        nCluster = length(ROI(iROI).clusters);

        for iCluster=1:nCluster
            
            iiCluster = iiCluster + 1;

            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);
            
            this_cluster_voxels = zeros(nVoxels,nTR);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

                this_cluster_voxels(iVoxel,:) = squeeze(all_runs(iRun).run(idxx,idxy,idxz,:));

            end
            
            func_clusters(i_global_Run).rest_run(:,iiCluster) = mean(this_cluster_voxels,1);

        end

    end

end

clear all_runs

iiRun = 0;
for iSet=6:8
    
    settings = all_settings(iSet).settings;
    
    disp(strcat('Subj:',int2str(iSet)));

    for iRun=1:nRuns
        
        iiRun = iiRun + 1;
        
        [all_runs(iiRun).run, settings] = real_get_data_HCP(settings,iRun);
    
    end
    
end

for iRun=1:12
    
    i_global_Run = i_global_Run + 1;
               
    func_clusters(i_global_Run).rest_run = zeros(nTR,nClusters);
    
    iiCluster = 0;
    for iROI=1:nROI

        nCluster = length(ROI(iROI).clusters);

        for iCluster=1:nCluster
            
            iiCluster = iiCluster + 1;

            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);
            
            this_cluster_voxels = zeros(nVoxels,nTR);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

                this_cluster_voxels(iVoxel,:) = squeeze(all_runs(iRun).run(idxx,idxy,idxz,:));

            end
            
            func_clusters(i_global_Run).rest_run(:,iiCluster) = mean(this_cluster_voxels,1);

        end

    end

end

for iRun=1:nTotalRuns
    
   func_clusters(iRun).rest_corr = corr(func_clusters(iRun).rest_run);
    
end
       
save('758-Functional-Clusters-HCP-Corr.mat','func_clusters','-v7.3');

end

function getCorrClustersDifferentClusterEverybodyAALPETRAWholeSegment

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nClusters = 758;
nROI = 90;
nRuns = 9;
MNI_size = [91 109 91];
nTR = 450;

run = 1;
all_settings = getAllSettingsPetra;

iRun = 0;
for iSet=1:length(all_settings)
    
    settings = all_settings(iSet).settings;

    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;

    [residual, settings] = real_get_data_FSL_Petra(settings,run,file,get_at_this_preprocessed_step);

    iRun = iRun + 1;

    Subject(iRun).residual(:,:,:,:) = residual(:,:,:,:);
 
end

for iRun=1:nRuns
               
    func_clusters(iRun).rest_run = zeros(nTR,nClusters);
    
    iiCluster = 0;
    for iROI=1:nROI

        nCluster = length(ROI(iROI).clusters);

        for iCluster=1:nCluster
            
            iiCluster = iiCluster + 1;

            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);
            
            this_cluster_voxels = zeros(nVoxels,nTR);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

                this_cluster_voxels(iVoxel,:) = squeeze(Subject(iRun).residual(idxx,idxy,idxz,:));

            end
            
            func_clusters(iRun).rest_run(:,iiCluster) = mean(this_cluster_voxels,1);

        end

    end

end


for iRun=1:nRuns
    
   func_clusters(iRun).rest_corr = corr(func_clusters(iRun).rest_run);
    
end
       
save('758-Functional-Clusters-Petra-Whole-Segment-Corr.mat','func_clusters','-v7.3');

end

function plotVarianceDistributionClustersDifferentClusterSameDifferentAAL

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat');
load('758-Functional-Clusters-Corr-corrected.mat');

nROI = 90;
nRuns = 32;

%%% SAME ALL - REST

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        for iiCluster=1:nClusters
            
            iPair = iPair + 1;
        
            for iRun=1:nRuns

                single_rho = func_clusters(iRun).rest_corr(iiiCluster,iiiCluster-iCluster+iiCluster);
                
                pairs_inside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

            end
        
        end

    end
    
end

my_var_inside = var(pairs_inside);

my_cv_inside = std(pairs_inside) ./ abs(mean(pairs_inside,1));

my_mean_inside = mean(pairs_inside,1);
my_std_inside = std(pairs_inside);

vector_one = my_std_inside;
vector_two = my_mean_inside;
vector_three = my_cv_inside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Same-AAL-Rest-corrected';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

nROI = 90;
nRuns = 32;

%%% SAME ALL - PASSIVE

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        for iiCluster=1:nClusters
            
            iPair = iPair + 1;
        
            for iRun=1:nRuns

                single_rho = func_clusters(iRun).passive_corr(iiiCluster,iiiCluster-iCluster+iiCluster);
                
                pairs_inside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

            end
        
        end

    end
    
end

my_var_inside = var(pairs_inside);

my_cv_inside = std(pairs_inside) ./ abs(mean(pairs_inside,1));

my_mean_inside = mean(pairs_inside,1);
my_std_inside = std(pairs_inside);

vector_one = my_std_inside;
vector_two = my_mean_inside;
vector_three = my_cv_inside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Same-AAL-Passive-corrected';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

nROI = 90;
nRuns = 32;

%%% SAME ALL - TRACK

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        for iiCluster=1:nClusters
            
            iPair = iPair + 1;
        
            for iRun=1:nRuns

                single_rho = func_clusters(iRun).track_corr(iiiCluster,iiiCluster-iCluster+iiCluster);
                
                pairs_inside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

            end
        
        end

    end
    
end

my_var_inside = var(pairs_inside);

my_cv_inside = std(pairs_inside) ./ abs(mean(pairs_inside,1));

my_mean_inside = mean(pairs_inside,1);
my_std_inside = std(pairs_inside);

vector_one = my_std_inside;
vector_two = my_mean_inside;
vector_three = my_cv_inside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Same-AAL-Track-corrected';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

%%% DIFFERENT AAL - REST

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        ivCluster = 0;
        for iiROI=1:nROI
            
            nnClusters = length(ROI(iiROI).clusters);
            
            for iiCluster=1:nnClusters
                
                ivCluster = ivCluster + 1;
                
                if iROI ~= iiROI
                    
                    iPair = iPair + 1;
                    
                    for iRun=1:nRuns

                        single_rho = func_clusters(iRun).rest_corr(iiiCluster,ivCluster);
                        
                        pairs_outside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

my_var_outside = var(pairs_outside);

my_cv_outside = std(pairs_outside) ./ abs(mean(pairs_outside,1));

my_mean_outside = mean(pairs_outside,1);
my_std_outside = std(pairs_outside);

vector_one = my_std_outside;
vector_two = my_mean_outside;
vector_three = my_cv_outside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Different-AAL-Rest-corrected';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

%%% DIFFERENT AAL - PASSIVE

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        ivCluster = 0;
        for iiROI=1:nROI
            
            nnClusters = length(ROI(iiROI).clusters);
            
            for iiCluster=1:nnClusters
                
                ivCluster = ivCluster + 1;
                
                if iROI ~= iiROI
                    
                    iPair = iPair + 1;
                    
                    for iRun=1:nRuns
                        
                        single_rho = func_clusters(iRun).passive_corr(iiiCluster,ivCluster);
                        
                        pairs_outside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

my_var_outside = var(pairs_outside);

my_cv_outside = std(pairs_outside) ./ abs(mean(pairs_outside,1));

my_mean_outside = mean(pairs_outside,1);
my_std_outside = std(pairs_outside);

vector_one = my_std_outside;
vector_two = my_mean_outside;
vector_three = my_cv_outside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Different-AAL-Passive-corrected';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

%%% DIFFERENT AAL- TRACK

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        ivCluster = 0;
        for iiROI=1:nROI
            
            nnClusters = length(ROI(iiROI).clusters);
            
            for iiCluster=1:nnClusters
                
                ivCluster = ivCluster + 1;
                
                if iROI ~= iiROI
                    
                    iPair = iPair + 1;
                    
                    for iRun=1:nRuns
                        
                        single_rho = func_clusters(iRun).track_corr(iiiCluster,ivCluster);
                        
                        pairs_outside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

my_var_outside = var(pairs_outside);

my_cv_outside = std(pairs_outside) ./ abs(mean(pairs_outside,1));

my_mean_outside = mean(pairs_outside,1);
my_std_outside = std(pairs_outside);

vector_one = my_std_outside;
vector_two = my_mean_outside;
vector_three = my_cv_outside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Different-AAL-Track-corrected';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotVarianceDistributionClustersDifferentClusterEverybodyAALPETRA

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Petra-Corr.mat');

nROI = 90;
% nRuns = 27;
nRuns = 24;
nClusters = 758;

iPair = 0;
for iCluster=1:nClusters

    for iiCluster=1:nClusters

        iPair = iPair + 1;

        for iRun=1:nRuns

            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end

end

my_var = nanvar(pairs);

my_cv = nanstd(pairs) ./ abs(nanmean(pairs,1));

my_mean = nanmean(pairs,1);
my_std = nanstd(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Petra-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotVarCluPARCELSpetraDATAmd

load('758-Functional-Clusters-Corr-PARCELS-PETRA-DATA-MD.mat');

nTotalRuns = 32;
nTotalClusters = 758;

iPair = 0;
for iCluster=1:nTotalClusters

    for iiCluster=1:nTotalClusters

        iPair = iPair + 1;

        for iRun=1:nTotalRuns

            single_rho = func_clusters(iRun).corr(iCluster,iiCluster);

            pairs(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end

end

my_var = nanvar(pairs);

my_cv = nanstd(pairs) ./ abs(nanmean(pairs,1));

my_mean = nanmean(pairs,1);
my_std = nanstd(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Func-Parcel-Petra-Data-MD';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotFCCluPARCELSpetraDATAmd

load('758-Functional-Clusters-Corr-PARCELS-PETRA-DATA-MD.mat');

nTotalRuns = 32;
nTotalClusters = 758;

iPair = 0;
for iCluster=1:nTotalClusters

    for iiCluster=1:nTotalClusters

        iPair = iPair + 1;

        for iRun=1:nTotalRuns

            single_rho = func_clusters(iRun).corr(iCluster,iiCluster);

            pairs(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end

end

my_var = nanvar(pairs);

my_cv = nanstd(pairs) ./ abs(nanmean(pairs,1));

my_mean = nanmean(pairs,1);
my_std = nanstd(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Func-Parcel-Petra-Data-MD';

iPair = 0;
jCluster = 0;
for iCluster=1:nTotalClusters
    
    jCluster = jCluster + 1;
    
    jjCluster = 0;
    for iiCluster=1:nTotalClusters
        
        jjCluster = jjCluster + 1;
        
            iPair = iPair + 1;
            
            if iCluster ~= iiCluster
        
                %if my_cv_aal(iPair) < cv_criterion & p_pairs(iPair) < p_criterion

                    mean_corr(jCluster,jjCluster) =  my_mean(iPair);
                    s_corr(jCluster,jjCluster) =  my_std(iPair);

                %end
            
            end
            
    end
    
end

plotFCgeneral(mean_corr,s_corr,title_label);

end

function plotVarCluEveryAALPETRACustomParcellation

load('PETRA-FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');
load('758-Functional-Clusters-Petra-Custom-Corr.mat');

nROI = 90;
nRuns = 27;
nClusters = 758;

iPair = 0;
for iCluster=1:nClusters

    for iiCluster=1:nClusters

        iPair = iPair + 1;

        for iRun=1:nRuns

            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end

end

my_var = nanvar(pairs);

my_cv = nanstd(pairs) ./ abs(nanmean(pairs,1));

my_mean = nanmean(pairs,1);
my_std = nanstd(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Petra-Custom-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotVarianceDistributionClustersDifferentClusterEverybodyAALHCP

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-HCP-Corr.mat');

nROI = 90;
nRuns = 32;
nClusters = 758;

iPair = 0;
for iCluster=1:nClusters

    for iiCluster=1:nClusters

        iPair = iPair + 1;

        for iRun=1:nRuns

            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end

end

my_var = nanvar(pairs);

my_cv = nanstd(pairs) ./ abs(nanmean(pairs,1));

my_mean = nanmean(pairs,1);
my_std = nanstd(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-HCP-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotVarianceDistributionClustersDifferentClusterEverybodyAALPETRAWholeSegment

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Petra-Whole-Segment-Corr.mat');

nROI = 90;
nRuns = 9;
nClusters = 758;

iPair = 0;
for iCluster=1:nClusters

    for iiCluster=1:nClusters

        iPair = iPair + 1;

        for iRun=1:nRuns

            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end

end

my_var = nanvar(pairs);

my_cv = nanstd(pairs) ./ abs(nanmean(pairs,1));

my_mean = nanmean(pairs,1);
my_std = nanstd(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Petra-Whole-Segment-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotVarianceDistributionClustersDifferentClusterEverybody

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat');
load('758-Functional-Clusters-Corr-corrected.mat');

nTotalCluster = 758;
nRuns = 32;

%%% EVERYBODY - REST

iPair = 0;
for iCluster=1:nTotalCluster
   
    for iiCluster=1:nTotalCluster

        iPair = iPair + 1;

        for iRun=1:nRuns

            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end
    
end

my_var = var(pairs);

my_cv = std(pairs) ./ abs(mean(pairs,1));

my_mean = mean(pairs,1);
my_std = std(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Everybody-Rest-corrected';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

fraction = FractionSignificantConsistent( vector_two, vector_one, z_threshold, cv_threshold )

end

function plotVarianceDistributionClustersDifferentClusterEverybodyONESubject(nStartRun,nEndRun)

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Corr.mat');

nTotalCluster = 758;
nRuns = 4;

%%% EVERYBODY - REST

iPair = 0;
for iCluster=1:nTotalCluster
   
    for iiCluster=1:nTotalCluster

        iPair = iPair + 1;

        iiRun = 0;
        for iRun=nStartRun:nEndRun
            
            iiRun = iiRun + 1;

            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs(iiRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end
    
end

my_var = var(pairs);

my_cv = std(pairs) ./ abs(mean(pairs,1));

my_mean = mean(pairs,1);
my_std = std(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = strcat('Functional-Everybody-Rest','-Runs-',int2str(nStartRun),'-',int2str(nEndRun));

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

fraction = FractionSignificantConsistent( vector_two, vector_one, z_threshold, cv_threshold )

end

function plotMeanFCFunctionalClusters(nStartCluster,nEndCluster,label)

showAAL = 0;

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat');
load('758-Functional-Clusters-Corr-corrected.mat');

nROI = 90;
nRuns = 32;
nTotalClusters = 758;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

%%% ALL PAIRS - REST

iPair = 0;
for iCluster=nStartCluster:nEndCluster
    
    for iiCluster=nStartCluster:nEndCluster
        
        iPair = iPair + 1;
        
        for iRun=1:nRuns
            
            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            pairs_rest_rho(iRun,iPair) = single_rho;
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ abs(mean(pairs_rest_z,1));

% mean_corr = zeros(nEndCluster-nStartCluster + 1,nEndCluster-nStartCluster + 1);

iPair = 0;
jCluster = 0;
for iCluster=nStartCluster:nEndCluster
    
    jCluster = jCluster + 1;
    
    jjCluster = 0;
    for iiCluster=nStartCluster:nEndCluster
        
        jjCluster = jjCluster + 1;
        
            iPair = iPair + 1;
            
            if iCluster ~= iiCluster
        
                %if my_cv_aal(iPair) < cv_criterion & p_pairs(iPair) < p_criterion

                    mean_corr(jCluster,jjCluster) =  m_pairs_z(iPair);
                    s_corr(jCluster,jjCluster) =  s_pairs_z(iPair);

                %end
            
            end
            
    end
    
end

% f = figure;
% 
% hold on
% 
% clmap = colormap('jet');
% clmap(1,:) = [1 1 1];
% 
% min_C = 0;
% max_C = 1;
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% if showAAL
%     
%     before_Angular_L = 545;
%     after_Angular_L = 6;
%     
%     before_Frontal_Mid_L = 78;
%     after_Frontal_Mid_L = 25;
%     
%     rectangle('Position',[before_Angular_L before_Angular_L after_Angular_L after_Angular_L],'Curvature',0.2);
%     rectangle('Position',[before_Frontal_Mid_L before_Frontal_Mid_L after_Frontal_Mid_L after_Frontal_Mid_L],'Curvature',0.2);
% 
% end
% 
% print(f,'-depsc',strcat('FC-Functional-Clusters','-',label,'.eps'));

new_label = strcat('FC-Functional-Clusters-corrected','-',label);
plotFCgeneral(mean_corr,s_corr,new_label);


end

function plotMeanFCFunctionalClustersONESubject(nStartRun,nEndRun,nStartCluster,nEndCluster,label)

showAAL = 0;

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Corr.mat');

nROI = 90;
nRuns = 32;
nTotalClusters = 758;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

%%% ALL PAIRS - REST

iPair = 0;
for iCluster=nStartCluster:nEndCluster
    
    for iiCluster=nStartCluster:nEndCluster
        
        iPair = iPair + 1;
        
        iiRun = 0;
        for iRun=nStartRun:nEndRun
            
            iiRun = iiRun + 1;
            
            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs_rest_z(iiRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            pairs_rest_rho(iiRun,iPair) = single_rho;
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ abs(mean(pairs_rest_z,1));

% mean_corr = zeros(nEndCluster-nStartCluster + 1,nEndCluster-nStartCluster + 1);

iPair = 0;
jCluster = 0;
for iCluster=nStartCluster:nEndCluster
    
    jCluster = jCluster + 1;
    
    jjCluster = 0;
    for iiCluster=nStartCluster:nEndCluster
        
        jjCluster = jjCluster + 1;
        
            iPair = iPair + 1;
            
            if iCluster ~= iiCluster
        
                %if my_cv_aal(iPair) < cv_criterion & p_pairs(iPair) < p_criterion

                    mean_corr(jCluster,jjCluster) =  m_pairs_z(iPair);
                    s_corr(jCluster,jjCluster) =  s_pairs_z(iPair);

                %end
            
            end
            
    end
    
end

% f = figure;
% 
% hold on
% 
% clmap = colormap('jet');
% clmap(1,:) = [1 1 1];
% 
% min_C = 0;
% max_C = 1;
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% if showAAL
%     
%     before_Angular_L = 545;
%     after_Angular_L = 6;
%     
%     before_Frontal_Mid_L = 78;
%     after_Frontal_Mid_L = 25;
%     
%     rectangle('Position',[before_Angular_L before_Angular_L after_Angular_L after_Angular_L],'Curvature',0.2);
%     rectangle('Position',[before_Frontal_Mid_L before_Frontal_Mid_L after_Frontal_Mid_L after_Frontal_Mid_L],'Curvature',0.2);
% 
% end
% 
% print(f,'-depsc',strcat('FC-Functional-Clusters','-',label,'.eps'));

new_label = strcat('FC-Functional-Clusters','-',label,'-Runs-',int2str(nStartRun),'-',int2str(nEndRun));
plotFCgeneral(mean_corr,s_corr,new_label);


end

function plotMeanFCFunctionalClustersPETRA(nStartCluster,nEndCluster,label)

showAAL = 0;

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Petra-Corr.mat');

nROI = 90;
nRuns = 27;
nTotalClusters = 758;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

%%% ALL PAIRS - REST

iPair = 0;
for iCluster=nStartCluster:nEndCluster
    
    for iiCluster=nStartCluster:nEndCluster
        
        iPair = iPair + 1;
        
        for iRun=1:nRuns
            
            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            pairs_rest_rho(iRun,iPair) = single_rho;
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ abs(mean(pairs_rest_z,1));

% mean_corr = zeros(nEndCluster-nStartCluster + 1,nEndCluster-nStartCluster + 1);

iPair = 0;
jCluster = 0;
for iCluster=nStartCluster:nEndCluster
    
    jCluster = jCluster + 1;
    
    jjCluster = 0;
    for iiCluster=nStartCluster:nEndCluster
        
        jjCluster = jjCluster + 1;
        
            iPair = iPair + 1;
            
            if iCluster ~= iiCluster
        
                %if my_cv_aal(iPair) < cv_criterion & p_pairs(iPair) < p_criterion

                    mean_corr(jCluster,jjCluster) =  m_pairs_z(iPair);
                    s_corr(jCluster,jjCluster) =  s_pairs_z(iPair);

                %end
            
            end
            
    end
    
end

% f = figure;
% 
% hold on
% 
% clmap = colormap('jet');
% clmap(1,:) = [1 1 1];
% 
% min_C = 0;
% max_C = 1;
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% if showAAL
%     
%     before_Angular_L = 545;
%     after_Angular_L = 6;
%     
%     before_Frontal_Mid_L = 78;
%     after_Frontal_Mid_L = 25;
%     
%     rectangle('Position',[before_Angular_L before_Angular_L after_Angular_L after_Angular_L],'Curvature',0.2);
%     rectangle('Position',[before_Frontal_Mid_L before_Frontal_Mid_L after_Frontal_Mid_L after_Frontal_Mid_L],'Curvature',0.2);
% 
% end
% 
% print(f,'-depsc',strcat('FC-Functional-Clusters','-',label,'.eps'));

new_label = strcat('FC-Functional-Clusters-Petra','-',label);
plotFCgeneral(mean_corr,s_corr,new_label);


end

function plotMeanFCFunctionalClustersPETRAcustom(nStartCluster,nEndCluster,label)

showAAL = 0;

load('PETRA-FC-Voxels-AAL-ROI-corr-custom-KMeans-Info.mat');
load('758-Functional-Clusters-Petra-Custom-Corr.mat');

nROI = 90;
nRuns = 27;
nTotalClusters = 758;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

%%% ALL PAIRS - REST

iPair = 0;
for iCluster=nStartCluster:nEndCluster
    
    for iiCluster=nStartCluster:nEndCluster
        
        iPair = iPair + 1;
        
        for iRun=1:nRuns
            
            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            pairs_rest_rho(iRun,iPair) = single_rho;
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ abs(mean(pairs_rest_z,1));

% mean_corr = zeros(nEndCluster-nStartCluster + 1,nEndCluster-nStartCluster + 1);

iPair = 0;
jCluster = 0;
for iCluster=nStartCluster:nEndCluster
    
    jCluster = jCluster + 1;
    
    jjCluster = 0;
    for iiCluster=nStartCluster:nEndCluster
        
        jjCluster = jjCluster + 1;
        
            iPair = iPair + 1;
            
            if iCluster ~= iiCluster
        
                %if my_cv_aal(iPair) < cv_criterion & p_pairs(iPair) < p_criterion

                    mean_corr(jCluster,jjCluster) =  m_pairs_z(iPair);
                    s_corr(jCluster,jjCluster) =  s_pairs_z(iPair);

                %end
            
            end
            
    end
    
end

% f = figure;
% 
% hold on
% 
% clmap = colormap('jet');
% clmap(1,:) = [1 1 1];
% 
% min_C = 0;
% max_C = 1;
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% if showAAL
%     
%     before_Angular_L = 545;
%     after_Angular_L = 6;
%     
%     before_Frontal_Mid_L = 78;
%     after_Frontal_Mid_L = 25;
%     
%     rectangle('Position',[before_Angular_L before_Angular_L after_Angular_L after_Angular_L],'Curvature',0.2);
%     rectangle('Position',[before_Frontal_Mid_L before_Frontal_Mid_L after_Frontal_Mid_L after_Frontal_Mid_L],'Curvature',0.2);
% 
% end
% 
% print(f,'-depsc',strcat('FC-Functional-Clusters','-',label,'.eps'));

new_label = strcat('FC-Functional-Clusters-Petra-custom','-',label);
plotFCgeneral(mean_corr,s_corr,new_label);


end

function plotMeanFCFunctionalClustersHCP(nStartCluster,nEndCluster,label)

showAAL = 0;

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-HCP-Corr.mat');

nROI = 90;
nRuns = 32;
nTotalClusters = 758;
nTR = 1200;
cv_criterion = 0.75;
p_criterion = 0.01;

%%% ALL PAIRS - REST

iPair = 0;
for iCluster=nStartCluster:nEndCluster
    
    for iiCluster=nStartCluster:nEndCluster
        
        iPair = iPair + 1;
        
        for iRun=1:nRuns
            
            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            pairs_rest_rho(iRun,iPair) = single_rho;
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ abs(mean(pairs_rest_z,1));

% mean_corr = zeros(nEndCluster-nStartCluster + 1,nEndCluster-nStartCluster + 1);

iPair = 0;
jCluster = 0;
for iCluster=nStartCluster:nEndCluster
    
    jCluster = jCluster + 1;
    
    jjCluster = 0;
    for iiCluster=nStartCluster:nEndCluster
        
        jjCluster = jjCluster + 1;
        
            iPair = iPair + 1;
            
            if iCluster ~= iiCluster
        
                %if my_cv_aal(iPair) < cv_criterion & p_pairs(iPair) < p_criterion

                    mean_corr(jCluster,jjCluster) =  m_pairs_z(iPair);
                    s_corr(jCluster,jjCluster) =  s_pairs_z(iPair);

                %end
            
            end
            
    end
    
end

% f = figure;
% 
% hold on
% 
% clmap = colormap('jet');
% clmap(1,:) = [1 1 1];
% 
% min_C = 0;
% max_C = 1;
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% if showAAL
%     
%     before_Angular_L = 545;
%     after_Angular_L = 6;
%     
%     before_Frontal_Mid_L = 78;
%     after_Frontal_Mid_L = 25;
%     
%     rectangle('Position',[before_Angular_L before_Angular_L after_Angular_L after_Angular_L],'Curvature',0.2);
%     rectangle('Position',[before_Frontal_Mid_L before_Frontal_Mid_L after_Frontal_Mid_L after_Frontal_Mid_L],'Curvature',0.2);
% 
% end
% 
% print(f,'-depsc',strcat('FC-Functional-Clusters','-',label,'.eps'));

new_label = strcat('FC-Functional-Clusters-HCP','-',label);
plotFCgeneral(mean_corr,s_corr,new_label);


end

function plotMeanFCFunctionalClustersPETRAWholeSegment(nStartCluster,nEndCluster,label)

showAAL = 0;

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Petra-Whole-Segment-Corr.mat');

nROI = 90;
nRuns = 9;
nTotalClusters = 758;
nTR = 450;
cv_criterion = 0.75;
p_criterion = 0.01;

%%% ALL PAIRS - REST

iPair = 0;
for iCluster=nStartCluster:nEndCluster
    
    for iiCluster=nStartCluster:nEndCluster
        
        iPair = iPair + 1;
        
        for iRun=1:nRuns
            
            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            pairs_rest_rho(iRun,iPair) = single_rho;
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ abs(mean(pairs_rest_z,1));

% mean_corr = zeros(nEndCluster-nStartCluster + 1,nEndCluster-nStartCluster + 1);

iPair = 0;
jCluster = 0;
for iCluster=nStartCluster:nEndCluster
    
    jCluster = jCluster + 1;
    
    jjCluster = 0;
    for iiCluster=nStartCluster:nEndCluster
        
        jjCluster = jjCluster + 1;
        
            iPair = iPair + 1;
            
            if iCluster ~= iiCluster
        
                %if my_cv_aal(iPair) < cv_criterion & p_pairs(iPair) < p_criterion

                    mean_corr(jCluster,jjCluster) =  m_pairs_z(iPair);
                    s_corr(jCluster,jjCluster) =  s_pairs_z(iPair);

                %end
            
            end
            
    end
    
end

% f = figure;
% 
% hold on
% 
% clmap = colormap('jet');
% clmap(1,:) = [1 1 1];
% 
% min_C = 0;
% max_C = 1;
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% if showAAL
%     
%     before_Angular_L = 545;
%     after_Angular_L = 6;
%     
%     before_Frontal_Mid_L = 78;
%     after_Frontal_Mid_L = 25;
%     
%     rectangle('Position',[before_Angular_L before_Angular_L after_Angular_L after_Angular_L],'Curvature',0.2);
%     rectangle('Position',[before_Frontal_Mid_L before_Frontal_Mid_L after_Frontal_Mid_L after_Frontal_Mid_L],'Curvature',0.2);
% 
% end
% 
% print(f,'-depsc',strcat('FC-Functional-Clusters','-',label,'.eps'));

new_label = strcat('FC-Functional-Clusters-Petra-Whole-Segment','-',label);
plotFCgeneral(mean_corr,s_corr,new_label);


end

function totalInformationContentPerSubject

load('758-Functional-Clusters-Corr-corrected.mat');

nSegments = 4;
nSubjects = 8;
r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

iCorr = 0;
for iSubject=1:nSubjects

    for iSegment=1:nSegments
        
        iCorr = iCorr + 1;
        
        subject(iSubject).corr(iSegment,:,:) = func_clusters(iCorr).rest_corr;
        
    end
    
end

for iSubject=1:nSubjects

    single_rho = subject(iSubject).corr;

    subject(iSubject).fisher = (1/2).*log((1 + single_rho)./(1 - single_rho));
    
end
         
for iSubject=1:nSubjects
    
    my(iSubject).mean = squeeze(mean(subject(iSubject).fisher,1));
    my(iSubject).std = squeeze(std(subject(iSubject).fisher));      
    
end

for iSubject=1:nSubjects
    
    info(iSubject).IM_total = TotalInformationContent( my(iSubject).mean(:), my(iSubject).std(:), z_threshold, cv_threshold );
    
end

end

function totalInformationContentPerSubjectPetra

load('758-Functional-Clusters-Petra-Corr.mat');

nSegments = 3;
nSubjects = 9;
r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

iCorr = 0;
for iSubject=1:nSubjects

    for iSegment=1:nSegments
        
        iCorr = iCorr + 1;
        
        subject(iSubject).corr(iSegment,:,:) = func_clusters(iCorr).rest_corr;
        
    end
    
end

for iSubject=1:nSubjects

    single_rho = subject(iSubject).corr;

    subject(iSubject).fisher = (1/2).*log((1 + single_rho)./(1 - single_rho));
    
end
         
for iSubject=1:nSubjects
    
    my(iSubject).mean = squeeze(mean(subject(iSubject).fisher,1));
    my(iSubject).std = squeeze(std(subject(iSubject).fisher));      
    
end

for iSubject=1:nSubjects
    
    info(iSubject).IM_total = TotalInformationContent( my(iSubject).mean(:), my(iSubject).std(:), z_threshold, cv_threshold );
    
end

end

function totalInformationContent

load('758-Functional-Clusters-Corr-corrected.mat');

nTotalRuns = 32;
r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

for iRun=1:nTotalRuns

	all_corr(iRun,:,:) = func_clusters(iRun).rest_corr;
        
end

all_corr_fisher = (1/2).*log((1 + all_corr)./(1 - all_corr));
             
m_corr_z = squeeze(mean(all_corr_fisher,1));
s_corr_z = squeeze(std(all_corr_fisher,0,1));      

TotalInformationContent(m_corr_z(:), s_corr_z(:), z_threshold, cv_threshold );

end

