function real_methods_paper_v1_anatomical

% groupAnatomicalSpatialClusters;

% getCorrAnatomicalDifferentClusterSameDifferentAAL;

% plotVarianceDistributionAnatomicalDifferentClusterSameDifferentAAL;

totalInformationContent;

end

%%% ANATOMICAL CLUSTER

function getClustersOfVoxelsPer3DCoordinatesPerROI

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nROI = 90;
MNI_img = [91 109 91];
nTotalVoxels = 160990;
nTotalClusters = 758;
nClusterSize = 200;

iiVoxel = 0;

idx = zeros(3,nTotalVoxels);
all_idx_voxels = [];

for iROI=1:nROI
    
    disp(strcat('ROI:',int2str(iROI)));

    nClusters = length(ROI(iROI).clusters);
    
    nVoxels = 0;
    
    all_idx_voxels = [];
    for iCluster=1:nClusters
       
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
        
        nVoxels = nVoxels + length(idx_voxels);
        
        all_idx_voxels = [all_idx_voxels; idx_voxels];
        
    end
    
    if nClusters == 1
        
        ROI_clusters(iROI).clusters(1).idx_voxels = ROI(iROI).clusters(1).idx_voxels;
        
    else
        
        idx = zeros(3,nVoxels);
        
        for iVoxel=1:nVoxels
            
           iiVoxel = iiVoxel + 1;
            
           [idx(1,iiVoxel) idx(2,iiVoxel) idx(3,iiVoxel)] = ind2sub(MNI_img,all_idx_voxels(iVoxel)); 
            
        end
        
        D = pdist(idx','euclidean');

        closest_voxels = zeros(nVoxels,nClusterSize);

        for iVoxel=1:nVoxels
    
            %D((I-1)*(M-I/2)+J-I);
    
            distances = zeros(1,nVoxels);
    
            if mod(iVoxel,20000) == 0; disp(int2str(iVoxel)); end

            distances((iVoxel+1):nVoxels) = D((iVoxel-1)*(nVoxels-iVoxel/2) + ((iVoxel+1):nVoxels)-iVoxel);

            if iVoxel > 1

                distances((iVoxel-1):-1:1) = D((((iVoxel-1):-1:1)-1).*(nVoxels-((iVoxel-1):-1:1)./2)+iVoxel-((iVoxel-1):-1:1));
    
            end

            distances(iVoxel) = 0;

            [s,i] = sort(distances,'ascend');

            closest_voxels(iVoxel,:) = i(1:nClusterSize);
    
        end
        
        idx_clusters = kmeans(closest_voxels,nClusters);
        
        for iCluster=1:nClusters
   
            idx_voxels_on_cluster = find(idx_clusters==iCluster);
    
            ROI_clusters(iROI).clusters(iCluster).idx_voxels = all_idx_voxels(idx_voxels_on_cluster);
    
        end

    end
    
end

save('Anatomical-Spatial-Clusters.mat','ROI_clusters');

end

function groupAnatomicalSpatialClusters
        
nROIs = 90;
        
%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for iROI=1:nROIs
        
    cc = load(strcat('Clusters-Voxels-Per-3D-Coordinate-v2-',strrep(AAL_ROI(iROI).Nom_L,'_','-'),'.mat'));
    
    ROI_clusters(iROI).clusters = cc.clusters;
    
    clear cc
    
end

save('Anatomical-Spatial-Clusters-v2.mat','ROI_clusters');

end

function getCorrAnatomicalDifferentClusterSameDifferentAAL

load('Anatomical-Spatial-Clusters-v2.mat');

nRuns = 4;
nTotalRuns = 32;
nROI = 90;
nTR = 150;
nSubjects = 8;
nTotalClusters = 758;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

iiRun = 0;
all_settings = getAllSettings;
for iSubject=1:nSubjects
    settings = all_settings(iSubject).settings;
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
    
    kind = 'RestingState';
    for irun=1:nRuns;
        
        iiRun = iiRun + 1;
        
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    
    end
    
end

iiRun = 0;
all_settings = getAllSettings;
for iSubject=1:nSubjects
    settings = all_settings(iSubject).settings;
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;

    kind = 'Passive';
    for irun=1:nRuns;
        
        iiRun = iiRun + 1;
        
        [Passive(iiRun).run, Passive(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    
    end

end

iiRun = 0;
all_settings = getAllSettings;
for iSubject=1:nSubjects
    settings = all_settings(iSubject).settings;
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;

    kind = 'Track';
    for irun=1:nRuns;
        
        iiRun = iiRun + 1;
        
        [Track(iiRun).run, Track(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    
    end
    
end

resting_state_clusters = zeros(nTotalClusters,nTR);
passive_clusters = zeros(nTotalClusters,nTR);
track_clusters = zeros(nTotalClusters,nTR);

for iRun=1:nTotalRuns
    
    iiCluster = 0;
    for iROI=1:nROI

        nClusters = length(ROI_clusters(iROI).clusters);

        for iCluster=1:nClusters

            iiCluster = iiCluster + 1;

            idx_voxels = ROI_clusters(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);

            resting_state_voxels = zeros(nTR,nVoxels);
            passive_voxels = zeros(nTR,nVoxels);
            track_voxels = zeros(nTR,nVoxels);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

                resting_state_voxels(:,iVoxel) = RestingState(iRun).run(idxx,idxy,idxz,:);
                passive_voxels(:,iVoxel) = Passive(iRun).run(idxx,idxy,idxz,:);
                track_voxels(:,iVoxel) = Track(iRun).run(idxx,idxy,idxz,:);

            end

            resting_state_clusters(iiCluster,:) = mean(resting_state_voxels',1);
            passive_clusters(iiCluster,:) = mean(passive_voxels',1);
            track_clusters(iiCluster,:) = mean(track_voxels',1);

        end

    end
    
    Spatial_Clusters(iRun).rest_corr = corr(resting_state_clusters');
    Spatial_Clusters(iRun).passive_corr = corr(passive_clusters');
    Spatial_Clusters(iRun).track_corr = corr(track_clusters');
    
end

save('Anatomical-Spatial-Clusters-v2-Corr.mat','Spatial_Clusters');

end

function plotVarianceDistributionAnatomicalDifferentClusterSameDifferentAAL

load('Anatomical-Spatial-Clusters-v2.mat');
load('Anatomical-Spatial-Clusters-v2-Corr.mat');

nROI = 90;
nRuns = 32;

%%% SAME ALL - RESTING STATE

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI_clusters(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        for iiCluster=1:nClusters
            
            iPair = iPair + 1;
        
            for iRun=1:nRuns

                single_rho = Spatial_Clusters(iRun).rest_corr(iiiCluster,iiiCluster-iCluster+iiCluster);
                
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
title_label = 'Anatomical-Same-AAL-Rest-v2';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

%%% SAME ALL - PASSIVE

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI_clusters(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        for iiCluster=1:nClusters
            
            iPair = iPair + 1;
        
            for iRun=1:nRuns

                single_rho = Spatial_Clusters(iRun).passive_corr(iiiCluster,iiiCluster-iCluster+iiCluster);
                
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
title_label = 'Anatomical-Same-AAL-Passive-v2';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

%%% SAME ALL - TRACK

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI_clusters(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        for iiCluster=1:nClusters
            
            iPair = iPair + 1;
        
            for iRun=1:nRuns

                single_rho = Spatial_Clusters(iRun).track_corr(iiiCluster,iiiCluster-iCluster+iiCluster);
                
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
title_label = 'Anatomical-Same-AAL-Track-v2';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

%%% DIFFERENT AAL - RESTING STATE

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI_clusters(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        ivCluster = 0;
        for iiROI=1:nROI
            
            nnClusters = length(ROI_clusters(iiROI).clusters);
            
            for iiCluster=1:nnClusters
                
                ivCluster = ivCluster + 1;
                
                if iROI ~= iiROI
                    
                    iPair = iPair + 1;
                    
                    for iRun=1:nRuns
                        
                        single_rho = Spatial_Clusters(iRun).rest_corr(iiiCluster,ivCluster);
                        
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
title_label = 'Anatomical-Different-AAL-Rest-v2';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

%%% DIFFERENT AAL - PASSIVE

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI_clusters(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        ivCluster = 0;
        for iiROI=1:nROI
            
            nnClusters = length(ROI_clusters(iiROI).clusters);
            
            for iiCluster=1:nnClusters
                
                ivCluster = ivCluster + 1;
                
                if iROI ~= iiROI
                    
                    iPair = iPair + 1;
                    
                    for iRun=1:nRuns
                        
                        single_rho = Spatial_Clusters(iRun).passive_corr(iiiCluster,ivCluster);
                        
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
title_label = 'Anatomical-Different-AAL-Passive-v2';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

%%% DIFFERENT AAL - TRACK

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI_clusters(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        ivCluster = 0;
        for iiROI=1:nROI
            
            nnClusters = length(ROI_clusters(iiROI).clusters);
            
            for iiCluster=1:nnClusters
                
                ivCluster = ivCluster + 1;
                
                if iROI ~= iiROI
                    
                    iPair = iPair + 1;
                    
                    for iRun=1:nRuns
                        
                        single_rho = Spatial_Clusters(iRun).track_corr(iiiCluster,ivCluster);
                        
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
title_label = 'Anatomical-Different-AAL-Track-v2';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);


end

function totalInformationContent

load('Anatomical-Spatial-Clusters-v2-Corr.mat');

nTotalRuns = 32;
r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

for iRun=1:nTotalRuns

	all_corr(iRun,:,:) = Spatial_Clusters(iRun).rest_corr;
        
end

all_corr_fisher = (1/2).*log((1 + all_corr)./(1 - all_corr));
             
m_corr_z = squeeze(mean(all_corr_fisher,1));
s_corr_z = squeeze(std(all_corr_fisher,0,1));      

TotalInformationContent(m_corr_z(:), s_corr_z(:), z_threshold, cv_threshold );

end


