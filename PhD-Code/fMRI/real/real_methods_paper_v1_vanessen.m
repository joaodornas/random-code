function real_methods_paper_v1_vanessen

%%% VAN ESSEN PARCELLATION

function convert180to360VanEssen

nClusters = 180;
MNI_size = [91 109 91];

file = 'HCP-MMP1_on_MNI152_ICBM2009a_nlin-2mm.nii';

img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file));
img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file);
vanEssen = img.dat(:,:,:);

for iCluster=1:nClusters
    
    idx_voxels = find(vanEssen==iCluster);
    
    nVoxels = length(idx_voxels);
    
    for iVoxel=1:nVoxels
        
        [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
        
        if idxx > MNI_size(1)/2
            
            vanEssen(idxx,idxy,idxz) = vanEssen(idxx,idxy,idxz) + nClusters;
            
        end
        
    end
    
end

nifti_file = img;
offset = img.dat.offset;
scl_slope = img.dat.scl_slope;
scl_inter = img.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = img.dat.dim;

descrip = 'VanEssen-LR';

fname = strcat('VanEssen-LR-2mm','.nii');
input_data = vanEssen; 
real_save_image;

end

function getCorrFromVanEssenParcellation

%nClusters = 180;
nClusters = 360;
nTR = 150;
nRuns = 4;
nTotalRuns = 32;
MNI_size = [91 109 91];

all_settings = getAllSettings;
nSubjects = length(all_settings);

mean_area_RestingState_voxel = zeros(nClusters,nTR);

iiRun = 0;
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

%%% 180
% file = 'HCP-MMP1_on_MNI152_ICBM2009a_nlin-2mm.nii';
% img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file));
% img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file);
% vanEssen = img.dat(:,:,:);

%%% 360
file = 'VanEssen-LR-2mm.nii';
img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file));
img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file);
vanEssen = img.dat(:,:,:);

for iiRun=1:nTotalRuns
    
    for iCluster=1:nClusters
    
        idx_voxels = find(vanEssen==iCluster); 
        nVoxels = length(idx_voxels);
   
        this_cluster_voxels = zeros(nVoxels,nTR);
        
        for iVoxel=1:nVoxels
        
            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
            
            this_cluster_voxels(iVoxel,:) = RestingState(iiRun).run(idxx,idxy,idxz,:);
            
        end
        
        m_this_cluster_voxels = mean(this_cluster_voxels,1);
        
        all_clusters(iiRun).run(:,iCluster) = m_this_cluster_voxels(:);
   
    end

end

for iiRun=1:nTotalRuns
   
    all_corr(iiRun,:,:) = corr(all_clusters(iiRun).run);
    
end

save('VanEssen-Corr.mat','all_corr');

end

function getCorrFromVanEssenParcellationInsideParcelVoxelLevel

%nClusters = 180;
nClusters = 360;
nTR = 150;
nRuns = 4;
nTotalRuns = 32;
MNI_size = [91 109 91];

all_settings = getAllSettings;
nSubjects = length(all_settings);

mean_area_RestingState_voxel = zeros(nClusters,nTR);

iiRun = 0;
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

%%% 180
% file = 'HCP-MMP1_on_MNI152_ICBM2009a_nlin-2mm.nii';
% img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file));
% img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file);
% vanEssen = img.dat(:,:,:);

%%% 360
file = 'VanEssen-LR-2mm.nii';
img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file));
img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file);
vanEssen = img.dat(:,:,:);

for iiRun=1:nTotalRuns
    
    all_corr(iiRun).corr = [];
    
    for iCluster=1:nClusters
        
        disp(int2str(iCluster));
    
        idx_voxels = find(vanEssen==iCluster); 
        nVoxels = length(idx_voxels);
   
        this_cluster_voxels = zeros(nTR,nVoxels);
        
        for iVoxel=1:nVoxels
        
            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
            
            this_cluster_voxels(:,iVoxel) = RestingState(iiRun).run(idxx,idxy,idxz,:);
            
        end
        
        corr_this_cluster_voxels = corr(this_cluster_voxels);
        
        all_corr(iiRun).corr = [all_corr(iiRun).corr, corr_this_cluster_voxels(:)'];
   
    end

end

save('VanEssen-Inside-Corr.mat','all_corr','-v7.3');

end

function getAndPlotVarianceDistributionFromVanEssenParcellation

load('VanEssen-360-Corr.mat');

% nClusters = 180;
nClusters = 360;
nTotalRuns = 32;

all_corr_z = (1/2).*log((1 + all_corr)./(1 - all_corr));
            
my_std = squeeze(std(all_corr_z,0,1));
my_mean = squeeze(mean(all_corr_z,1));

vector_one = my_std;
vector_two = my_mean;
label_one = 'STD';
label_two = 'Mean';
title_label = 'VanEssen';

MyContourCorrelation_v5(vector_two(:),vector_one(:),title_label,1);

end

function getAndPlotVarianceDistributionFromVanEssenParcellationInsideParcelVoxelLevel

load('VanEssen-Inside-Corr.mat');

% nClusters = 180;
nClusters = 360;
nTotalRuns = 32;

all_corr_struct = all_corr;
clear all_corr;

all_corr = zeros(nTotalRuns,length(all_corr_struct(1).corr));

for iRun=1:nTotalRuns
    
    all_corr(iRun,:) = all_corr_struct(iRun).corr(:);
    
end

all_corr_z = (1/2).*log((1 + all_corr)./(1 - all_corr));
            
my_std = squeeze(std(all_corr_z,0,1));
my_mean = squeeze(mean(all_corr_z,1));

vector_one = my_std;
vector_two = my_mean;
label_one = 'STD';
label_two = 'Mean';
title_label = 'VanEssen-Inside';

MyContourCorrelation_v5(vector_two(:),vector_one(:),title_label,1);

end

function plotMeanFCVanEssen

load('VanEssen-360-Corr.mat');

% nClusters = 180;
nClusters = 360;
nTotalRuns = 32;

all_corr_z = (1/2).*log((1 + all_corr)./(1 - all_corr));
            
my_std = squeeze(std(all_corr_z,0,1));
my_mean = squeeze(mean(all_corr_z,1));


new_label = strcat('FC-VanEssen');
plotFCgeneral(my_mean,my_std,new_label);

end

function getInformationVanEssen

load('VanEssen-360-Corr.mat');

% nClusters = 180;
nClusters = 360;
nTotalRuns = 32;

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

all_corr_z = (1/2).*log((1 + all_corr)./(1 - all_corr));
            
my_std = squeeze(nanstd(all_corr_z,0,1));
my_mean = squeeze(nanmean(all_corr_z,1));

new_label = strcat('Information-VanEssen');
biInfo = TotalInformationContent( my_mean, my_std , z_threshold, cv_threshold );

save('VanEssen-biInfo.mat','biInfo');

end

end

