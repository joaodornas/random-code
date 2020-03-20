function real_methods_paper_v1_c400

getCorrFromC400Parcellation;

end

%%% CC400 PARCELLATION

function getCorrFromC400Parcellation

nClusters = 400;
nTR = 150;
nRuns = 4;
nTotalRuns = 32;
MNI_size = [91 109 91];
MNI_size_mm = [47 56 46];
mm = 4;

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
        % [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL_mm(settings,kind,irun,file,mask,get_at_this_preprocessed_step,mm);
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL_MAC(settings,kind,irun,file,mask,get_at_this_preprocessed_step,mm);
    
    end
    
end

%%% C400
% file = 'C400_4mm_tcorr_2level.nii';
file = 'tcorr05_2level_all_400.nii';
% img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\craddock_2011_parcellations\',file));
img = nifti(strcat('/Volumes/dropbox/_TOOLBOX/craddock_2011_parcellations/',file));
% img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\craddock_2011_parcellations\',file);
img.dat.fname = strcat('/Volumes/dropbox/_TOOLBOX/craddock_2011_parcellations/',file);

C400 = img.dat(:,:,:);

for iiRun=1:nTotalRuns
    
    for iCluster=1:nClusters
    
        idx_voxels = find(C400==iCluster); 
        nVoxels = length(idx_voxels);
   
        this_cluster_voxels = zeros(nVoxels,nTR);
        
        for iVoxel=1:nVoxels
        
            % [idxx,idxy,idxz] = ind2sub(MNI_size_mm,idx_voxels(iVoxel));
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

% save('C400-tcorr-2level-Corr.mat','all_corr');
save('tcorr05_2level_all_400-Corr.mat','all_corr');

end

function getCorrFromC400ParcellationInsideParcelVoxelLevel

nClusters = 400;
nTR = 150;
nRuns = 4;
nTotalRuns = 32;
MNI_size = [91 109 91];
MNI_size_mm = [47 56 46];
mm = 4;

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
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL_mm(settings,kind,irun,file,mask,get_at_this_preprocessed_step,mm);
    end
    
end

%%% C400
file = 'C400_4mm_tcorr_2level.nii';
img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\craddock_2011_parcellations\',file));
img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\craddock_2011_parcellations\',file);
C400 = img.dat(:,:,:);

for iiRun=1:nTotalRuns
    
    disp(int2str(iiRun));
    
    all_corr(iiRun).corr = [];
    
    for iCluster=1:nClusters
       
        idx_voxels = find(C400==iCluster); 
        nVoxels = length(idx_voxels);
   
        this_cluster_voxels = zeros(nTR,nVoxels);
        
        for iVoxel=1:nVoxels
        
            [idxx,idxy,idxz] = ind2sub(MNI_size_mm,idx_voxels(iVoxel));
            
            this_cluster_voxels(:,iVoxel) = RestingState(iiRun).run(idxx,idxy,idxz,:);
            
        end
        
        if ~isempty(idx_voxels); 
            
            corr_this_cluster_voxels = corr(this_cluster_voxels);
        
            all_corr(iiRun).corr = [all_corr(iiRun).corr, corr_this_cluster_voxels(:)'];
   
        end
        
    end

end

save('C400-4mm-2level-Inside-Corr.mat','all_corr','-v7.3');

end

function getCorrFromC400ParcellationInsideParcelVoxelLevelsavepercluster

nClusters = 400;
nTR = 150;
nRuns = 4;
nTotalRuns = 32;
MNI_size = [91 109 91];
MNI_size_4mm = [47 56 46];
mm = 4;

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
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL_mm(settings,kind,irun,file,mask,get_at_this_preprocessed_step,mm);
    end
    
end

%%% C400
file = 'C400_4mm_tcorr_2level.nii';
img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\craddock_2011_parcellations\',file));
img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\craddock_2011_parcellations\',file);
C400 = img.dat(:,:,:);
    
for iCluster=1:nClusters

    idx_voxels = find(C400==iCluster); 
    nVoxels = length(idx_voxels);

    clusters_corr(iCluster).corr = zeros(nTotalRuns,nVoxels,nVoxels);
    
    if ~isempty(idx_voxels)

        all_idx_zeros = [];
        for iiRun=1:nTotalRuns

            this_cluster_voxels = zeros(nTR,nVoxels);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_size_4mm,idx_voxels(iVoxel));

                this_cluster_voxels(:,iVoxel) = RestingState(iiRun).run(idxx,idxy,idxz,:);

            end

            s_voxels = sum(this_cluster_voxels,1);
            idx_zeros = find(s_voxels==0);
            all_idx_zeros = [all_idx_zeros,idx_zeros];

            clusters_corr(iCluster).corr(iiRun,:,:) = corr(this_cluster_voxels);

        end

        if ~isempty(all_idx_zeros)

            clusters_corr(iCluster).corr(:,unique(all_idx_zeros),:) = [];
            clusters_corr(iCluster).corr(:,:,unique(all_idx_zeros)) = [];

        end

        if iCluster < 5; plotFCVoxelsInsideAClusterC400(clusters_corr(iCluster).corr,iCluster); end
    
    end
    
end

save('C400-4mm-tcorr-2level-Inside-Corr.mat','clusters_corr','-v7.3');

end

function getAndPlotVarianceDistributionFromC400Parcellation

load('C400-4mm-tcorr-2level-Corr.mat');

nClusters = 400;
nTotalRuns = 32;

all_corr_z = (1/2).*log((1 + all_corr)./(1 - all_corr));
            
my_std = squeeze(std(all_corr_z,0,1));
my_mean = squeeze(mean(all_corr_z,1));

vector_one = my_std;
vector_two = my_mean;
label_one = 'STD';
label_two = 'Mean';
title_label = 'C400-tcorr-2level';

MyContourCorrelation_v5(vector_two(:),vector_one(:),title_label,1);

end

function getAndPlotVarianceDistributionFromC400ParcellationInsideParcelVoxelLevel

load('C400-4mm-tcorr-2level-Inside-Corr.mat');

nClusters = 400;
nTotalRuns = 32;

for iRun=1:nTotalRuns
   
    iPair = 0;
    for iCluster=1:nClusters
        
        nVoxels = size(clusters_corr(iCluster).corr,2);
    
        for iVoxel=1:nVoxels
            
            for iiVoxel=1:nVoxels
                
                iPair = iPair + 1;
                
                all_pairs(iRun,iPair) = clusters_corr(iCluster).corr(iRun,iVoxel,iiVoxel);
                
            end
            
        end
        
    end
    
end

all_corr_z = (1/2).*log((1 + all_pairs)./(1 - all_pairs));
            
my_std = squeeze(std(all_corr_z,0,1));
my_mean = squeeze(mean(all_corr_z,1));

vector_one = my_std;
vector_two = my_mean;
label_one = 'STD';
label_two = 'Mean';
title_label = 'C400-4mm-tcorr-2level-Inside';

MyContourCorrelation_v5(vector_two(:),vector_one(:),title_label,1);

end

function plotFCVoxelsInsideAClusterC400(voxels_corr,idx_cluster)

nRuns = 32;
nTR = 150;
nVoxels = size(voxels_corr,2);
nTotalRuns = 32;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;
for iVoxel=1:nVoxels
    
    for iiVoxel=1:nVoxels
    
        if iVoxel ~= iiVoxel
            
            iPair = iPair + 1;
            for iRun=1:nTotalRuns

                single_rho = squeeze(voxels_corr(iRun,iVoxel,iiVoxel));

                pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

                pairs_rest_rho(iRun,iPair) = single_rho;

            end
        
        end
    
    end

end
            
m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ mean(pairs_rest_z,1);


mean_corr = zeros(nVoxels,nVoxels);
s_corr = zeros(nVoxels,nVoxels);

iPair = 0;
for iVoxel=1:nVoxels
    
    for iiVoxel=1:nVoxels
        
        if iVoxel ~= iiVoxel
            
            iPair = iPair + 1;
        
            %if my_cv_aal(iPair) < cv_criterion & my_cv_aal(iPair) > 0 & p_pairs(iPair) < p_criterion
                
                mean_corr(iVoxel,iiVoxel) =  m_pairs_z(iPair);
                s_corr(iVoxel,iiVoxel) =  s_pairs_z(iPair);
                
            %end
            
        end
        
    end
    
end

new_label = strcat('FC-C400-Voxels-Cluster-',int2str(idx_cluster),'-','Voxels','-',int2str(nVoxels));
plotFCgeneral(mean_corr,s_corr,new_label);

close all;

end

function plotMeanFCC400

load('C400-4mm-tcorr-2level-Corr.mat');

all_corr_z = (1/2).*log((1 + all_corr)./(1 - all_corr));
            
my_std = squeeze(std(all_corr_z,0,1));
my_mean = squeeze(mean(all_corr_z,1));


new_label = strcat('FC-C400');
plotFCgeneral(my_mean,my_std,new_label);

end

function getInformationC400

load('C400-4mm-tcorr-2level-Corr.mat');

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

all_corr_z = (1/2).*log((1 + all_corr)./(1 - all_corr));
            
my_std = squeeze(nanstd(all_corr_z,0,1));
my_mean = squeeze(nanmean(all_corr_z,1));

new_label = strcat('Information-C400');
biInfo = TotalInformationContent( my_mean, my_std , z_threshold, cv_threshold );

save('C400-biInfo.mat','biInfo');

end



