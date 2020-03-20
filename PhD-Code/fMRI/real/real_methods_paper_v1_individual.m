function real_methods_paper_v1_individual

% getCorrVoxelsInsideALLClusterSpatial;

% plotVarianceDistributionClustersVoxelsInsideALLClusterSpatial;

% parcellation_label = 'SP758';
% ROI = getSPATIALClusters;
% getCorrFromFromVoxelsFromASampleOf20VoxelsPerClusterV2(ROI,parcellation_label);
% 
% parcellation_label = 'MD758';
% ROI = getMD758Clusters;
% getCorrFromFromVoxelsFromASampleOf20VoxelsPerClusterV2(ROI,parcellation_label);

% condition_label = 'rest';
% parcellation_label = 'SP758';
% ROI = getSPATIALClusters;
% plotVarianceDistributionVoxelsFromASampleOf20VoxelsPerClusterV2(ROI,condition_label,parcellation_label);

condition_label = 'passive';
parcellation_label = 'SP758';
ROI = getSPATIALClusters;
plotVarianceDistributionVoxelsFromASampleOf20VoxelsPerClusterV2(ROI,condition_label,parcellation_label);

condition_label = 'track';
parcellation_label = 'SP758';
ROI = getSPATIALClusters;
plotVarianceDistributionVoxelsFromASampleOf20VoxelsPerClusterV2(ROI,condition_label,parcellation_label);

% condition_label = 'rest';
% parcellation_label = 'MD758';
% ROI = getMD758Clusters;
% plotVarianceDistributionVoxelsFromASampleOf20VoxelsPerClusterV2(ROI,condition_label,parcellation_label);

condition_label = 'passive';
parcellation_label = 'MD758';
ROI = getMD758Clusters;
plotVarianceDistributionVoxelsFromASampleOf20VoxelsPerClusterV2(ROI,condition_label,parcellation_label);

condition_label = 'track';
parcellation_label = 'MD758';
ROI = getMD758Clusters;
plotVarianceDistributionVoxelsFromASampleOf20VoxelsPerClusterV2(ROI,condition_label,parcellation_label);

end

%%% INDIVIDUAL VOXELS

function getCorrFromFromVoxelsFromASampleOf20VoxelsPerClusterV2(ROI,parcellation_label)


nClusters = 758;
nROI = 90;
nRuns = 4;
MNI_img = [91 109 91];
nTotalVoxels = 160990;
nTotalClusters = 758;
nClusterSize = 200;
nSamples = 20;
nTR = 150;

all_idx_voxels = [];

for iROI=1:nROI
   
    nClusters = length(ROI(iROI).clusters);
    
    for iCluster=1:nClusters
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
        nVoxels = length(idx_voxels);
        
        if nVoxels < nSamples; nSamples_effective = nVoxels; else nSamples_effective = nSamples; end
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels(1:nSamples_effective);
        
        all_idx_voxels = [all_idx_voxels;idx_voxels];
        
    end
    
end

all_settings = getAllSettings;

nSubjects = length(all_settings);

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\_ATLAS\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

% %%% REST
% 
% voxels_corr_rest = zeros(nRuns*nSubjects,length(all_idx_voxels),length(all_idx_voxels));
% 
% iiRun = 0;
% 
% for iSubject=1:nSubjects
%     
%     settings = all_settings(iSubject).settings;
%     
%     %% LOAD DATA
%     get_at_this_preprocessed_step = settings.FSL.folders.custom;
%     file = settings.FSL.files.functional.custom.residual_voxel;
%     mask = settings.FSL.files.mask.custom;
%     
%     kind = 'RestingState';
%     for irun=1:nRuns;
%         [RestingState(irun).run, RestingState(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
%     end
%     
%     for irun=1:nRuns
%         
%         iiRun = iiRun + 1;
%         
%         rest_sample_voxels = zeros(nTR,length(all_idx_voxels));
%         
%         disp(strcat('Run:',int2str(iiRun)));
%         
%         for iVoxel=1:length(all_idx_voxels)
% 
%             [idxx,idxy,idxz] = ind2sub(size(AAL_img),all_idx_voxels(iVoxel));
%             
%             rest_sample_voxels(:,iVoxel) = RestingState(irun).run(idxx,idxy,idxz,:);
%             
%         end
%         
%         voxels_corr_rest(iiRun,:,:) = corr(rest_sample_voxels);
%         
%     end
% 
% end
% 
% disp('I will save right now...');
% 
% save(strcat('Samples-Of-Voxels-Corr-Rest-',parcellation_label,'.mat'),'voxels_corr_rest','-v7.3');

%%% PASSIVE

voxels_corr_passive = zeros(nRuns*nSubjects,length(all_idx_voxels),length(all_idx_voxels));

iiRun = 0;

for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
        
    kind = 'Passive';
    for irun=1:nRuns;
        [Passive(irun).run, Passive(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
    for irun=1:nRuns
        
        iiRun = iiRun + 1;
        
        passive_sample_voxels = zeros(nTR,length(all_idx_voxels));
        
        disp(strcat('Run:',int2str(iiRun)));
        
        for iVoxel=1:length(all_idx_voxels)

            [idxx,idxy,idxz] = ind2sub(size(AAL_img),all_idx_voxels(iVoxel));
            
            passive_sample_voxels(:,iVoxel) = Passive(irun).run(idxx,idxy,idxz,:);
            
        end
        
        voxels_corr_passive(iiRun,:,:) = corr(passive_sample_voxels);
        
    end

end

disp('I will save right now...');

save(strcat('Samples-Of-Voxels-Corr-Passive-',parcellation_label,'.mat'),'voxels_corr_passive','-v7.3');

%%% TRACK

voxels_corr_track = zeros(nRuns*nSubjects,length(all_idx_voxels),length(all_idx_voxels));

iiRun = 0;

for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
        
    kind = 'Track';
    for irun=1:nRuns;
        [Track(irun).run, Track(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
    for irun=1:nRuns
        
        iiRun = iiRun + 1;
        
        track_sample_voxels = zeros(nTR,length(all_idx_voxels));
        
        disp(strcat('Run:',int2str(iiRun)));
        
        for iVoxel=1:length(all_idx_voxels)

            [idxx,idxy,idxz] = ind2sub(size(AAL_img),all_idx_voxels(iVoxel));
            
            track_sample_voxels(:,iVoxel) = Track(irun).run(idxx,idxy,idxz,:);
            
        end
        
        voxels_corr_track(iiRun,:,:) = corr(track_sample_voxels);
        
    end

end

disp('I will save right now...');

save(strcat('Samples-Of-Voxels-Corr-Track-',parcellation_label,'.mat'),'voxels_corr_track','-v7.3');

end

function plotVarianceDistributionVoxelsFromASampleOf20VoxelsPerClusterV2(ROI,condition_label,parcellation_label)

load(strcat('Samples-Of-Voxels-Corr','-',condition_label,'-',parcellation_label,'.mat'));

eval(strcat('voxels_corr = voxels_corr','_',lower(condition_label),';'));

nROI = 90;
nRuns = 32;
nSamplesOfVoxels = 20;

nTotalVoxels = size(voxels_corr,2);

iPair_same_cluster = 0;
iPair_same_aal = 0;
iPair_diff_aal = 0;
IVoxel = 0;
for iROI=1:nROI
    nClusters = length(ROI(iROI).clusters);
    for iCluster=1:nClusters
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
        nVoxels = length(idx_voxels);
        
        if nVoxels < nSamplesOfVoxels; i_nSamples_effective = nVoxels; else i_nSamples_effective = nSamplesOfVoxels; end
        
       for iVoxel=1:i_nSamples_effective
          IVoxel = IVoxel + 1;
          if mod(IVoxel,1000) == 0; disp(int2str(IVoxel)); end
          IIVoxel = 0;
          for iiROI=1:nROI
              nnClusters = length(ROI(iiROI).clusters);
              for iiCluster=1:nnClusters
                  
                idx_voxels = ROI(iiROI).clusters(iiCluster).idx_voxels;
                nnVoxels = length(idx_voxels);
        
                if nnVoxels < nSamplesOfVoxels; ii_nSamples_effective = nnVoxels; else ii_nSamples_effective = nSamplesOfVoxels; end
        
                  for iiVoxel=1:ii_nSamples_effective
                      IIVoxel = IIVoxel + 1;
                      if IIVoxel >= IVoxel
                          if iROI == iiROI & iCluster == iiCluster 
                              iPair_same_cluster = iPair_same_cluster + 1;
                          end
                          if iROI == iiROI & iCluster ~= iiCluster
                              iPair_same_aal = iPair_same_aal + 1;
                          end
                          if iROI ~= iiROI
                              iPair_diff_aal = iPair_diff_aal + 1;
                          end
                      end
                  end
              end
          end
       end
    end
end

pairs_same_cluster = zeros(nRuns,iPair_same_cluster);
pairs_same_aal = zeros(nRuns,iPair_same_aal);
pairs_diff_aal = zeros(nRuns,iPair_diff_aal);

iPair_same_cluster = 0;
iPair_same_aal = 0;
iPair_diff_aal = 0;
IVoxel = 0;
for iROI=1:nROI
    
    nClusters = length(ROI(iROI).clusters);
   
    for iCluster=1:nClusters
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
        nVoxels = length(idx_voxels);
        
        if nVoxels < nSamplesOfVoxels; i_nSamples_effective = nVoxels; else i_nSamples_effective = nSamplesOfVoxels; end
       
       for iVoxel=1:i_nSamples_effective
           
          IVoxel = IVoxel + 1;
          
          if mod(IVoxel,1000) == 0; disp(int2str(IVoxel)); end
          
          IIVoxel = 0;
          
          for iiROI=1:nROI
              
              nnClusters = length(ROI(iiROI).clusters);
              
              for iiCluster=1:nnClusters
                  
                idx_voxels = ROI(iiROI).clusters(iiCluster).idx_voxels;
                nnVoxels = length(idx_voxels);
        
                if nnVoxels < nSamplesOfVoxels; ii_nSamples_effective = nnVoxels; else ii_nSamples_effective = nSamplesOfVoxels; end
        
                  for iiVoxel=1:ii_nSamples_effective
                      
                      IIVoxel = IIVoxel + 1;
                      
                      if IIVoxel >= IVoxel
                    
                          if iROI == iiROI & iCluster == iiCluster 

                              iPair_same_cluster = iPair_same_cluster + 1;

                               single_rho = squeeze(voxels_corr(:,IVoxel,IIVoxel));
            
                               pairs_same_cluster(:,iPair_same_cluster) = (1/2)*log((1 + single_rho)./(1 - single_rho));

                          end

                          if iROI == iiROI & iCluster ~= iiCluster

                              iPair_same_aal = iPair_same_aal + 1;

                              single_rho = squeeze(voxels_corr(:,IVoxel,IIVoxel));
                              
                              pairs_same_aal(:,iPair_same_aal) = (1/2)*log((1 + single_rho)./(1 - single_rho));

                          end

                          if iROI ~= iiROI

                              iPair_diff_aal = iPair_diff_aal + 1;

                              single_rho = squeeze(voxels_corr(:,IVoxel,IIVoxel));
                              
                              pairs_diff_aal(:,iPair_diff_aal) = (1/2)*log((1 + single_rho)./(1 - single_rho));

                          end
                      
                      end

                  end
                  
              end
              
          end
           
       end
        
    end
    
end

cv_threshold = 0.5;
% r_threshold = 0.186728;  % corresponds to p-value = 0.01 for N=150
r_threshold = 0.3;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );

my_var_same_cluster = nanvar(pairs_same_cluster);

my_cv_same_cluster = nanstd(pairs_same_cluster) ./ abs(nanmean(pairs_same_cluster,1));

my_mean_same_cluster = nanmean(pairs_same_cluster,1);
my_std_same_cluster = nanstd(pairs_same_cluster);

s = sum(pairs_same_cluster,1);
idx_zeros = find(s==0);

my_var_same_cluster(idx_zeros) = [];
my_cv_same_cluster(idx_zeros) = [];
my_mean_same_cluster(idx_zeros) = [];
my_std_same_cluster(idx_zeros) = [];

vector_one = my_std_same_cluster;
vector_two = my_mean_same_cluster;
vector_three = my_cv_same_cluster;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = strcat('Voxels-Same-Cluster','-',condition_label,'-',parcellation_label);

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);
fraction = FractionSignificantConsistent( vector_two, vector_one, z_threshold, cv_threshold );

my_var_same_aal = var(pairs_same_aal);

my_cv_same_aal = std(pairs_same_aal) ./ abs(mean(pairs_same_aal,1));

my_mean_same_aal = mean(pairs_same_aal,1);
my_std_same_aal = std(pairs_same_aal);

s = sum(pairs_same_aal,1);
idx_zeros = find(s==0);

my_var_same_aal(idx_zeros) = [];
my_cv_same_aal(idx_zeros) = [];
my_mean_same_aal(idx_zeros) = [];
my_std_same_aal(idx_zeros) = [];

vector_one = my_std_same_aal;
vector_two = my_mean_same_aal;
vector_three = my_cv_same_aal;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = strcat('Voxels-Same-AAL','-',condition_label,'-',parcellation_label);

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

my_var_diff_aal = var(pairs_diff_aal);

my_cv_diff_aal = std(pairs_diff_aal) ./ abs(mean(pairs_diff_aal,1));

my_mean_diff_aal = mean(pairs_diff_aal,1);
my_std_diff_aal = std(pairs_diff_aal);

s = sum(pairs_diff_aal,1);
idx_zeros = find(s==0);

my_var_diff_aal(idx_zeros) = [];
my_cv_diff_aal(idx_zeros) = [];
my_mean_diff_aal(idx_zeros) = [];
my_std_diff_aal(idx_zeros) = [];

vector_one = my_std_diff_aal;
vector_two = my_mean_diff_aal;
vector_three = my_cv_diff_aal;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = strcat('Voxels-Diff-AAL','-',condition_label,'-',parcellation_label);

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function getCorrVoxelsInsideACluster(idx_cluster)

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

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
        
    kind = 'RestingState';
    for irun=1:nRuns;
        iiRun = iiRun + 1;
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
end

for iROI=1:nROI
    
    nClusters = length(ROI(iROI).clusters);
    
    for iCluster=1:nClusters
        
        iiClusters = iiClusters + 1;
        
        if iiClusters == idx_cluster
            
            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
            
            nVoxels = length(idx_voxels);
            
            for iRun=1:nTotalRuns
            
                voxels_ts = zeros(nTR,nVoxels);
            
                for iVoxel=1:nVoxels
               
                    [idxx,idxy,idxz] = ind2sub(MNI_img,idx_voxels(iVoxel));
                
                    voxels_ts(:,iVoxel) = squeeze(RestingState(iRun).run(idxx,idxy,idxz,:));
                    
                end
                
                voxels_corr(iRun,:,:) = corr(voxels_ts);
                
            end
            
        end
        
    end

end

save(strcat('Voxels-Inside-Cluster-Corr-',int2str(idx_cluster),'.mat'),'voxels_corr');

end

function getCorrVoxelsInsideAClusterPetra

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nROI = 90;
iiClusters = 0;
nTR = 150;
MNI_img = [91 109 91];
nTotalRuns = 9;
nSegments = 3;

all_settings = getAllSettingsPetra;
nSubjects = length(all_settings);

for iSubject=1:nSubjects
    
    run = 1;
    
    settings = all_settings(iSubject).settings;
    
    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;

    [RestingState(iSubject).run, settings] = real_get_data_FSL_Petra(settings,run,file,get_at_this_preprocessed_step);
    
end

iiClusters = 0;
for iROI=1:nROI
    
    nClusters = length(ROI(iROI).clusters);
    
    for iCluster=1:nClusters
        
        iiClusters = iiClusters + 1;
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

        nVoxels = length(idx_voxels);

        iiRunSeg = 0;
        for iRun=1:nTotalRuns

            voxels_ts = zeros(nTR,nVoxels);

            for iSeg=1:nSegments
                
                iiRunSeg = iiRunSeg + 1;
                
                for iVoxel=1:nVoxels
                    
                    [idxx,idxy,idxz] = ind2sub(MNI_img,idx_voxels(iVoxel));

                    voxels_ts(:,iVoxel) = squeeze(RestingState(iRun).run(idxx,idxy,idxz,(1+(iSeg-1)*nTR):(nTR+(iSeg-1)*nTR)));

                end
                
                voxels_corr(iiRunSeg,:,:) = corr(voxels_ts);

            end
  

        end

        save(strcat('Voxels-Inside-Cluster-Petra-Corr-',int2str(iiClusters),'.mat'),'voxels_corr');
            
        clear voxels_corr
        
   end
        
end

end

function getCorrVoxelsInsideALLCluster

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat');

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
        
    kind = 'RestingState';
    for irun=1:nRuns;
        iiRun = iiRun + 1;
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
end

for iROI=1:nROI
    
    nClusters = length(ROI(iROI).clusters);
    
    for iCluster=1:nClusters
        
        iiClusters = iiClusters + 1;
        
        iiClusters
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

        nVoxels = length(idx_voxels);

        for iRun=1:nTotalRuns

            voxels_ts = zeros(nTR,nVoxels);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_img,idx_voxels(iVoxel));

                voxels_ts(:,iVoxel) = squeeze(RestingState(iRun).run(idxx,idxy,idxz,:));

            end

            voxels_corr(iRun,:,:) = corr(voxels_ts);

        end
        
        save(strcat('Voxels-Inside-Cluster-Corr-',int2str(iiClusters),'.mat'),'voxels_corr');
        
        clear voxels_corr

    end

end



end

function getCorrVoxelsInsideALLClusterSpatial

load('Anatomical-Spatial-Clusters-v2.mat');

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
        
    kind = 'RestingState';
    for irun=1:nRuns;
        iiRun = iiRun + 1;
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
end

for iROI=1:nROI
    
    nClusters = length(ROI_clusters(iROI).clusters);
    
    for iCluster=1:nClusters
        
        iiClusters = iiClusters + 1;
        
        iiClusters
        
        idx_voxels = ROI_clusters(iROI).clusters(iCluster).idx_voxels;

        nVoxels = length(idx_voxels);

        for iRun=1:nTotalRuns

            voxels_ts = zeros(nTR,nVoxels);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_img,idx_voxels(iVoxel));

                voxels_ts(:,iVoxel) = squeeze(RestingState(iRun).run(idxx,idxy,idxz,:));

            end

            voxels_corr(iRun,:,:) = corr(voxels_ts);

        end
        
        save(strcat('Voxels-Inside-Cluster-Corr-Spatial-v2-',int2str(iiClusters),'.mat'),'voxels_corr');
        
        clear voxels_corr

    end

end



end

function plotVarianceDistributionClustersVoxelsInsideALLCluster

nTotalClusters = 758;
nRuns = 32;

pairs = [];
for iCluster=1:nTotalClusters
    
    load(strcat('Voxels-Inside-Cluster-Corr-',int2str(iCluster)));
    nVoxels = size(voxels_corr,2);
    v = reshape(voxels_corr,[nRuns nVoxels*nVoxels]);
    clear voxels_corr
    pairs = [pairs,v];

end

pairs_z = (1/2).*log((1 + pairs)./(1 - pairs));

%%% EVERYBODY - VOXELS INSIDE - REST

my_var = var(pairs_z);

my_cv = std(pairs_z) ./ abs(mean(pairs_z,1));

my_mean = mean(pairs_z,1);
my_std = std(pairs_z);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Everybody-Voxels-Inside-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotVarianceDistributionClustersVoxelsInsideALLClusterSpatial

nTotalClusters = 758;
nRuns = 32;

pairs = [];
for iCluster=1:nTotalClusters
    
    load(strcat('Voxels-Inside-Cluster-Corr-Spatial-v2-',int2str(iCluster)));
    nVoxels = size(voxels_corr,2);
    v = reshape(voxels_corr,[nRuns nVoxels*nVoxels]);
    clear voxels_corr
    pairs = [pairs,v];

end

pairs_z = (1/2).*log((1 + pairs)./(1 - pairs));

%%% EVERYBODY - VOXELS INSIDE - REST

my_var = var(pairs_z);

my_cv = std(pairs_z) ./ abs(mean(pairs_z,1));

my_mean = mean(pairs_z,1);
my_std = std(pairs_z);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Spatial-Everybody-Voxels-Inside-Rest-v2';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotVarianceDistributionClustersVoxelsInsideALLClusterPetra

nTotalClusters = 758;
nRuns = 27;

pairs = [];
for iCluster=1:nTotalClusters
    
    load(strcat('Voxels-Inside-Cluster-Petra-Corr-',int2str(iCluster)));
    nVoxels = size(voxels_corr,2);
    v = reshape(voxels_corr,[nRuns nVoxels*nVoxels]);
    clear voxels_corr
    pairs = [pairs,v];

end

pairs_z = (1/2).*log((1 + pairs)./(1 - pairs));

%%% EVERYBODY - VOXELS INSIDE - REST

my_var = var(pairs_z);

my_cv = std(pairs_z) ./ abs(mean(pairs_z,1));

my_mean = mean(pairs_z,1);
my_std = std(pairs_z);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Everybody-Voxels-Inside-Petra-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotFCVoxelsInsideACluster(idx_cluster)

load(strcat('Voxels-Inside-Cluster-Corr-',int2str(idx_cluster),'.mat'));

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

new_label = strcat('FC-Voxels-Cluster-',int2str(idx_cluster),'-','Voxels','-',int2str(nVoxels));
plotFCgeneral(mean_corr,s_corr,new_label);

close all;

end

function getInformationRunAALVoxelLevel


nROIs = 90;
nSubjects = 8;
nRuns = 4;
nTotalRuns = 32;

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

all_settings = getAllSettings;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

biInfo = zeros(nSubjects,nROIs);
for iROI=1:nROIs
   
    area_Label = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    
    disp(strcat(area_Label,'-',datestr(now)));
    
    iiRun = 0;
    for iSubject=1:nSubjects
       
        settings = all_settings(iSubject).settings;
        
        load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-RestingState-',area_Label,'.mat'));
        
        for iRun=1:nRuns
            
            iiRun = iiRun + 1;
            
            single_rho = FC_Voxels.run(iRun).rho_RestingState;
            
            fisher(iRun,:,:) = (1/2).*log((1 + single_rho)./(1 - single_rho));
            
        end
        
        [l,c,z] = size(fisher);
        
        fisher_pairs = reshape(fisher,[l c*z]);
        
        mu_z = mean(fisher_pairs,1);
        sigma_z = std(fisher_pairs);
   
        biInfo(iSubject,iROI) = TotalInformationContent( mu_z, sigma_z, z_threshold, cv_threshold );
            
        clear FC_Voxels subject file area_label fisher fisher_pairs mu_z sigma_z
        
    end
    
end

save('LHR-biInfo.mat','biInfo');

end

function getInformationRunAALVoxelLevelPETRA


nROIs = 90;
nSubjects = 9;
nRuns = 3;
nTotalRuns = 27;

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

all_settings = getAllSettingsPetra;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

biInfo = zeros(nSubjects,nROIs);
for iROI=1:nROIs
   
    area_Label = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    
    disp(strcat(area_Label,'-',datestr(now)));
    
    iiRun = 0;
    for iSubject=1:nSubjects
       
        settings = all_settings(iSubject).settings;
        
        load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-Petra-RestingState-',area_Label,'.mat'));
        
        for iRun=1:nRuns
            
            iiRun = iiRun + 1;
            
            single_rho = FC_Voxels.run(iRun).rho_RestingState;
            
            fisher(iRun,:,:) = (1/2).*log((1 + single_rho)./(1 - single_rho));
            
        end
        
        [l,c,z] = size(fisher);
        
        fisher_pairs = reshape(fisher,[l c*z]);
        
        mu_z = mean(fisher_pairs,1);
        sigma_z = std(fisher_pairs);
   
        biInfo(iSubject,iROI) = TotalInformationContent( mu_z, sigma_z, z_threshold, cv_threshold );
            
        clear FC_Voxels subject file area_label fisher fisher_pairs mu_z sigma_z
        
    end
    
end

save('PETRA-biInfo.mat','biInfo');

end

function getInformationRunAALVoxelLevelHCP


nROIs = 90;
nSubjects = 8;
nRuns = 4;
nTotalRuns = 32;

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

biInfo = zeros(nSubjects,nROIs);
for iROI=1:nROIs
   
    area_Label = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    
    disp(strcat(area_Label,'-',datestr(now)));
    
    iiRun = 0;
    for iSubject=1:nSubjects

        for iRun=1:nRuns
            
            load(strcat('HCP','-','SUBJ',int2str(iSubject),'-','FC-Voxels-AAL-ROI-corr-HCP-RestingState-',area_Label,'-',int2str(iRun),'.mat'));
            
            iiRun = iiRun + 1;
            
            single_rho = FC_Voxels.run(iRun).rho_RestingState;
            
            fisher(iRun,:,:) = (1/2).*log((1 + single_rho)./(1 - single_rho));
            
        end
        
        [l,c,z] = size(fisher);
        
        fisher_pairs = reshape(fisher,[l c*z]);
        
        mu_z = mean(fisher_pairs,1);
        sigma_z = std(fisher_pairs);
   
        biInfo(iSubject,iROI) = TotalInformationContent( mu_z, sigma_z, z_threshold, cv_threshold );
            
        clear FC_Voxels subject file area_label fisher fisher_pairs mu_z sigma_z
        
    end
    
end

save('HCP-biInfo.mat','biInfo');

end

function plotInfoLHRPETRA

nROI = 90;

LHR = load('LHR-biInfo.mat');
PETRA = load('PETRA-biInfo.mat');

biInfo = [LHR.biInfo;PETRA.biInfo];

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for iROI=1:nROI
    
   idx_ROI  = AAL_ROI(iROI).ID; 
   idx_voxels = find(AAL_img == idx_ROI);
   
   nVoxels = length(idx_voxels);
   
   %biInfo(:,iROI) = biInfo(:,iROI) ./ ( nVoxels * ( nVoxels - 1 ) / 2 );
   
   LHR.biInfo(:,iROI) = LHR.biInfo(:,iROI) ./ ( nVoxels * ( nVoxels - 1 ) / 2 );
   PETRA.biInfo(:,iROI) = PETRA.biInfo(:,iROI) ./ ( nVoxels * ( nVoxels - 1 ) / 2 );
   
   
end

[H P] = ttest2(LHR.biInfo,PETRA.biInfo);

m_LHR = mean(LHR.biInfo,1);
m_PETRA = mean(PETRA.biInfo,1);

most_LHR = (m_LHR > m_PETRA) & H;
most_PETRA = (m_PETRA > m_LHR) & H;

for iROI=1:nROI
    
    if most_LHR(iROI)
        
        disp(strrep(AAL_ROI(iROI).Nom_L,'_','-'));
        
    end
    
end

imagesc(biInfo);

hold on
plot([0 90],[8.5 8.5],'w-');

for i=1:90
    my_labels{i} = AAL_ROI(i).Nom_L;
end

set(gca,'XTick',1:90);
set(gca,'XTickLabel',my_labels);

xticklabel_rotate;

colorbar;

end

function plotInfoLHRHCP

nROI = 90;

LHR = load('LHR-biInfo.mat');
HCP = load('HCP-biInfo.mat');

biInfo = [LHR.biInfo;HCP.biInfo];

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for iROI=1:nROI
    
   idx_ROI  = AAL_ROI(iROI).ID; 
   idx_voxels = find(AAL_img == idx_ROI);
   
   nVoxels = length(idx_voxels);
   
   %biInfo(:,iROI) = biInfo(:,iROI) ./ ( nVoxels * ( nVoxels - 1 ) / 2 );
   
   LHR.biInfo(:,iROI) = LHR.biInfo(:,iROI) ./ ( nVoxels * ( nVoxels - 1 ) / 2 );
   HCP.biInfo(:,iROI) = HCP.biInfo(:,iROI) ./ ( nVoxels * ( nVoxels - 1 ) / 2 );
   
   
end

[H P] = ttest2(LHR.biInfo,HCP.biInfo);

m_LHR = mean(LHR.biInfo,1);
m_HCP = mean(HCP.biInfo,1);

most_LHR = (m_LHR > m_HCP) & H;
most_HCP = (m_HCP > m_LHR) & H;

for iROI=1:nROI
    
    if most_LHR(iROI)
        
        disp(strrep(AAL_ROI(iROI).Nom_L,'_','-'));
        
    end
    
end

imagesc(biInfo);

hold on
plot([0 90],[8.5 8.5],'w-');

for i=1:90
    my_labels{i} = AAL_ROI(i).Nom_L;
end

set(gca,'XTick',1:90);
set(gca,'XTickLabel',my_labels);

xticklabel_rotate;

colorbar;

end

function totalInformationContentVoxelsInsideSpatial

nTotalClusters = 758;
nRuns = 32;

pairs = [];
for iCluster=1:nTotalClusters
    
    load(strcat('Voxels-Inside-Cluster-Corr-Spatial-v2-',int2str(iCluster)));
    nVoxels = size(voxels_corr,2);
    v = reshape(voxels_corr,[nRuns nVoxels*nVoxels]);
    clear voxels_corr
    pairs = [pairs,v];

end

pairs_z = (1/2).*log((1 + pairs)./(1 - pairs));

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;
             
m_corr_z = squeeze(mean(pairs_z,1));
s_corr_z = squeeze(std(pairs_z,0,1));      

TotalInformationContent(m_corr_z(:), s_corr_z(:), z_threshold, cv_threshold );

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

function ROI = getMD758Clusters

load('Z:\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat');

end

