
%% NATIVE PRISMA SPACE
load_aal = nifti('AAL-Prisma.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_DATA\Parcellation\AAL\',load_aal.dat.fname);

%% MNI SPACE
% load_aal = nifti('ROI_MNI_V4.nii');
% load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_Temporal_Poles = [83 84 87 88];
nROI = length(idx_Temporal_Poles);

%% LOAD DATA

%% MOT 4 - RUN 1
%load_nii = nifti('20150508_100617BOLDfMRI2x2x303p3run1s003a001.nii');
%load_nii.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\Jan-08-05-2015\preprocessed\T2-Stimulus-MOT-4-Balls\Run-1-1\FSL\4D-original\',load_nii.dat.fname);

%% MOT 4 - RUN 1 - MCFLIRT
%load_nii = nifti('MOT4-Run1_mcf.nii');
%load_nii.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\Jan-08-05-2015\preprocessed\T2-Stimulus-MOT-4-Balls\Run-1-1\FSL\4D-original\',load_nii.dat.fname);

%% MOT 4 - RUN 1 - AFTER MELODIC - BACK TO NATIVE SPACE
load_nii = nifti('filtered_func_data_back.nii');
load_nii.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\Jan-08-05-2015\preprocessed\T2-Stimulus-MOT-4-Balls\Run-1-1\FSL\Melodic-Fieldmap.ica\',load_nii.dat.fname);

%% MOT 4 - RUN 1 - MCFLIRT + UNWARP
%load_nii = nifti('20150508_100617BOLDfMRI2x2x303p3run1s003a001_unwarp.nii');
%load_nii.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\Jan-08-05-2015\preprocessed\T2-Stimulus-MOT-4-Balls\Run-1-1\FSL\4D-original\',load_nii.dat.fname);

%% MOT 4 - RUN 1 - AFTER MELODIC - BEFORE WARP
%load_nii = nifti('filtered_func_data.nii');
%load_nii.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\Jan-08-05-2015\preprocessed\T2-Stimulus-MOT-4-Balls\Run-1-1\FSL\Melodic-Fieldmap.ica\',load_nii.dat.fname);

Run = load_nii.dat(:,:,:,:);

nVoxels = zeros(1,nROI);
nTR = size(Run,4);

for iROI=1:nROI
    
    idx_ROI = AAL_ROI(idx_Temporal_Poles(iROI)).ID;
    label = AAL_ROI(idx_Temporal_Poles(iROI)).Nom_L;
    labels(iROI).label = strrep(label,'_','-');
    
    disp(labels(iROI).label);
    
    idx_voxels(iROI).idx = find(AAL_img == idx_ROI);
    
    nVoxels(iROI) = length(idx_voxels(iROI).idx);
    
    area(iROI).voxels = zeros(nVoxels(iROI),nTR);
    
    for iVoxel=1:nVoxels(iROI)
        
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iROI).idx(iVoxel));

        area(iROI).voxels(iVoxel,:) = Run(idxx,idxy,idxz,:);
        
    end
    
    all_voxels = area(iROI).voxels;
    
    area(iROI).mean = mean(all_voxels,1);
    
    clear all_voxels

end

f = figure;

for iROI=1:nROI
    
subplot(2,2,iROI);

plot(area(iROI).mean);

hold on

xlim([0 nTR]);

title(labels(iROI).label);

end