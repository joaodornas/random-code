

%% MNI SPACE
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_nodes_frontal = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];
idx_nodes_occipital = [49 50 51 52 53 54];
idx_nodes_parietal = [59 60 61 62];
idx_nodes_temporal = [81 82 83 84 85 86 87 88 89 90];

idx_nodes = [idx_nodes_frontal, idx_nodes_occipital, idx_nodes_parietal, idx_nodes_temporal];

nROI = length(idx_nodes);

%% LOAD DATA

%% JAN-08-05-2015 - MOT 4 - RUN 2 - CUSTOM FILTERED
load_nii = nifti('filtered_func_data_mcf_unwarp2standard.nii');
load_nii.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\Jan-08-05-2015\preprocessed\T2-Stimulus-MOT-4-Balls\Run-2-4\FSL\custom\',load_nii.dat.fname);

filtered = load_nii.dat(:,:,:,:);

clear load_nii

%% JAN-08-05-2015 - MOT 4 - RUN 2 - CUSTOM FILTERED WAVELET VOXEL LEVEL
load_nii = nifti('filtered_func_data_mcf_unwarp2standard-clean-voxel.nii');
load_nii.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\Jan-08-05-2015\preprocessed\T2-Stimulus-MOT-4-Balls\Run-2-4\FSL\custom\',load_nii.dat.fname);

waveletvoxel = load_nii.dat(:,:,:,:);

clear load_nii

% %% JAN-08-05-2015 - MOT 4 - RUN 2 - D8
% load_nii = nifti('filtered_func_data2standard_despike_wds_d8.nii');
% load_nii.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\Jan-08-05-2015\preprocessed\T2-Stimulus-MOT-4-Balls\Run-2-4\FSL\warped\',load_nii.dat.fname);
% 
% despiked8 = load_nii.dat(:,:,:,:);
% 
% clear load_nii

nVoxels = zeros(1,nROI);
nTR = 331;

for iROI=1:nROI
     
    idx_ROI = AAL_ROI(idx_nodes(iROI)).ID;
    label = AAL_ROI(idx_nodes(iROI)).Nom_L;
    labels(iROI).label = strrep(label,'_','-');
    
    disp(labels(iROI).label);
    
    idx_voxels(iROI).idx = find(AAL_img == idx_ROI);
    
    nVoxels(iROI) = length(idx_voxels(iROI).idx);
    
    melodic_area(iROI).voxels = zeros(nVoxels(iROI),nTR);
    
    %despiked8_area(iROI).voxels = zeros(nVoxels(iROI),nTR);
    
    despiked_area(iROI).voxels = zeros(nVoxels(iROI),nTR);
    
    for iVoxel=1:nVoxels(iROI)
        
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iROI).idx(iVoxel));

        melodic_area(iROI).voxels(iVoxel,:) = filtered(idxx,idxy,idxz,:);
    
        %despiked8_area(iROI).voxels(iVoxel,:) = despiked8(idxx,idxy,idxz,:); 
        
        despiked_area(iROI).voxels(iVoxel,:) = waveletvoxel(idxx,idxy,idxz,:);
        
    end
    
    all_voxels = melodic_area(iROI).voxels;
    melodic_area(iROI).mean = mean(all_voxels,1);
    clear all_voxels
    
%     all_voxels = despiked8_area(iROI).voxels;
%     despiked8_area(iROI).mean = mean(all_voxels,1);
%     clear all_voxels
    
    all_voxels = despiked_area(iROI).voxels;
    despiked_area(iROI).mean = mean(all_voxels,1);
    clear all_voxels

end

f = figure;

%for iROI=1:nROI
 for iROI=1:20
     
subplot(5,4,iROI);

plot(zscore(melodic_area(iROI).mean),'r');

%hold on

%plot(despiked8_area(iROI).mean,'r');
plot(zscore(despiked_area(iROI).mean),'b');

xlim([0 nTR]);

title(labels(iROI).label);

%legend({'Melodic','Wavelet D8'});
%legend({'Melodic','Wavelet'});

end