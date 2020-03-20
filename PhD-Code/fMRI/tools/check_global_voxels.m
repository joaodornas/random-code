
idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26 27 28];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [17 18 19 20 57 58 59 60 61 62 63 64 65 66 67 68 69 70];

idx_all_rois = [idx_frontal,idx_occipital,idx_parietal];

nROI = length(idx_all_rois);

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

for iROI=1:nROI
    
    label_ROI{iROI} = AAL_ROI(idx_all_rois(iROI)).Nom_L;
    idx_ROI = AAL_ROI(idx_all_rois(iROI)).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    nVoxels(iROI) = length(idx_voxels);
    
end

total_voxels = sum(nVoxels);

