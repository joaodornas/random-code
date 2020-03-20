
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_Temporal_Poles = [83 84 87 88];

for iArea=1:length(idx_Temporal_Poles)
    
    ROI(iArea).idx_region = AAL_ROI(idx_Temporal_Poles(iArea)).ID;

    ROI(iArea).idx_voxels = find(AAL_img == ROI(iArea).idx_region);
    nVoxels = length(ROI(iArea).idx_voxels);
   
end

col = colormap(jet);

AAL_img(:,:,:) = 0;

for iArea=1:length(ROI)
   
    idx_voxels = ROI(iArea).idx_voxels;
    
    for iVoxel=1:length(idx_voxels)
        
        AAL_img(idx_voxels(iVoxel)) = iArea*100;
        
    end
        
end

fname = strcat('Temporal-RGB','.nii');
scl_slope = 1;
scl_inter = 0; 
dim = [size(AAL_img,1),size(AAL_img,2),size(AAL_img,3)];
dtype = 'FLOAT32';
offset = 0;
descrip = 'AAL - RGB';
nifti_file.mat = load_aal.mat;
nifti_file.mat_intent = load_aal.mat_intent;
nifti_file.mat0 = load_aal.mat0;
nifti_file.mat0_intent = load_aal.mat0_intent;

input_data = AAL_img; 

lowhigh_save_image;