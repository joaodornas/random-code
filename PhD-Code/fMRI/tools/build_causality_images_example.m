
%% AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

disp('Frontal');
idx_nodes_frontal = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];

disp('Occipital');
idx_nodes_occipital = [49 50 51 52 53 54];

disp('Parietal');
idx_nodes_parietal = [59 60 61 62];

disp('Temporal');
idx_nodes_temporal = [81 82 83 84 85 86 87 88 89 90];

MOT4 = zeros(size(AAL_img));
MOT2 = zeros(size(AAL_img));
RestingState = zeros(size(AAL_img));


for iNode=1:length(idx_nodes_frontal)
    
    idx_ROI = AAL_ROI(idx_nodes_frontal(iNode)).ID;
    
    idx_voxels = find(AAL_img == idx_ROI);
    
    for iVoxel=1:length(idx_voxels)
        
        MOT4(idx_voxels(iVoxel)) = randi([31,40],1);
        MOT2(idx_voxels(iVoxel)) = randi([1,10],1);
        RestingState(idx_voxels(iVoxel)) = 5;
        
    end
    
end

for iNode=1:length(idx_nodes_occipital)
    
    idx_ROI = AAL_ROI(idx_nodes_occipital(iNode)).ID;
    
    idx_voxels = find(AAL_img == idx_ROI);
    
    for iVoxel=1:length(idx_voxels)
        
        MOT4(idx_voxels(iVoxel)) = 50;
        MOT2(idx_voxels(iVoxel)) = 50;
        RestingState(idx_voxels(iVoxel)) = 50;
        
    end
    
end

for iNode=1:length(idx_nodes_parietal)
    
    idx_ROI = AAL_ROI(idx_nodes_parietal(iNode)).ID;
    
    idx_voxels = find(AAL_img == idx_ROI);
    
    for iVoxel=1:length(idx_voxels)
        
        MOT4(idx_voxels(iVoxel)) = randi([1,10],1);
        MOT2(idx_voxels(iVoxel)) = randi([1,10],1);
        RestingState(idx_voxels(iVoxel)) = 5;
        
    end
    
end

for iNode=1:length(idx_nodes_temporal)
    
    idx_ROI = AAL_ROI(idx_nodes_temporal(iNode)).ID;
    
    idx_voxels = find(AAL_img == idx_ROI);
    
    for iVoxel=1:length(idx_voxels)
        
        MOT4(idx_voxels(iVoxel)) = randi([21,30],1);
        MOT2(idx_voxels(iVoxel)) = randi([21,30],1);
        RestingState(idx_voxels(iVoxel)) = 5;
        
    end
    
end

rendfile = 'K:\Dropbox (Uni Magdeburg)\_TOOLBOX\spm8\canonical\cortex_20484.surf.gii';

brt = 1;
min_value = 1;
max_value = 50;

setGlobalMinMax(min_value,max_value);

fname = strcat('Causality-Example-MOT4','.nii');
scl_slope = 1;
scl_inter = 0; 
dim = [size(AAL_img,1),size(AAL_img,2),size(AAL_img,3)];
dtype = 'FLOAT32';
offset = 0;
descrip = 'Causality';
nifti_file.mat = load_aal.mat;
nifti_file.mat_intent = load_aal.mat_intent;
nifti_file.mat0 = load_aal.mat0;
nifti_file.mat0_intent = load_aal.mat0_intent;

input_data = MOT4; 

lowhigh_save_image;

img_dat = fromImageToSPMdat(fname);

snapshot_label = 'Causality-Example-MOT4';
setRenderingSnapshot(true,snapshot_label);

spm_render_my(img_dat,brt,rendfile,fname);

fname = strcat('Causality-Example-MOT2','.nii');
scl_slope = 1;
scl_inter = 0; 
dim = [size(AAL_img,1),size(AAL_img,2),size(AAL_img,3)];
dtype = 'FLOAT32';
offset = 0;
descrip = 'Causality';
nifti_file.mat = load_aal.mat;
nifti_file.mat_intent = load_aal.mat_intent;
nifti_file.mat0 = load_aal.mat0;
nifti_file.mat0_intent = load_aal.mat0_intent;

input_data = MOT2; 

lowhigh_save_image;

img_dat = fromImageToSPMdat(fname);

snapshot_label = 'Causality-Example-MOT2';
setRenderingSnapshot(true,snapshot_label);

spm_render_my(img_dat,brt,rendfile,fname);

fname = strcat('Causality-Example-RestingState','.nii');
scl_slope = 1;
scl_inter = 0; 
dim = [size(AAL_img,1),size(AAL_img,2),size(AAL_img,3)];
dtype = 'FLOAT32';
offset = 0;
descrip = 'Causality';
nifti_file.mat = load_aal.mat;
nifti_file.mat_intent = load_aal.mat_intent;
nifti_file.mat0 = load_aal.mat0;
nifti_file.mat0_intent = load_aal.mat0_intent;

input_data = RestingState; 

lowhigh_save_image;

img_dat = fromImageToSPMdat(fname);

snapshot_label = 'Causality-Example-RestingState';
setRenderingSnapshot(true,snapshot_label);

spm_render_my(img_dat,brt,rendfile,fname);

 