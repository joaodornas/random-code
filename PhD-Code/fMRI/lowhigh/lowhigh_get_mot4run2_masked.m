

settings_jan_0805;

lowhigh_load_all_data_FSL;

idx_space = find(mask_MOT4Run2 ~= 1);

nSpaceVoxels = length(idx_space);
nTR = size(MOT4Run2,4);

for iSpaceVoxel=1:nSpaceVoxels
    
    [idxx,idxy,idxz] = ind2sub(size(mask_MOT4Run2),idx_space(iSpaceVoxel));
    
    MOT4Run2(idxx,idxy,idxz,:) = 0;
    
end

run = 1;

dim = [size(MOT4Run2,1), size(MOT4Run2,2), size(MOT4Run2,3), size(MOT4Run2,4)];

nifti_file = settings.mot4.FSL.run(run).nifti_file;
dtype = settings.mot4.FSL.run(run).dtype;
offset = settings.mot4.FSL.run(run).offset;
scl_slope = settings.mot4.FSL.run(run).scl_slope;
scl_inter = settings.mot4.FSL.run(run).scl_inter;
    
dtype = 'FLOAT32';
offset = 0;

descrip = 'MOT4Run2 masked';
 
fname = strcat('Low-High-',settings.folders.subject,'-MOT4-Run2-','masked.nii');

input_data = MOT4Run2; 

lowhigh_save_image;     