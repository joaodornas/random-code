function lowhigh_std_time_FSL

%settings_jan_0805;
settings_elena_2905;

doTheMath(settings);



return


function doTheMath(settings)

lowhigh_load_all_data_FSL;

MOT4 = ( MOT4Run1 + MOT4Run2 ) ./ 2;
MOT2 = ( MOT2Run1 + MOT2Run2 ) ./ 2;
RestingState = ( RestingStateRun1 + RestingStateRun2 ) ./ 2;

masks(1).mask = mask_MOT4Run1;
masks(2).mask = mask_MOT4Run2;
masks(3).mask = mask_MOT2Run1;
masks(4).mask = mask_MOT2Run2;
masks(5).mask = mask_RestingStateRun1;
masks(6).mask =  mask_RestingStateRun2;

common_mask = get_common_mask(masks);

idx_space = find(common_mask == 0);

MOT4std = getStd(MOT4,common_mask);
MOT2std = getStd(MOT2,common_mask);
RestingStatestd = getStd(RestingState,common_mask);

MOT4std(idx_space) = 0;
MOT2std(idx_space) = 0;
RestingStatestd(idx_space) = 0;

saveAllStd(MOT4std,MOT2std,RestingStatestd,settings);

return

function img_std = getStd(img,common_mask)

idx_brain = find(common_mask);

nVoxels = length(idx_brain);
nTR = size(img,4);
voxels_time_series = zeros(nTR,nVoxels);

for iVoxel=1:nVoxels
    
    [idxx(iVoxel) idxy(iVoxel) idxz(iVoxel)] = ind2sub(size(common_mask),idx_brain(iVoxel));

    time_series = img(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
    
    if ~isempty(find(time_series == 0)), time_series = zeros(1,length(time_series)); end
    
    voxels_time_series(:,iVoxel) = time_series;
        
end

voxels_std = std(voxels_time_series);

%take_max = max(voxels_std);
take_max = 600;

voxels_std(voxels_std > take_max/2) = 0;

img_std = zeros(size(common_mask));

img_std(idx_brain) = voxels_std(:);

return

function saveAllStd(MOT4std,MOT2std,RestingStatestd,settings)

run = 1;

dim = [settings.mot4.FSL.run(run).dim(1), settings.mot4.FSL.run(run).dim(2),settings.mot4.FSL.run(run).dim(3)];

nifti_file = settings.mot4.FSL.run(run).nifti_file;
dtype = settings.mot4.FSL.run(run).dtype;
offset = settings.mot4.FSL.run(run).offset;
scl_slope = settings.mot4.FSL.run(run).scl_slope;
scl_inter = settings.mot4.FSL.run(run).scl_inter;
    
dtype = 'FLOAT32';
offset = 0;

descrip = 'Std through Time';
 
fname = strcat('Low-High-',settings.folders.subject,'-MOT4-AllRuns-','std-time','.nii');

input_data = MOT4std; 

lowhigh_save_image;     

fname = strcat('Low-High-',settings.folders.subject,'-MOT2-AllRuns-','std-time','.nii');

input_data = MOT2std; 

lowhigh_save_image;   

fname = strcat('Low-High-',settings.folders.subject,'-RestingState-AllRuns-','std-time','.nii');

input_data = RestingStatestd; 

lowhigh_save_image;   


return
