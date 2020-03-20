
function lowhigh_percentage_increase


settings_jan_0805;

doTheMath(settings);


return

function doTheMath(settings)

lowhigh_load_all_data_FSL;

disp('get Common Mask');

masks(1).mask = mask_MOT4Run1;
masks(2).mask = mask_MOT4Run2;
masks(3).mask = mask_MOT2Run1;
masks(4).mask = mask_MOT2Run2;
masks(5).mask = mask_RestingStateRun1;
masks(6).mask =  mask_RestingStateRun2;

common_mask = get_common_mask(masks);

MOT4 = getOverallMean(MOT4Run1,MOT4Run2,common_mask);
MOT2 = getOverallMean(MOT2Run1,MOT2Run2,common_mask);
RestingState = getOverallMean(RestingStateRun1,RestingStateRun2,common_mask);

saveOverallMean(MOT4,MOT2,RestingState,settings);

MOT4RestingState = getRelativePercentage(MOT4,RestingState,common_mask);
MOT2RestingState = getRelativePercentage(MOT2,RestingState,common_mask);
MOT4MOT2 = getRelativePercentage(MOT4,MOT2,common_mask);

MOT4RestingState_labels = getLabels(MOT4RestingState);
MOT2RestingState_labels = getLabels(MOT2RestingState);
MOT4MOT2_labels = getLabels(MOT4MOT2);

save(strcat('Low-High-',settings.folders.subject,'-Relative-Percentage-Parcellation-Labels.mat'),'MOT4RestingState_labels','MOT2RestingState_labels','MOT4MOT2_labels');

MOT4RestingState = removeAreas(MOT4RestingState);
MOT2RestingState = removeAreas(MOT2RestingState);
MOT4MOT2 = removeAreas(MOT4MOT2);

Thre = 100;

MOT4RestingState(MOT4RestingState>Thre) = 0;
MOT2RestingState(MOT2RestingState>Thre) = 0;
MOT4MOT2(MOT4MOT2>Thre) = 0;

saveRelativePercentage(MOT4RestingState,MOT2RestingState,MOT4MOT2,settings);

return

function Condition_img = getOverallMean(Run1,Run2,common_mask)

mean_Run1 = mean(Run1,4,'native');
mean_Run2 = mean(Run2,4,'native');

All_Runs = zeros(size(Run1,1),size(Run1,2),size(Run1,3),2);

All_Runs(:,:,:,1) = mean_Run1;
All_Runs(:,:,:,2) = mean_Run2;

Condition_img = mean(All_Runs,4,'native');

Condition_img(common_mask == 0) = 0;

return

function percentage_img = getRelativePercentage(img1,img2,common_mask)


relative = img1./img2;

relative = relative - 1;

percentage_img = relative.*100;

percentage_img(common_mask == 0) = 0;

percentage_img(isnan(percentage_img)) = 0;
percentage_img(isinf(percentage_img)) = 0;

return

function saveOverallMean(MOT4,MOT2,RestingState,settings)

run = 1;
nifti_file = settings.restingstate.FSL.run(run).nifti_file;
dtype = settings.restingstate.FSL.run(run).dtype;
offset = settings.restingstate.FSL.run(run).offset;
scl_slope = settings.restingstate.FSL.run(run).scl_slope;
scl_inter = settings.restingstate.FSL.run(run).scl_inter;
    
dim = [settings.restingstate.FSL.run(run).dim(1), settings.restingstate.FSL.run(run).dim(2),settings.restingstate.FSL.run(run).dim(3)];

dtype = 'FLOAT32';
offset = 0;

descrip = 'Overall Mean';

fname = strcat('Low-High-',settings.folders.subject,'-MOT4-Overall-Mean','.nii');

input_data = MOT4;

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT2-Overall-Mean','.nii');

input_data = MOT2;

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-RestingState-Overall-Mean','.nii');

input_data = RestingState;

lowhigh_save_image;

return

function saveRelativePercentage(MOT4RestingState,MOT2RestingState,MOT4MOT2,settings)

run = 1;
nifti_file = settings.restingstate.FSL.run(run).nifti_file;
dtype = settings.restingstate.FSL.run(run).dtype;
offset = settings.restingstate.FSL.run(run).offset;
scl_slope = settings.restingstate.FSL.run(run).scl_slope;
scl_inter = settings.restingstate.FSL.run(run).scl_inter;
    
dim = [settings.restingstate.FSL.run(run).dim(1), settings.restingstate.FSL.run(run).dim(2),settings.restingstate.FSL.run(run).dim(3)];

dtype = 'FLOAT32';
offset = 0;

descrip = 'Relative Percentage';

fname = strcat('Low-High-',settings.folders.subject,'-MOT4-RestingState-Relative-Percentage','.nii');

input_data = MOT4RestingState;

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT2-RestingState-Relative-Percentage','.nii');

input_data = MOT2RestingState;

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT4-MOT2-Relative-Percentage','.nii');

input_data = MOT4MOT2;

lowhigh_save_image;

return

function labels = getLabels(img)

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_labels = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_labels.ROI;

nStructures = length(AAL_ROI);

labels = cell.empty;

for iStructure=1:nStructures
   
    idx_structure = AAL_ROI(1,iStructure).ID;
    label_structure = AAL_ROI(1,iStructure).Nom_L;
    
    idx_voxels_structure = find(AAL_img == idx_structure);
    n_voxels = length(idx_voxels_structure);
    
    get_voxels_structure = img(idx_voxels_structure);
    mean_voxels_structure = mean(get_voxels_structure);
    
%     get_max_voxel = max(get_voxels_structure);
%     idx_all_max_voxels = find(get_voxels_structure == get_max_voxel);
%     idx_one_max_voxel = idx_all_max_voxels(1);
%     
%     idx_max_voxel = idx_voxels_structure(idx_one_max_voxel);
%     
%     [max_voxel_x, max_voxel_y, max_voxel_z] = ind2sub(size(img),idx_max_voxel);
    
    labels{iStructure,1} = label_structure;
    
    labels{iStructure,2} = n_voxels;
    
    labels{iStructure,3} = mean_voxels_structure;
    
end

return

function img = removeAreas(img)

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

idx_ROI = [4111 4112 9001 9002 9011 9012 9021 9022 9031 9032 9041 9042 9051 9052 9061 9062 9071 9072 9081 9082 9100 9110 9120 9130 9140 9150 9160 9170];

nROI = length(idx_ROI);

for iROI=1:nROI
    
    img(AAL_img == idx_ROI(iROI)) = 0;
    
end
    
return


