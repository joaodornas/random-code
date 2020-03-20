
settings_jan_0805;
%settings_elena_2905;

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

lowhigh_load_all_data_FSL;

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

% disp('Frontal');
% idx_nodes = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];
% 
% disp('Occipital');
% idx_nodes = [49 50 51 52 53 54];
% 
% disp('Parietal');
% idx_nodes = [59 60 61 62];
% 
% disp('Temporal');
% idx_nodes = [81 82 83 84 85 86 87 88 89 90];

idx_nodes = [13 14 49 50 65 66];

iside = 3;
iiVoxels = iside^3;

for iNode=1:length(idx_nodes)
   
    idx = idx_nodes(iNode);
    
    idx_AAL = AAL_ROI(idx).ID;
    
    label_AAL = AAL_ROI(idx).Nom_L;
    
    idx_voxels = find(AAL_img==idx_AAL);
    
    nVoxels = length(idx_voxels);
    
    disp(strcat(label_AAL,':',num2str(nVoxels),'(voxels)'));
    
    for iVoxels=1:nVoxels
        
        idx = idx_voxels(iVoxels);
        
        [idxx(iVoxels), idxy(iVoxels), idxz(iVoxels)] = ind2sub(size(AAL_img),idx);
        
    end
    
    nTR = size(MOT4Run1,4);
    
    MOT4Run1_cube = zeros(iside,iside,iside,nTR);
    MOT2Run1_cube = zeros(iside,iside,iside,nTR);
    RestingStateRun1_cube = zeros(iside,iside,iside,nTR);
   
    MOT4Run2_cube = zeros(iside,iside,iside,nTR);
    MOT2Run2_cube = zeros(iside,iside,iside,nTR);
    RestingStateRun2_cube = zeros(iside,iside,iside,nTR);
    
    for iTR=1:nTR
        
       for iVoxels=1:iiVoxels
           
           MOT4Run1_voxel(iVoxels) = MOT4Run1(idxx(iVoxels),idxy(iVoxels),idxz(iVoxels),iTR);
           MOT2Run1_voxel(iVoxels) = MOT2Run1(idxx(iVoxels),idxy(iVoxels),idxz(iVoxels),iTR);
           RestingStateRun1_voxel(iVoxels) = RestingStateRun1(idxx(iVoxels),idxy(iVoxels),idxz(iVoxels),iTR);
           
           MOT4Run2_voxel(iVoxels) = MOT4Run2(idxx(iVoxels),idxy(iVoxels),idxz(iVoxels),iTR);
           MOT2Run2_voxel(iVoxels) = MOT2Run2(idxx(iVoxels),idxy(iVoxels),idxz(iVoxels),iTR);
           RestingStateRun2_voxel(iVoxels) = RestingStateRun2(idxx(iVoxels),idxy(iVoxels),idxz(iVoxels),iTR);
        
       end
        
       ROI(iNode).MOT4Run1_cube(:,:,:,iTR) = reshape(MOT4Run1_voxel,[iside,iside,iside]);
       ROI(iNode).MOT2Run1_cube(:,:,:,iTR) = reshape(MOT2Run1_voxel,[iside,iside,iside]);
       ROI(iNode).RestingStateRun1_cube(:,:,:,iTR) = reshape(RestingStateRun1_voxel,[iside,iside,iside]);
       
       ROI(iNode).MOT4Run2_cube(:,:,:,iTR) = reshape(MOT4Run2_voxel,[iside,iside,iside]);
       ROI(iNode).MOT2Run2_cube(:,:,:,iTR) = reshape(MOT2Run2_voxel,[iside,iside,iside]);
       ROI(iNode).RestingStateRun2_cube(:,:,:,iTR) = reshape(RestingStateRun2_voxel,[iside,iside,iside]);
       
       ROI(iNode).label_AAL = label_AAL;
       
    end
    
 end

