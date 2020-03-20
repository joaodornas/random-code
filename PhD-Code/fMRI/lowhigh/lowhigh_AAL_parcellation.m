
%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

%%% GET CONDITIONS

MOT4 = ( MOT4Run1 + MOT4Run2 ) ./ 2;
MOT2 = ( MOT2Run1 + MOT2Run2 ) ./ 2;
RestingState = ( RestingStateRun1 + RestingStateRun2 ) ./ 2;

% MOT4 = cat(4,MOT4Run1,MOT4Run2);
% MOT2 = cat(4,MOT2Run1,MOT2Run2);
% RestingState = cat(4,RestingStateRun1,RestingStateRun2);

%%% GET ROI MEAN TIME SERIES

nTR = size(MOT4,4);

mean_area_MOT4_voxel = zeros(nNodes,nTR);
mean_area_MOT2_voxel = zeros(nNodes,nTR);
mean_area_RestingState_voxel = zeros(nNodes,nTR);

for iNode=1:nNodes
    
   
    idx_ROI = AAL_ROI(iNode).ID;
    idx_ROI_Label = AAL_ROI(iNode).Nom_L;
    idx_ROI_Label = strrep(idx_ROI_Label,'_','-');
    
    idx_voxels = find(AAL_img == idx_ROI);
    
    nVoxels = length(idx_voxels);
   
    area_MOT4 = zeros(nVoxels,nTR);
    area_MOT2 = zeros(nVoxels,nTR);
    area_RestingState = zeros(nVoxels,nTR);
    
    for iVoxel=1:nVoxels
       
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
            
        area_MOT4(iVoxel,:) = MOT4(idxx,idxy,idxz,:);
        
        area_MOT2(iVoxel,:) = MOT2(idxx,idxy,idxz,:);
        
        area_RestingState(iVoxel,:) = RestingState(idxx,idxy,idxz,:);
        
    end
    
    mn_MOT4 = mean(area_MOT4,1);
    mn_MOT2 = mean(area_MOT2,1);
    mn_RestingState = mean(area_RestingState,1);
    
    mean_area_MOT4_voxel(iNode,:) = mn_MOT4(:);
    mean_area_MOT2_voxel(iNode,:) = mn_MOT2(:);
    mean_area_RestingState_voxel(iNode,:) = mn_RestingState(:);
    
end