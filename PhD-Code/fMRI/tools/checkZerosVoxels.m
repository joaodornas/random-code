
nROI = length(idx_ROI);

maxVoxels = 4000;

zeros_MOT4Run1 = zeros(maxVoxels,nROI);
zeros_MOT4Run2 = zeros(maxVoxels,nROI);

zeros_MOT2Run1 = zeros(maxVoxels,nROI);
zeros_MOT2Run2 = zeros(maxVoxels,nROI);

zeros_RestingStateRun1 = zeros(maxVoxels,nROI);
zeros_RestingStateRun2 = zeros(maxVoxels,nROI);

for iROI=1:nROI
    
    idx_AAL = AAL_ROI(idx_ROI(iROI)).ID;
    
    idx_voxels = find(AAL_img==idx_AAL);
    
    for iVoxel=1:length(idx_voxels)
        
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
        
        voxel = MOT4Run1(idxx,idxy,idxz,:);
        if sum(voxel) == 0
            zeros_MOT4Run1(iVoxel,iROI) = 1;
        end
        
        voxel = MOT4Run2(idxx,idxy,idxz,:);
        if sum(voxel) == 0
            zeros_MOT4Run2(iVoxel,iROI) = 1;
        end
        
        voxel = MOT2Run1(idxx,idxy,idxz,:);
        if sum(voxel) == 0
            zeros_MOT2Run1(iVoxel,iROI) = 1;
        end
        
        voxel = MOT2Run2(idxx,idxy,idxz,:);
        if sum(voxel) == 0
            zeros_MOT2Run2(iVoxel,iROI) = 1;
        end
        
        voxel = RestingStateRun1(idxx,idxy,idxz,:);
        if sum(voxel) == 0
            zeros_RestingStateRun1(iVoxel,iROI) = 1;
        end
        
        voxel = RestingStateRun2(idxx,idxy,idxz,:);
        if sum(voxel) == 0
            zeros_RestingStateRun2(iVoxel,iROI) = 1;
        end
        
    end
    
end