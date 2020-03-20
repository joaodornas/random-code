
nTR = size(MOT4Run1,4);

MOT4Run1_voxels = zeros(nTR,size(MOT4Run1,1)*size(MOT4Run1,2)*size(MOT4Run1,3));

RestingStateRun1_voxels = zeros(nTR,size(MOT4Run1,1)*size(MOT4Run1,2)*size(MOT4Run1,3));

for iTR=1:nTR
    
   whole_brain = squeeze(MOT4Run1(:,:,:,iTR));
   
   voxels = reshape(whole_brain,[1,size(whole_brain,1)*size(whole_brain,2)*size(whole_brain,3)]);
   
   MOT4Run1_voxels(iTR,:) = voxels;
   
   whole_brain = squeeze(MOT4Run2(:,:,:,iTR));
   
   voxels = reshape(whole_brain,[1,size(whole_brain,1)*size(whole_brain,2)*size(whole_brain,3)]);
   
   MOT4Run2_voxels(iTR,:) = voxels;
   
   whole_brain = squeeze(RestingStateRun1(:,:,:,iTR));
   
   voxels = reshape(whole_brain,[1,size(whole_brain,1)*size(whole_brain,2)*size(whole_brain,3)]);
   
   RestingStateRun1_voxels(iTR,:) = voxels;
   
end

contrastMOT4Resting = MOT4Run1_voxels - RestingStateRun1_voxels;

[H_MOT4Resting,P_MOT4Resting] = ttest(contrastMOT4Resting);

contrastMOT4Run12 = MOT4Run1_voxels - MOT4Run2_voxels;

[H_MOT4Run12,P_MOT4Run12] = ttest(contrastMOT4Run12);

