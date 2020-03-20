
function real_check_Martin_dataset

nRuns = 4;
settings_martin;

for iRun=1:nRuns
    
    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.warped;

    [run(iRun).warped, settings] = real_get_data_FSL_Martin(settings,iRun,file,get_at_this_preprocessed_step);

    file = settings.FSL.files.functional.residual;

    [run(iRun).residual, settings] = real_get_data_FSL_Martin(settings,iRun,file,get_at_this_preprocessed_step);

end

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for iRun=1:nRuns

    run(iRun).warped_mean = getConditionMeanROI(run(iRun).warped,AAL_img,AAL_ROI);
    
    run(iRun).residual_mean = getConditionMeanROI(run(iRun).residual,AAL_img,AAL_ROI);

end

idx_ROI = 3; %%% Frontal Sup L
plotCondition(run,idx_ROI,AAL_ROI);

idx_ROI = 49; %%% Occipital Sup L
plotCondition(run,idx_ROI,AAL_ROI);

idx_ROI = 59; %%% Parietal Sup L
plotCondition(run,idx_ROI,AAL_ROI);

end

function meanCondition = getConditionMeanROI(Condition,AAL_img,AAL_ROI)

nNodes = 90;
for iNode=1:nNodes
    
   idx_region = AAL_ROI(iNode).ID;
   idx_voxels_structures = find(AAL_img == idx_region);
   nVoxels = length(idx_voxels_structures);
   
   idxx = zeros(1,nVoxels);
   idxy = zeros(1,nVoxels);
   idxz = zeros(1,nVoxels);
   
   for iVoxel=1:nVoxels
        
        [idxx(iVoxel), idxy(iVoxel), idxz(iVoxel)] = ind2sub(size(AAL_img),idx_voxels_structures(iVoxel));
    
   end  
   
   nTR = size(Condition,4);
    
   Run1_voxel_time_series = zeros(nVoxels,nTR);
    
   for iVoxel=1:nVoxels
        
        Run1_voxel_time_series(iVoxel,1:nTR) = Condition(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        
   end
    
   Run1_AAL = mean(Run1_voxel_time_series,1);

   meanCondition(iNode,:) = Run1_AAL(:);
 
end

end

function plotCondition(Subject,idx_ROI,AAL_ROI)

nSubject = length(Subject);

f = figure;

for iSubject=1:nSubject
    
    subplot(6,2,iSubject);
    
    plot(Subject(iSubject).warped_mean(idx_ROI,:),'b');
    hold on
    plot(Subject(iSubject).residual_mean(idx_ROI,:),'k');
    
    title(strcat(AAL_ROI(idx_ROI).Nom_L,':SUBJ:',int2str(iSubject)));

end

legend('Warped','Residual');

print(f,'-depsc',strcat('Check-Petra-Data-',AAL_ROI(idx_ROI).Nom_L,'.eps'));

end