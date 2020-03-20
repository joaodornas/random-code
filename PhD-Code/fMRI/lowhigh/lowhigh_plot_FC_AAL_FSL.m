
%%% SETTINGS

%settings_jan_0805;
settings_elena_2905;

%%% LOAD DATA

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

colors{1} = 'b';
colors{2} = 'r';
colors{3} = 'y';

nTR = size(MOT4Run1,4);

MOT4_AAL_Run1 = zeros(nNodes,nTR);
MOT2_AAL_Run1 = zeros(nNodes,nTR);
RestingState_AAL_Run1 = zeros(nNodes,nTR);

MOT4_AAL_Run2 = zeros(nNodes,nTR);
MOT2_AAL_Run2 = zeros(nNodes,nTR);
RestingState_AAL_Run2 = zeros(nNodes,nTR);

for iNode=1:nNodes
    
   idx_region = AAL_ROI(iNode).ID;
   label = AAL_ROI(iNode).Nom_L;
   
   label = strrep(label,'_','-');
   
   idx_voxels_structures = find(AAL_img == idx_region);
   nVoxels = length(idx_voxels_structures);
   
   for iVoxel=1:nVoxels
        
        [idxx(iVoxel), idxy(iVoxel), idxz(iVoxel)] = ind2sub(size(AAL_img),idx_voxels_structures(iVoxel));
    
   end  
    
   MOT4Run1_voxel_time_series = zeros(nVoxels,nTR);
   MOT2Run1_voxel_time_series = zeros(nVoxels,nTR);
   RestingStateRun1_voxel_time_series = zeros(nVoxels,nTR);
    
   for iVoxel=1:nVoxels
        
        MOT4Run1_voxel_time_series(iVoxel,1:nTR) = MOT4Run1(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        MOT2Run1_voxel_time_series(iVoxel,1:nTR) = MOT2Run1(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        RestingStateRun1_voxel_time_series(iVoxel,1:nTR) = RestingStateRun1(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
    
   end
    
   MOT4Run1_AAL_ROI = mean(MOT4Run1_voxel_time_series,1);
   MOT2Run1_AAL_ROI = mean(MOT2Run1_voxel_time_series,1);
   RestingStateRun1_AAL_ROI = mean(RestingStateRun1_voxel_time_series,1);
    
   MOT4Run2_voxel_time_series = zeros(nVoxels,nTR);
   MOT2Run2_voxel_time_series = zeros(nVoxels,nTR);
   RestingStateRun2_voxel_time_series = zeros(nVoxels,nTR);
    
   for iVoxel=1:nVoxels
        
        MOT4Run2_voxel_time_series(iVoxel,1:nTR) = MOT4Run2(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        MOT2Run2_voxel_time_series(iVoxel,1:nTR) = MOT2Run2(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        RestingStateRun2_voxel_time_series(iVoxel,1:nTR) = RestingStateRun2(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
    
   end
    
   MOT4Run2_AAL_ROI = mean(MOT4Run2_voxel_time_series,1);
   MOT2Run2_AAL_ROI = mean(MOT2Run2_voxel_time_series,1);
   RestingStateRun2_AAL_ROI = mean(RestingStateRun2_voxel_time_series,1);
   
   MOT4_AAL_Run1(iNode,:) = MOT4Run1_AAL_ROI(:);
   MOT2_AAL_Run1(iNode,:) = MOT2Run1_AAL_ROI(:);
   RestingState_AAL_Run1(iNode,:) = RestingStateRun1_AAL_ROI(:);
   
   MOT4_AAL_Run2(iNode,:) = MOT4Run2_AAL_ROI(:);
   MOT2_AAL_Run2(iNode,:) = MOT2Run2_AAL_ROI(:);
   RestingState_AAL_Run2(iNode,:) = RestingStateRun2_AAL_ROI(:);
   
end

MOT4_AAL_Run1 = MOT4_AAL_Run1';
MOT2_AAL_Run1 = MOT2_AAL_Run1';
RestingState_AAL_Run1 = RestingState_AAL_Run1';

MOT4_AAL_Run2 = MOT4_AAL_Run2';
MOT2_AAL_Run2 = MOT2_AAL_Run2';
RestingState_AAL_Run2 = RestingState_AAL_Run2';

[MOT4Run1_rho, MOT4Run1_pval] = corr(MOT4_AAL_Run1);
[MOT2Run1_rho, MOT2Run1_pval] = corr(MOT2_AAL_Run1);
[RestingStateRun1_rho, RestingStateRun1_pval] = corr(RestingState_AAL_Run1);

[MOT4Run2_rho, MOT4Run2_pval] = corr(MOT4_AAL_Run2);
[MOT2Run2_rho, MOT2Run2_pval] = corr(MOT2_AAL_Run2);
[RestingStateRun2_rho, RestingStateRun2_pval] = corr(RestingState_AAL_Run2);

CLOW = min([min(min(MOT4Run1_rho)), min(min(MOT4Run2_rho)), min(min(MOT2Run1_rho)), min(min(MOT2Run2_rho)), min(min(RestingStateRun1_rho)), min(min(RestingStateRun2_rho))]);
CHIGH = max([max(max(MOT4Run1_rho)), max(max(MOT4Run2_rho)), max(max(MOT2Run1_rho)), max(max(MOT2Run2_rho)), max(max(RestingStateRun1_rho)), max(max(RestingStateRun2_rho))]);

CLIM = [CLOW CHIGH];

f = figure;

imagesc(MOT4Run1_rho,CLIM);
colorbar;
title('FC matrix');
xlabel('AAL regions');
ylabel('AAL regions');

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-HighAttention-Run-1','.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-HighAttention-Run-1','.eps'));

g = figure;

imagesc(MOT2Run1_rho,CLIM);
colorbar;
title('FC matrix');
xlabel('AAL regions');
ylabel('AAL regions');

print(g,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-LowAttention-Run-1','.jpg'));
print(g,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-LowAttention-Run-1','.eps'));

h = figure;

imagesc(RestingStateRun1_rho,CLIM);
colorbar;
title('FC matrix');
xlabel('AAL regions');
ylabel('AAL regions');

print(h,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-RestingState-Run-1','.jpg'));
print(h,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-RestingState-Run-1','.eps'));

f = figure;

imagesc(MOT4Run2_rho,CLIM);
colorbar;
title('FC matrix');
xlabel('AAL regions');
ylabel('AAL regions');

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-HighAttention-Run-2','.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-HighAttention-Run-2','.eps'));

g = figure;

imagesc(MOT2Run2_rho,CLIM);
colorbar;
title('FC matrix');
xlabel('AAL regions');
ylabel('AAL regions');

print(g,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-LowAttention-Run-2','.jpg'));
print(g,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-LowAttention-Run-2','.eps'));

h = figure;

imagesc(RestingStateRun2_rho,CLIM);
colorbar;
title('FC matrix');
xlabel('AAL regions');
ylabel('AAL regions');

print(h,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-RestingState-Run-2','.jpg'));
print(h,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-RestingState-Run-2','.eps'));

    
