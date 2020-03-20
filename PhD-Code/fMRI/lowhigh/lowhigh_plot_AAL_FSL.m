
function lowhigh_plot_AAL_FSL

%%% SETTINGS

settings_jan_0805;

doTheMath(settings);

clear settings
settings_elena_2905;

doTheMath(settings);

return

function doTheMath(settings)

%%% LOAD DATA

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

for iNode=1:nNodes
    
   idx_region = AAL_ROI(iNode).ID;
   label = AAL_ROI(iNode).Nom_L;
   
   label = strrep(label,'_','-');
   
   idx_voxels_structures = find(AAL_img == idx_region);
   nVoxels = length(idx_voxels_structures);
   
   for iVoxel=1:nVoxels
        
        [idxx(iVoxel), idxy(iVoxel), idxz(iVoxel)] = ind2sub(size(AAL_img),idx_voxels_structures(iVoxel));
    
   end  
   
   nTR = size(MOT4Run1,4);
    
   MOT4Run1_voxel_time_series = zeros(nVoxels,nTR);
   MOT2Run1_voxel_time_series = zeros(nVoxels,nTR);
   RestingStateRun1_voxel_time_series = zeros(nVoxels,nTR);
    
   for iVoxel=1:nVoxels
        
        MOT4Run1_voxel_time_series(iVoxel,1:nTR) = MOT4Run1(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        MOT2Run1_voxel_time_series(iVoxel,1:nTR) = MOT2Run1(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        RestingStateRun1_voxel_time_series(iVoxel,1:nTR) = RestingStateRun1(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
    
   end
    
   MOT4Run1_AAL = mean(MOT4Run1_voxel_time_series,1);
   MOT2Run1_AAL = mean(MOT2Run1_voxel_time_series,1);
   RestingStateRun1_AAL = mean(RestingStateRun1_voxel_time_series,1);
    
   MOT4Run2_voxel_time_series = zeros(nVoxels,nTR);
   MOT2Run2_voxel_time_series = zeros(nVoxels,nTR);
   RestingStateRun2_voxel_time_series = zeros(nVoxels,nTR);
    
   for iVoxel=1:nVoxels
        
        MOT4Run2_voxel_time_series(iVoxel,1:nTR) = MOT4Run2(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        MOT2Run2_voxel_time_series(iVoxel,1:nTR) = MOT2Run2(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        RestingStateRun2_voxel_time_series(iVoxel,1:nTR) = RestingStateRun2(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
    
   end
    
   MOT4Run2_AAL = mean(MOT4Run2_voxel_time_series,1);
   MOT2Run2_AAL = mean(MOT2Run2_voxel_time_series,1);
   RestingStateRun2_AAL = mean(RestingStateRun2_voxel_time_series,1);
   
   f = figure;
   
   plot(MOT4Run1_AAL,colors{1});
   hold on;
   plot(MOT2Run1_AAL,colors{2});
   hold on;
   plot(RestingStateRun1_AAL,colors{3});
   hold on;
   
   title(strcat(label,'-','Run-1'));
   
   xlabel('TR');
   ylabel('BOLD Activity');
   legend('High Attention','Low Attention','Resting State');
   
   print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','AAL-time-series','-',label,'-Run-1','.jpg'));
   print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','AAL-time-series','-',label,'-Run-1','.eps'));
   
   close all;
   
   f = figure;
   
   plot(MOT4Run2_AAL,colors{1});
   hold on;
   plot(MOT2Run2_AAL,colors{2});
   hold on;
   plot(RestingStateRun2_AAL,colors{3});
   hold on;
   
   title(strcat(label,'-','Run-2'));
   
   xlabel('TR');
   ylabel('BOLD Activity');
   legend('High Attention','Low Attention','Resting State');
   
   print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','AAL-time-series','-',label,'-Run-2','.jpg'));
   print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','AAL-time-series','-',label,'-Run-2','.eps'));
   
   close all;
   
end

return