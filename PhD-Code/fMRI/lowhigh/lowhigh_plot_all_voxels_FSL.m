function lowhigh_plot_all_voxels_FSL

%settings_jan_0805;
settings_elena_2905;

lowhigh_load_all_data_FSL;

run = 1;
plotAllVoxels(settings,MOT4Run1,MOT2Run1,RestingStateRun1,run);
run = 2;
plotAllVoxels(settings,MOT4Run2,MOT2Run2,RestingStateRun2,run);

return


function plotAllVoxels(settings,MOT4,MOT2,RestingState,run)

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

%%% Plot Voxels

nRegions = 90;

colors{1} = 'b';
colors{2} = 'r';
colors{3} = 'y';

for iRegion=1:nRegions
    
    idx_region = AAL_ROI(iRegion).ID;
    
    label_AAL = AAL_ROI(iRegion).Nom_L;
   
    idx_voxels_structures = find(AAL_img == idx_region);
    nVoxels = length(idx_voxels_structures);
    
    disp(strcat(label_AAL,':',num2str(nVoxels),'(voxels)'));
    
    for iVoxel=1:nVoxels
        
        [idxx(iVoxel), idxy(iVoxel), idxz(iVoxel)] = ind2sub(size(AAL_img),idx_voxels_structures(iVoxel));
    
    end
    
    nTR = size(MOT4,4);
    
    MOT4_voxel_time_series = zeros(nVoxels,nTR);
    MOT2_voxel_time_series = zeros(nVoxels,nTR);
    RestingState_voxel_time_series = zeros(nVoxels,nTR);
    
    for iVoxel=1:nVoxels
        
        MOT4_voxel_time_series(iVoxel,1:nTR) = MOT4(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        MOT2_voxel_time_series(iVoxel,1:nTR) = MOT2(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        RestingState_voxel_time_series(iVoxel,1:nTR) = RestingState(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
    
    end
    
    MOT4_mean_time = mean(MOT4_voxel_time_series,2);
    MOT2_mean_time = mean(MOT2_voxel_time_series,2);
    RestingState_mean_time = mean(RestingState_voxel_time_series,2);
        
    [m, I] = sort(MOT4_mean_time);
    MOT4_idx_voxel = I(end);
    
    [m, I] = sort(MOT2_mean_time);
    MOT2_idx_voxel = I(end);
    
    [m, I] = sort(RestingState_mean_time);
    RestingState_idx_voxel = I(end);
    
    disp(strcat('idx_voxels:',int2str(MOT4_idx_voxel),':',int2str(MOT2_idx_voxel),':',int2str(RestingState_idx_voxel)));
    
    voxel(iRegion).MOT4_voxel = squeeze(MOT4_voxel_time_series(MOT4_idx_voxel,:));
    voxel(iRegion).MOT2_voxel = squeeze(MOT2_voxel_time_series(MOT2_idx_voxel,:));
    voxel(iRegion).RestingState_voxel = squeeze(RestingState_voxel_time_series(RestingState_idx_voxel,:));
    
end

for iRegion=1:nRegions
    
    max_MOT4(iRegion) = max(voxel(iRegion).MOT4_voxel(:));
    max_MOT2(iRegion) = max(voxel(iRegion).MOT2_voxel(:));
    max_RestingState(iRegion) = max(voxel(iRegion).RestingState_voxel(:));
    
    min_MOT4(iRegion) = min(voxel(iRegion).MOT4_voxel(:));
    min_MOT2(iRegion) = min(voxel(iRegion).MOT2_voxel(:));
    min_RestingState(iRegion) = min(voxel(iRegion).RestingState_voxel(:));
    
end

for iRegion=1:nRegions

    f = figure;

    time_series_MOT4 = voxel(iRegion).MOT4_voxel(:);
    plot(time_series_MOT4,colors{1});

    hold on

    time_series_MOT2 = voxel(iRegion).MOT2_voxel(:);
    plot(time_series_MOT2,colors{2});

    time_series_RestingState = voxel(iRegion).RestingState_voxel(:);
    plot(time_series_RestingState,colors{3});

    legend({'High Attention', 'Low Attention', 'Resting State'});

    label = AAL_ROI(iRegion).Nom_L;
   
    label = strrep(label,'_','-');
    
    title(strcat(label,'-','Run','-',int2str(run)));

    ylabel('BOLD activity');

    xlabel('TR');

    axis = gca;
    axis.xtick = 1:nTR;

    ylim([min(min([time_series_MOT4, time_series_MOT2, time_series_RestingState])) max(max([time_series_MOT4, time_series_MOT2, time_series_RestingState]))]);

    print(f,'-djpeg',strcat('Low-High-',settings.folders.subject,'-',label,'-voxels-Run-',int2str(run),'.jpeg'));
    print(f,'-depsc',strcat('Low-High-',settings.folders.subject,'-',label,'-voxels-Run-',int2str(run),'.eps'));
    
    close all

end

return



