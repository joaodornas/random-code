function real_FC_voxel_AAL_ROI(settings)

% settings_subj1_2210;
% 
% idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26 27 28];
% idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
% idx_parietal = [17 18 19 20 57 58 59 60 61 62 63 64 65 66 67 68 69 70];
% 
% doRegionByRegionAllROI(settings);
% 
% end
% 
% function doRegionByRegionAllROI(settings)

subject = settings.folders.subject;

%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

real_load_all_data_FSL;

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

for iROI=1:nNodes
    
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(strcat(int2str(iROI),':',label_ROI{iROI},':',int2str(length(idx_voxels)),':voxels'));
    
    start = now;
    disp(strcat('High Attention:',datestr(start)));
    [FC_Voxels.run(1).rho_Track, FC_Voxels.run(1).pval_Track] = getVoxelCorrelations_gpu(Track(1).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(2).rho_Track, FC_Voxels.run(2).pval_Track] = getVoxelCorrelations_gpu(Track(2).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(3).rho_Track, FC_Voxels.run(3).pval_Track] = getVoxelCorrelations_gpu(Track(3).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(4).rho_Track, FC_Voxels.run(4).pval_Track] = getVoxelCorrelations_gpu(Track(4).run,AAL_img,AAL_ROI(iROI).ID);
    finish = now;
    disp(strcat('lasted:',datestr(finish-start)));
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Track-',area_label,'.mat'),'subject','file','area_label','FC_Voxels');
    
    clear FC_Voxels
    
    start = now;
    disp(strcat('Low Attention:',datestr(start)));
    [FC_Voxels.run(1).rho_Passive, FC_Voxels.run(1).pval_Passive] = getVoxelCorrelations_gpu(Passive(1).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(2).rho_Passive, FC_Voxels.run(2).pval_Passive] = getVoxelCorrelations_gpu(Passive(2).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(3).rho_Passive, FC_Voxels.run(3).pval_Passive] = getVoxelCorrelations_gpu(Passive(3).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(4).rho_Passive, FC_Voxels.run(4).pval_Passive] = getVoxelCorrelations_gpu(Passive(4).run,AAL_img,AAL_ROI(iROI).ID);
    finish = now;
    disp(strcat('lasted:',datestr(finish-start)));
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Passive-',area_label,'.mat'),'subject','file','area_label','FC_Voxels');
    
    clear FC_Voxels
    
    start = now;
    disp(strcat('Trials:',datestr(start)));
    [FC_Voxels.run(1).rho_Trials, FC_Voxels.run(1).pval_Trials] = getVoxelCorrelations_gpu(Trials(1).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(2).rho_Trials, FC_Voxels.run(2).pval_Trials] = getVoxelCorrelations_gpu(Trials(2).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(3).rho_Trials, FC_Voxels.run(3).pval_Trials] = getVoxelCorrelations_gpu(Trials(3).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(4).rho_Trials, FC_Voxels.run(4).pval_Trials] = getVoxelCorrelations_gpu(Trials(4).run,AAL_img,AAL_ROI(iROI).ID);
    finish = now;
    disp(strcat('lasted:',datestr(finish-start)));
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Trials-',area_label,'.mat'),'subject','file','area_label','FC_Voxels');
    
    clear FC_Voxels
    
    start = now;
    disp(strcat('Resting State:',datestr(start)));
    [FC_Voxels.run(1).rho_RestingState, FC_Voxels.run(1).pval_RestingState] = getVoxelCorrelations_gpu(RestingState(1).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(2).rho_RestingState, FC_Voxels.run(2).pval_RestingState] = getVoxelCorrelations_gpu(RestingState(2).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(3).rho_RestingState, FC_Voxels.run(3).pval_RestingState] = getVoxelCorrelations_gpu(RestingState(3).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(4).rho_RestingState, FC_Voxels.run(4).pval_RestingState] = getVoxelCorrelations_gpu(RestingState(4).run,AAL_img,AAL_ROI(iROI).ID);
    finish = now;
    disp(strcat('lasted:',datestr(finish-start)));
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-RestingState-',area_label,'.mat'),'subject','file','area_label','FC_Voxels');

    clear FC_Voxels
    
end

end
