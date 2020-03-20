function real_FC_voxel_AAL_ROI_densities

% settings_subj1_2210;
% doTheMath(settings,idx_ROI);

% idx_ROI = 19:90;
% settings_subj2_2610;
% doTheMath(settings,idx_ROI);

nROI = 90;
idx_ROI = 1:nROI;

% settings_subj3_0311;
% doTheMath(settings,idx_ROI);
% 
% settings_subj4_0211;
% doTheMath(settings,idx_ROI);
 
% settings_subj5_0211;
% doTheMath(settings,idx_ROI);
%  
% settings_subj6_2411;
% doTheMath(settings,idx_ROI);
% 
% settings_subj7_1401;
% doTheMath(settings,idx_ROI);

settings_subj8_1401;
doTheMath(settings,idx_ROI);

end

function doTheMath(settings,idx_ROI)

subject = settings.folders.subject;

disp(strcat(subject,'-',datestr(now)));

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;
nRun = 4;

for iROI=idx_ROI
    
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(strcat(int2str(iROI),':',label_ROI{iROI},':',int2str(length(idx_voxels)),':voxels'));
    
    start = now;
    disp(strcat('High Attention:',datestr(start)));
    load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-Track-',area_label,'.mat'));
    for irun=1:nRun
        [dAll_mat, dAll_pos_mat, dAll_neg_mat, FC_Density.run(irun).voxels_all_density, FC_Density.run(irun).voxels_pos_density, FC_Density.run(irun).voxels_neg_density] = computeDenstitiesOfCorrelations(FC_Voxels.run(irun).rho_Track,FC_Voxels.run(irun).pval_Track);
    end
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-density-Track-',area_label,'.mat'),'FC_Density');
    clear FC_Voxels
    clear FC_Density
    
    start = now;
    disp(strcat('Low Attention:',datestr(start)));
    load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-Passive-',area_label,'.mat'));
    for irun=1:nRun
        [dAll_mat, dAll_pos_mat, dAll_neg_mat, FC_Density.run(irun).voxels_all_density, FC_Density.run(irun).voxels_pos_density, FC_Density.run(irun).voxels_neg_density] = computeDenstitiesOfCorrelations(FC_Voxels.run(irun).rho_Passive,FC_Voxels.run(irun).pval_Passive);
    end
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-density-Passive-',area_label,'.mat'),'FC_Density');
    clear FC_Voxels
    clear FC_Density
    
    start = now;
    disp(strcat('Trials:',datestr(start)));
    load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-Trials-',area_label,'.mat'));
    for irun=1:nRun
        [dAll_mat, dAll_pos_mat, dAll_neg_mat, FC_Density.run(irun).voxels_all_density, FC_Density.run(irun).voxels_pos_density, FC_Density.run(irun).voxels_neg_density] = computeDenstitiesOfCorrelations(FC_Voxels.run(irun).rho_Trials,FC_Voxels.run(irun).pval_Trials);
    end
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-density-Trials-',area_label,'.mat'),'FC_Density');
    clear FC_Voxels
    clear FC_Density
    
    start = now;
    disp(strcat('RestingState:',datestr(start)));
    load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-RestingState-',area_label,'.mat'));
    for irun=1:nRun
        [dAll_mat, dAll_pos_mat, dAll_neg_mat, FC_Density.run(irun).voxels_all_density, FC_Density.run(irun).voxels_pos_density, FC_Density.run(irun).voxels_neg_density] = computeDenstitiesOfCorrelations(FC_Voxels.run(irun).rho_RestingState,FC_Voxels.run(irun).pval_RestingState);
    end
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-density-RestingState-',area_label,'.mat'),'FC_Density');
    clear FC_Voxels
    clear FC_Density
 
end

end
