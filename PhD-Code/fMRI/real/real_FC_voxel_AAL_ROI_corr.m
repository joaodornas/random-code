function real_FC_voxel_AAL_ROI_corr

nROI = 90;
% idx_ROI = 1:nROI;

% settings_subj1_2210;
% doTheMath(settings,idx_ROI);

% idx_ROI = 19;
% settings_subj2_2610;
% doTheMath(settings,idx_ROI);

% idx_ROI = 7:nROI;
% settings_subj3_0311;
% doTheMath(settings,idx_ROI);

% idx_ROI = 67:nROI;
% settings_subj4_0211;
% doTheMath(settings,idx_ROI);
% 
% idx_ROI = 1:nROI;
% settings_subj5_0211;
% doTheMath(settings,idx_ROI);
% 
% idx_ROI = 1:nROI;
% settings_subj6_2411;
% doTheMath(settings,idx_ROI);

% idx_ROI = 1:nROI;
% settings_subj7_1401;
% doTheMath(settings,idx_ROI);

% idx_ROI = 1:nROI;
% settings_subj8_1401;
% doTheMath(settings,idx_ROI);

% all_settings = getAllSettingsPetra;
% nSubjects = length(all_settings);
% idx_ROI = 1:nROI;
% for iSubject=1:nSubjects
%     for iROI=idx_ROI
%         doTheMathPETRA(all_settings(iSubject).settings,iROI);
%     end
% end

% all_settings = getAllSettingsHCP;
% nSubjects = length(all_settings);
% for iSubject=35:nSubjects
%     doTheMathHCP(all_settings(iSubject).settings);
% end

settings_martin;
doTheMathMartin(settings);

end

function doTheMath(settings,idx_ROI)

subject = settings.folders.subject;

disp(strcat(subject,'-',datestr(now)));

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

for iROI=idx_ROI
    
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(strcat(int2str(iROI),':',label_ROI{iROI},':',int2str(length(idx_voxels)),':voxels'));
    
    start = now;
    disp(strcat('High Attention:',datestr(start)));
    [FC_Voxels.run(1).rho_Track, FC_Voxels.run(1).pval_Track] = getVoxelCorrelations_corr(Track(1).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(2).rho_Track, FC_Voxels.run(2).pval_Track] = getVoxelCorrelations_corr(Track(2).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(3).rho_Track, FC_Voxels.run(3).pval_Track] = getVoxelCorrelations_corr(Track(3).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(4).rho_Track, FC_Voxels.run(4).pval_Track] = getVoxelCorrelations_corr(Track(4).run,AAL_img,AAL_ROI(iROI).ID);
    finish = now;
    disp(strcat('lasted:',datestr(finish-start)));
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-Track-',area_label,'.mat'),'subject','file','area_label','FC_Voxels');
    
    clear FC_Voxels
    
    start = now;
    disp(strcat('Low Attention:',datestr(start)));
    [FC_Voxels.run(1).rho_Passive, FC_Voxels.run(1).pval_Passive] = getVoxelCorrelations_corr(Passive(1).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(2).rho_Passive, FC_Voxels.run(2).pval_Passive] = getVoxelCorrelations_corr(Passive(2).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(3).rho_Passive, FC_Voxels.run(3).pval_Passive] = getVoxelCorrelations_corr(Passive(3).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(4).rho_Passive, FC_Voxels.run(4).pval_Passive] = getVoxelCorrelations_corr(Passive(4).run,AAL_img,AAL_ROI(iROI).ID);
    finish = now;
    disp(strcat('lasted:',datestr(finish-start)));
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-Passive-',area_label,'.mat'),'subject','file','area_label','FC_Voxels');
    
    clear FC_Voxels
    
    start = now;
    disp(strcat('Trials:',datestr(start)));
    [FC_Voxels.run(1).rho_Trials, FC_Voxels.run(1).pval_Trials] = getVoxelCorrelations_corr(Trials(1).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(2).rho_Trials, FC_Voxels.run(2).pval_Trials] = getVoxelCorrelations_corr(Trials(2).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(3).rho_Trials, FC_Voxels.run(3).pval_Trials] = getVoxelCorrelations_corr(Trials(3).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(4).rho_Trials, FC_Voxels.run(4).pval_Trials] = getVoxelCorrelations_corr(Trials(4).run,AAL_img,AAL_ROI(iROI).ID);
    finish = now;
    disp(strcat('lasted:',datestr(finish-start)));
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-Trials-',area_label,'.mat'),'subject','file','area_label','FC_Voxels');
    
    clear FC_Voxels
    
    start = now;
    disp(strcat('Resting State:',datestr(start)));
    [FC_Voxels.run(1).rho_RestingState, FC_Voxels.run(1).pval_RestingState] = getVoxelCorrelations_corr(RestingState(1).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(2).rho_RestingState, FC_Voxels.run(2).pval_RestingState] = getVoxelCorrelations_corr(RestingState(2).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(3).rho_RestingState, FC_Voxels.run(3).pval_RestingState] = getVoxelCorrelations_corr(RestingState(3).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(4).rho_RestingState, FC_Voxels.run(4).pval_RestingState] = getVoxelCorrelations_corr(RestingState(4).run,AAL_img,AAL_ROI(iROI).ID);
    finish = now;
    disp(strcat('lasted:',datestr(finish-start)));
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-RestingState-',area_label,'.mat'),'subject','file','area_label','FC_Voxels');

    clear FC_Voxels
    
end

end

function doTheMathPETRA(settings,idx_ROI)

subject = settings.folders.subject;
disp(strcat(subject,'-',datestr(now)));

%% LOAD DATA
run = 1;
nSegments = 3;
nTR = 150;
get_at_this_preprocessed_step = settings.FSL.folders.melodic;
file = settings.FSL.files.functional.residual;
[RestingState_run, settings] = real_get_data_FSL_Petra(settings,run,file,get_at_this_preprocessed_step);

for iSeg=1:nSegments
        
    all_segments(iSeg).run(:,:,:,:) = squeeze(RestingState_run(:,:,:,(1+(iSeg-1)*nTR):(nTR+(iSeg-1)*nTR)));

end

%%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for iROI=idx_ROI
    
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(strcat(int2str(iROI),':',label_ROI{iROI},':',int2str(length(idx_voxels)),':voxels'));
    
    start = now;
    disp(strcat('Resting State:',datestr(start)));
    [FC_Voxels.run(1).rho_RestingState, FC_Voxels.run(1).pval_RestingState] = getVoxelCorrelations_corr(all_segments(1).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(2).rho_RestingState, FC_Voxels.run(2).pval_RestingState] = getVoxelCorrelations_corr(all_segments(2).run,AAL_img,AAL_ROI(iROI).ID);
    [FC_Voxels.run(3).rho_RestingState, FC_Voxels.run(3).pval_RestingState] = getVoxelCorrelations_corr(all_segments(3).run,AAL_img,AAL_ROI(iROI).ID);
    finish = now;
    disp(strcat('lasted:',datestr(finish-start)));
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-Petra-RestingState-',area_label,'.mat'),'subject','file','area_label','FC_Voxels');

    clear FC_Voxels
    
end

end

function doTheMathHCP(settings)

subject = settings.subject;
disp(strcat(subject,'-',datestr(now)));

%% LOAD DATA
nRuns = 4;
nTR = 1200;

for iRun=1:nRuns
           
    [run, settings] = real_get_data_HCP(settings,iRun);
    

%%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for iROI=1:90
    
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(strcat(int2str(iROI),':',label_ROI{iROI},':',int2str(length(idx_voxels)),':voxels'));
    
    start = now;
    disp(strcat('Resting State:',datestr(start)));
    [FC_Voxels.run(iRun).rho_RestingState, FC_Voxels.run(iRun).pval_RestingState] = getVoxelCorrelations_corr(run,AAL_img,AAL_ROI(iROI).ID);
    finish = now;
    disp(strcat('lasted:',datestr(finish-start)));
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-HCP-RestingState-',area_label,'-',int2str(iRun),'.mat'),'FC_Voxels');

    clear FC_Voxels
    
end

clear run

end

end

function doTheMathMartin(settings)

subject = settings.subject;
disp(strcat(subject,'-',datestr(now)));

%% LOAD DATA
nRuns = 4;
nTR = 150;

for iRun=1:nRuns
    
    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;
           
    [run, settings] = real_get_data_FSL_Martin(settings,iRun,file,get_at_this_preprocessed_step);
    
%%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for iROI=1:90
    
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(strcat(int2str(iROI),':',label_ROI{iROI},':',int2str(length(idx_voxels)),':voxels'));
    
    start = now;
    disp(strcat('Resting State:',datestr(start)));
    [FC_Voxels.run(iRun).rho_RestingState, FC_Voxels.run(iRun).pval_RestingState] = getVoxelCorrelations_corr(run,AAL_img,AAL_ROI(iROI).ID);
    finish = now;
    disp(strcat('lasted:',datestr(finish-start)));
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-Martin-RestingState-',area_label,'-',int2str(iRun),'.mat'),'FC_Voxels');

    clear FC_Voxels
    
end

clear run

end

end
