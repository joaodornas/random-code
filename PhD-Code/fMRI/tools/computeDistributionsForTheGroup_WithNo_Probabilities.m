
function [allTrackRho,allTrackPval,allPassiveRho,allPassivePval,allRestRho,allRestPval] = computeDistributionsForTheGroup_WithNo_Probabilities
        
%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

all_settings = getAllSettings;

nRuns = 4;

idx_ROI = 1:90;

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    iAllRuns = 0;
    
    for iset=1:length(all_settings)
    
        TrackRHOs = load(strcat(all_settings(iset).settings.codes.experiment,'-',all_settings(iset).settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'.mat'));
        PassiveRHOs = load(strcat(all_settings(iset).settings.codes.experiment,'-',all_settings(iset).settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'.mat'));
        RestingStateRHOs = load(strcat(all_settings(iset).settings.codes.experiment,'-',all_settings(iset).settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'.mat'));

        for irun=1:nRuns

            disp(strcat('irun=',int2str(irun)));
            
            iAllRuns = iAllRuns + 1;

            allTrackRho(iROI,iAllRuns).rho = TrackRHOs.FC_Voxels.run(irun).rho_Track;
            allPassiveRho(iROI,iAllRuns).rho = PassiveRHOs.FC_Voxels.run(irun).rho_Passive;
            allRestRho(iROI,iAllRuns).rho = RestingStateRHOs.FC_Voxels.run(irun).rho_RestingState;

            allTrackPval(iROI,iAllRuns).pval = TrackRHOs.FC_Voxels.run(irun).pval_Track;
            allPassivePval(iROI,iAllRuns).pval = PassiveRHOs.FC_Voxels.run(irun).pval_Passive;
            allRestPval(iROI,iAllRuns).pval = RestingStateRHOs.FC_Voxels.run(irun).pval_RestingState;

        end
    
    clear TrakRHOs
    clear PassiveRHOs
    clear RestingStateRHOs
    
    end
    
end

end