function real_FC_voxel_AAL_ROI_mean_everybody

all_settings = getAllSettings;

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_ROI = 1:90;
nRuns = 4;
nROI = length(idx_ROI);

for iROI=1:nROI
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    for iset=1:length(all_settings)
        
        disp(all_settings(iset).settings.codes.subject);
    
        TrackRHOs(iset).FC = load(strcat(all_settings(iset).settings.codes.experiment,'-',all_settings(iset).settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'.mat'));
        PassiveRHOs(iset).FC = load(strcat(all_settings(iset).settings.codes.experiment,'-',all_settings(iset).settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'.mat'));
        RestingStateRHOs(iset).FC = load(strcat(all_settings(iset).settings.codes.experiment,'-',all_settings(iset).settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'.mat'));
        
    end
    
    all_Track = zeros(size(TrackRHOs(1).FC.FC_Voxels.run(1).rho_Track));
    all_Passive = zeros(size(PassiveRHOs(1).FC.FC_Voxels.run(1).rho_Passive));
    all_RestingState = zeros(size(RestingStateRHOs(1).FC.FC_Voxels.run(1).rho_RestingState));
    
    for iset=1:length(all_settings)
       
        for irun=1:nRuns
            
           all_Track = all_Track + TrackRHOs(iset).FC.FC_Voxels.run(irun).rho_Track;
           all_Passive = all_Passive + PassiveRHOs(iset).FC.FC_Voxels.run(irun).rho_Passive;
           all_RestingState = all_RestingState + RestingStateRHOs(iset).FC.FC_Voxels.run(irun).rho_RestingState;
           
        end
        
    end

    all_Track = all_Track ./ nRuns*length(all_settings);
    all_Passive = all_Passive ./ nRuns*length(all_settings);
    all_RestingState = all_RestingState ./ nRuns*length(all_settings);
    
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'.mat'),'all_Track');
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'.mat'),'all_Passive');
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'.mat'),'all_RestingState');
    
end

end

