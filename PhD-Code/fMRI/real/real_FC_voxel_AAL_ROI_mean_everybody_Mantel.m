
%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_ROI = 1:90;
nROI = length(idx_ROI);

u = 100;

for iROI=1:nROI
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);

    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'.mat'));
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'.mat'));
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'.mat'));

    disp('Mantel - TP');
    z_TP(iROI) = qap(all_Track,all_Passive,u);
 
    disp('Mantel - TR');
    z_TR(iROI) = qap(all_Track,all_RestingState,u);
    
    disp('Mantel - PR');
    z_PR(iROI) = qap(all_Passive,all_RestingState,u);
    
end

save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Mantel','.mat'),'z_PR','z_TP','z_TR');

