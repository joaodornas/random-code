
%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_ROI = 1:90;

for iROI=1:length(idx_ROI)

    label_ROI = AAL_ROI(idx_ROI(iROI)).Nom_L;

    area_label = strrep(label_ROI,'_','-');       

    disp(area_label);

    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-','distribution','.mat'));
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-','distribution','.mat'));
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-','distribution','.mat'));

    [H_TP(iROI),P_TP(iROI)] = ttest(TrackDistribution,PassiveDistribution);
    [H_TR(iROI),P_TR(iROI)] = ttest(TrackDistribution,RestingStateDistribution);
    [H_PR(iROI),P_PR(iROI)] = ttest(PassiveDistribution,RestingStateDistribution);
    
end

pcriterion = 0.01/3;

P_TP_sig = P_TP < pcriterion;
P_TR_sig = P_TR < pcriterion;
P_PR_sig = P_PR < pcriterion;

if size(P_TP_sig,1) == 1, P_TP_sig = P_TP_sig'; end
if size(P_TR_sig,1) == 1, P_TR_sig = P_TR_sig'; end
if size(P_PR_sig,1) == 1, P_PR_sig = P_PR_sig'; end

all_P_sig = [P_TP_sig,P_TR_sig,P_PR_sig];

P_OnlyTrack_sig = P_TP_sig & ~P_PR_sig;