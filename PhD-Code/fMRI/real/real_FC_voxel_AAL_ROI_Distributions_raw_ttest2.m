
%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

[allTrackRho,allTrackPval,allPassiveRho,allPassivePval,allRestRho,allRestPval] = computeDistributionsForTheGroup_WithNo_Probabilities;
    
Track_all_rhos = getRawDistribution(allTrackRho,allTrackPval);
Passive_all_rhos = getRawDistribution(allPassiveRho,allPassivePval);
Rest_all_rhos = getRawDistribution(allRestRho,allRestPval);
    
idx_ROI = 1:90;

for iROI=1:length(idx_ROI)
    
    [H_TP2(iROI),P_TP2(iROI)] = ttest2(Track_all_rhos(iROI).rhos,Passive_all_rhos(iROI).rhos);
    [H_TR2(iROI),P_TR2(iROI)] = ttest2(Track_all_rhos(iROI).rhos,Rest_all_rhos(iROI).rhos);
    [H_PR2(iROI),P_PR2(iROI)] = ttest2(Passive_all_rhos(iROI).rhos,Rest_all_rhos(iROI).rhos);
    
end

pcriterion = 0.01/3;
%pcriterin = 0.05/3;

P_TP2_sig = P_TP2 < pcriterion;
P_TR2_sig = P_TR2 < pcriterion;
P_PR2_sig = P_PR2 < pcriterion;

if size(P_TP2_sig,1) == 1, P_TP2_sig = P_TP2_sig'; end
if size(P_TR2_sig,1) == 1, P_TR2_sig = P_TR2_sig'; end
if size(P_PR2_sig,1) == 1, P_PR2_sig = P_PR2_sig'; end

all_P_sig = [P_TP2_sig,P_TR2_sig,P_PR2_sig];

P_OnlyTrack_sig = P_TP2_sig & ~P_PR2_sig;