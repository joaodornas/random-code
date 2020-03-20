
settings_subj1_2210;
% settings_subj2_2610;
% settings_subj3_0311;
% settings_subj4_0211;
% settings_subj5_0211;
% settings_subj6_2411;

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRuns = 4;

idx_ROI = 1:90;

step = 0.001;
interval = -1:step:1;

mid = round(length(interval)/2);

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    Track = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-','distribution','.mat'));
    Passive = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-','distribution','.mat'));
    RestingState = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-','distribution','.mat'));

    for irun=1:nRuns
        
        disp('Track-Passive');
        dist1 = Track.TrackDistribution.run(irun).probability_distribution;
        dist2 = Passive.PassiveDistribution.run(irun).probability_distribution;
        distance = getHelingerDistance(dist1,dist2);
        
        pos1 = dist1(mid:end);
        int_pos1 = sum(pos1);
        pos2 = dist2(mid:end);
        int_pos2 = sum(pos2);
        
        neg1 = dist1(1:mid);
        int_neg1 = sum(neg1);
        neg2 = dist2(1:mid);
        int_neg2 = sum(neg2);
        
        shift_pos = int_pos1-int_pos2;
        shift_neg = int_neg1-int_neg2;
        
        TrackPassiveHelinger(iROI,irun) = distance;
        TrackPassiveShiftPos(iROI,irun) = shift_pos;
        TrackPassiveShiftNeg(iROI,irun) = shift_neg;
        
        disp('Track-RestingState');
        dist1 = Track.TrackDistribution.run(irun).probability_distribution;
        dist2 = RestingState.RestingStateDistribution.run(irun).probability_distribution;
        distance = getHelingerDistance(dist1,dist2);
        
        pos1 = dist1(mid:end);
        int_pos1 = sum(pos1);
        pos2 = dist2(mid:end);
        int_pos2 = sum(pos2);
        
        neg1 = dist1(1:mid);
        int_neg1 = sum(neg1);
        neg2 = dist2(1:mid);
        int_neg2 = sum(neg2);
        
        shift_pos = int_pos1-int_pos2;
        shift_neg = int_neg1-int_neg2;
        
        TrackRestingStateHelinger(iROI,irun) = distance;
        TrackRestingStateShiftPos(iROI,irun) = shift_pos;
        TrackRestingStateShiftNeg(iROI,irun) = shift_neg;
        
        disp('Passive-RestingState');
        dist1 = Passive.PassiveDistribution.run(irun).probability_distribution;
        dist2 = RestingState.RestingStateDistribution.run(irun).probability_distribution;
        distance = getHelingerDistance(dist1,dist2);
        
        pos1 = dist1(mid:end);
        int_pos1 = sum(pos1);
        pos2 = dist2(mid:end);
        int_pos2 = sum(pos2);
        
        neg1 = dist1(1:mid);
        int_neg1 = sum(neg1);
        neg2 = dist2(1:mid);
        int_neg2 = sum(neg2);
        
        shift_pos = int_pos1-int_pos2;
        shift_neg = int_neg1-int_neg2;
        
        PassiveRestingStateHelinger(iROI,irun) = distance;
        PassiveRestingStateShiftPos(iROI,irun) = shift_pos;
        PassiveRestingStateShiftNeg(iROI,irun) = shift_neg;
        
    end
    
end