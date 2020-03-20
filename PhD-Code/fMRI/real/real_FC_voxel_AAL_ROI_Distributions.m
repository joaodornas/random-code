function real_FC_voxel_AAL_ROI_Distributions


% settings_subj1_2210;
% computeDistributions(settings);
% settings_subj2_2610;
% computeDistributions(settings);
% settings_subj3_0311;
% computeDistributions(settings);
% settings_subj4_0211;
% computeDistributions(settings);
% settings_subj5_0211;
% computeDistributions(settings);
% settings_subj6_2411;
% computeDistributions(settings);

% settings_subj1_2210;
% computeDistributionsZFisher(settings);
% settings_subj2_2610;
% computeDistributionsZFisher(settings);
% settings_subj3_0311;
% computeDistributionsZFisher(settings);
% settings_subj4_0211;
% computeDistributionsZFisher(settings);
% settings_subj5_0211;
% computeDistributionsZFisher(settings);
% settings_subj6_2411;
% computeDistributionsZFisher(settings);

% compute3DHelingerDistance(settings);

%%% #62-Parietal-Inf-R
% idx_ROI = 62;
% plotDistributions(settings,idx_ROI);

% settings_subj1_2210;
% compute3DHelingerDistanceAllRuns(settings);
% settings_subj2_2610;
% compute3DHelingerDistanceAllRuns(settings);
% settings_subj3_0311;
% compute3DHelingerDistanceAllRuns(settings);
% settings_subj4_0211;
% compute3DHelingerDistanceAllRuns(settings);
% settings_subj5_0211;
% compute3DHelingerDistanceAllRuns(settings);
% settings_subj6_2411;
% compute3DHelingerDistanceAllRuns(settings);

%computeDistributionsForTheGroup;
computeDistributionsForTheGroupZFisher;

end

function computeDistributions(settings)
        
%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRuns = 4;

idx_ROI = 1:90;

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    TrackRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'.mat'));
    PassiveRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'.mat'));
    RestingStateRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'.mat'));
      
    for irun=1:nRuns
       
        disp(strcat('irun=',int2str(irun)));
        
        track_rho = TrackRHOs.FC_Voxels.run(irun).rho_Track;
        passive_rho = PassiveRHOs.FC_Voxels.run(irun).rho_Passive;
        rest_rho = RestingStateRHOs.FC_Voxels.run(irun).rho_RestingState;
        
        track_pval = TrackRHOs.FC_Voxels.run(irun).pval_Track;
        passive_pval = PassiveRHOs.FC_Voxels.run(irun).pval_Passive;
        rest_pval = RestingStateRHOs.FC_Voxels.run(irun).pval_RestingState;
        
        disp('probability Distribution - Track');
        probability_distribution = getProbabilityDistribution(track_rho,track_pval);
        TrackDistribution.run(irun).probability_distribution = probability_distribution;
        
        disp('probability Distribution - Passive');
        probability_distribution = getProbabilityDistribution(passive_rho,passive_pval);
        PassiveDistribution.run(irun).probability_distribution = probability_distribution;
        
        disp('probability Distribution - RestingState');
        probability_distribution = getProbabilityDistribution(rest_rho,rest_pval);
        RestingStateDistribution.run(irun).probability_distribution = probability_distribution;
        
    end
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-','distribution','.mat'),'TrackDistribution');
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-','distribution','.mat'),'PassiveDistribution');
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-','distribution','.mat'),'RestingStateDistribution');
    
    clear TrackDistribution
    clear PassiveDistribution
    clear RestingStateDistribution
         
end

end

function computeDistributionsZFisher(settings)
        
%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRuns = 4;

idx_ROI = 1:90;

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    TrackRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'.mat'));
    PassiveRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'.mat'));
    RestingStateRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'.mat'));
      
    for irun=1:nRuns
       
        disp(strcat('irun=',int2str(irun)));
        
        track_rho = TrackRHOs.FC_Voxels.run(irun).rho_Track;
        passive_rho = PassiveRHOs.FC_Voxels.run(irun).rho_Passive;
        rest_rho = RestingStateRHOs.FC_Voxels.run(irun).rho_RestingState;
        
        track_pval = TrackRHOs.FC_Voxels.run(irun).pval_Track;
        passive_pval = PassiveRHOs.FC_Voxels.run(irun).pval_Passive;
        rest_pval = RestingStateRHOs.FC_Voxels.run(irun).pval_RestingState;
        
        disp('probability Distribution - Track');
        [interval, probability_distribution] = getProbabilityDistributionZFisher(track_rho,track_pval);
        TrackDistribution.run(irun).probability_distribution = probability_distribution;
        TrackDistribution.run(irun).interval = interval;
        
        disp('probability Distribution - Passive');
        [interval, probability_distribution] = getProbabilityDistributionZFisher(passive_rho,passive_pval);
        PassiveDistribution.run(irun).probability_distribution = probability_distribution;
        PassiveDistribution.run(irun).interval = interval;
        
        disp('probability Distribution - RestingState');
        [interval, probability_distribution] = getProbabilityDistributionZFisher(rest_rho,rest_pval);
        RestingStateDistribution.run(irun).probability_distribution = probability_distribution;
        RestingStateDistribution.run(irun).interval = interval;
        
    end
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-','distribution-ZFisher','.mat'),'TrackDistribution');
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-','distribution-ZFisher','.mat'),'PassiveDistribution');
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-','distribution-ZFisher','.mat'),'RestingStateDistribution');
    
    clear TrackDistribution
    clear PassiveDistribution
    clear RestingStateDistribution
         
end

end

function computeDistributionsForTheGroup
        
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

            allTrackRho(iAllRuns).rho = TrackRHOs.FC_Voxels.run(irun).rho_Track;
            allPassiveRho(iAllRuns).rho = PassiveRHOs.FC_Voxels.run(irun).rho_Passive;
            allRestRho(iAllRuns).rho = RestingStateRHOs.FC_Voxels.run(irun).rho_RestingState;

            allTrackPval(iAllRuns).pval = TrackRHOs.FC_Voxels.run(irun).pval_Track;
            allPassivePval(iAllRuns).pval = PassiveRHOs.FC_Voxels.run(irun).pval_Passive;
            allRestPval(iAllRuns).pval = RestingStateRHOs.FC_Voxels.run(irun).pval_RestingState;

        end
    
    end
    
    TrackDistribution = getProbabilityDistributionGroup(allTrackRho,allTrackPval);
    PassiveDistribution = getProbabilityDistributionGroup(allPassiveRho,allPassivePval);
    RestingStateDistribution = getProbabilityDistributionGroup(allRestRho,allRestPval);
    
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-','distribution','.mat'),'TrackDistribution');
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-','distribution','.mat'),'PassiveDistribution');
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-','distribution','.mat'),'RestingStateDistribution');
    
    clear TrackDistribution
    clear PassiveDistribution
    clear RestingStateDistribution
         
end

end

function computeDistributionsForTheGroupZFisher
        
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

            allTrackRho(iAllRuns).rho = TrackRHOs.FC_Voxels.run(irun).rho_Track;
            allPassiveRho(iAllRuns).rho = PassiveRHOs.FC_Voxels.run(irun).rho_Passive;
            allRestRho(iAllRuns).rho = RestingStateRHOs.FC_Voxels.run(irun).rho_RestingState;

            allTrackPval(iAllRuns).pval = TrackRHOs.FC_Voxels.run(irun).pval_Track;
            allPassivePval(iAllRuns).pval = PassiveRHOs.FC_Voxels.run(irun).pval_Passive;
            allRestPval(iAllRuns).pval = RestingStateRHOs.FC_Voxels.run(irun).pval_RestingState;

        end
    
    end
    
    [interval, TrackDistribution] = getProbabilityDistributionGroupZFisher(allTrackRho,allTrackPval);
    [interval, PassiveDistribution] = getProbabilityDistributionGroupZFisher(allPassiveRho,allPassivePval);
    [interval, RestingStateDistribution] = getProbabilityDistributionGroupZFisher(allRestRho,allRestPval);
    
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-','distribution-ZFisher','.mat'),'TrackDistribution','interval');
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-','distribution-ZFisher','.mat'),'PassiveDistribution','interval');
    save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-','distribution-ZFisher','.mat'),'RestingStateDistribution','interval');
    
    clear TrackDistribution
    clear PassiveDistribution
    clear RestingStateDistribution
         
end

end

function computeZ3DHelingerDistanceSignificantCorrelationsShift(settings)

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRuns = 4;

idx_ROI = 1:90;

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
%     TrackRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'.mat'));
%     PassiveRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'.mat'));
%     RestingStateRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'.mat'));
      
    TrackDist = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-','distribution','.mat'));
    PassiveDist = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-','distribution','.mat'));
    RestingStateDist = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-','distribution','.mat'));

    for irun=1:nRuns
       
        disp(strcat('irun=',int2str(irun)));
        
%         track_rho = TrackRHOs.FC_Voxels.run(irun).rho_Track;
%         passive_rho = PassiveRHOs.FC_Voxels.run(irun).rho_Passive;
%         rest_rho = RestingStateRHOs.FC_Voxels.run(irun).rho_RestingState;
%         
%         track_pval = TrackRHOs.FC_Voxels.run(irun).pval_Track;
%         passive_pval = PassiveRHOs.FC_Voxels.run(irun).pval_Passive;
%         rest_pval = RestingStateRHOs.FC_Voxels.run(irun).pval_RestingState;
        
        disp('Track-Passive');
        dist1 = TrackDist.TrackDistribution.run(irun).probability_distribution;
        dist2 = PassiveDist.PassiveDistribution.run(irun).probability_distribution;
        distance = getHelingerDistance(dist1,dist2);
        
        %TrackPassiveShift(iROI,irun) = getShiftSignificant(track_rho,track_pval,passive_rho,passive_pval);
        TrackPassiveDistance(iROI,irun) = distance;
        
        disp('Track-RestingState');
        dist1 = TrackDist.TrackDistribution.run(irun).probability_distribution;
        dist2 = RestingStateDist.RestingStateDistribution.run(irun).probability_distribution;
        distance = getHelingerDistance(dist1,dist2);
        
        %TrackRestingStateShift(iROI,irun) = getShiftSignificant(track_rho,track_pval,rest_rho,rest_pval);
        TrackRestingStateDistance(iROI,irun) = distance;
        
        disp('Passive-RestingState');
        dist1 = PassiveDist.PassiveDistribution.run(irun).probability_distribution;
        dist2 = RestingStateDist.RestingStateDistribution.run(irun).probability_distribution;
        distance = getHelingerDistance(dist1,dist2);
        
        %PassiveRestingStateShift(iROI,irun) = getShiftSignificant(passive_rho,passive_pval,rest_rho,rest_pval);
        PassiveRestingStateDistance(iROI,irun) = distance;    
       
    end
   
end


end

function compute3DHelingerDistanceAllSubjects(all_settings)

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRuns = 4;

nSettings = length(all_settings);

idx_ROI = 1:90;

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
         
    for iset=1:nSettings
        
        TrackDist(iset).dist = load(strcat(all_settings(iset).settings.codes.experiment,'-',all_settings(iset).settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-','distribution','.mat'));
        PassiveDist(iset).dist = load(strcat(all_settings(iset).settings.codes.experiment,'-',all_settings(iset).settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-','distribution','.mat'));
        RestingStateDist(iset).dist = load(strcat(all_settings(iset).settings.codes.experiment,'-',all_settings(iset).settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-','distribution','.mat'));

    end
    
    iirun = 0;
    
    for iset=1:nSettings
    
        for irun=1:nRuns
            
            iirun = iirun + 1;

            disp(strcat('irun=',int2str(irun)));

            disp('Track-Passive');
            dist1 = TrackDist(iset).dist.TrackDistribution.run(irun).probability_distribution;
            dist2 = PassiveDist(iset).dist.PassiveDistribution.run(irun).probability_distribution;
            distance = getHelingerDistance(dist1,dist2);

            TrackPassiveDistance(iROI,iirun) = distance;

            disp('Track-RestingState');
            dist1 = TrackDist(iset).dist.TrackDistribution.run(irun).probability_distribution;
            dist2 = RestingStateDist(iset).dist.RestingStateDistribution.run(irun).probability_distribution;
            distance = getHelingerDistance(dist1,dist2);

            TrackRestingStateDistance(iROI,iirun) = distance;

            disp('Passive-RestingState');
            dist1 = PassiveDist(iset).dist.PassiveDistribution.run(irun).probability_distribution;
            dist2 = RestingStateDist(iset).dist.RestingStateDistribution.run(irun).probability_distribution;
            distance = getHelingerDistance(dist1,dist2);

            PassiveRestingStateDistance(iROI,iirun) = distance;    

        end
    
    end
   
end


end

function compute3DHelingerDistance(settings)

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRuns = 4;

idx_ROI = 1:90;

TrackPassive = zeros(length(idx_ROI),nRuns);
PassiveRestingState = zeros(length(idx_ROI),nRuns);
TrackRestingState = zeros(length(idx_ROI),nRuns);

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
        [shift_pos_tp(iROI,irun), shift_neg_tp(iROI,irun)] = getShift(dist1,dist2);
        
        TrackPassive(iROI,irun) = distance;
        
        disp('Track-RestingState');
        dist1 = Track.TrackDistribution.run(irun).probability_distribution;
        dist2 = RestingState.RestingStateDistribution.run(irun).probability_distribution;
        distance = getHelingerDistance(dist1,dist2);
        [shift_pos_tr(iROI,irun), shift_neg_tr(iROI,irun)] = getShift(dist1,dist2);
        
        TrackRestingState(iROI,irun) = distance;
        
        disp('Passive-RestingState');
        dist1 = Passive.PassiveDistribution.run(irun).probability_distribution;
        dist2 = RestingState.RestingStateDistribution.run(irun).probability_distribution;
        distance = getHelingerDistance(dist1,dist2);
        [shift_pos_pr(iROI,irun), shift_neg_pr(iROI,irun)] = getShift(dist1,dist2);
        
        PassiveRestingState(iROI,irun) = distance;
        
    end
    
end

for irun=1:nRuns
    
    TrackPassive_img = zeros(size(AAL_img));
    PassiveRestingState_img = zeros(size(AAL_img));
    TrackRestingState_img = zeros(size(AAL_img));

    for iROI=1:length(idx_ROI)

        AAL_ID = AAL_ROI(idx_ROI(iROI)).ID;

        idx_voxels = find(AAL_img==AAL_ID);

        nVoxels = length(idx_voxels);

        for iVoxel=1:nVoxels

            [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

            if shift_pos_tp(iROI,irun)>shift_neg_tp(iROI,irun)
            
                signal = 1;
            
                if shift_pos_tp(iROI,irun) > 0, increase = 1; else increase = 0; end
            
            else
            
                signal = -1;
            
                if shift_neg_tp(iROI,irun) > 0, increase = 1; else increase = 0; end
            
            end
        
            TrackPassive_img(idxx,idxy,idxz) = increase.*signal.*TrackPassive(iROI,irun);
            
            if shift_pos_pr(iROI,irun)>shift_neg_pr(iROI,irun)
            
                signal = 1;
            
                if shift_pos_pr(iROI,irun) > 0, increase = 1; else increase = 0; end
            
            else
            
                signal = -1;
            
                if shift_neg_pr(iROI,irun) > 0, increase = 1; else increase = 0; end
            
             end
        
            PassiveRestingState_img(idxx,idxy,idxz) = increase.*signal.*PassiveRestingState(iROI,irun);
            
             if shift_pos_tr(iROI,irun)>shift_neg_tr(iROI,irun)
            
                signal = 1;
            
                if shift_pos_tr(iROI,irun) > 0, increase = 1; else increase = 0; end
            
             else
            
                signal = -1;
            
                if shift_neg_tr(iROI,irun) > 0, increase = 1; else increase = 0; end
            
             end
        
            TrackRestingState_img(idxx,idxy,idxz) = increase.*signal.*TrackRestingState(iROI,irun);

        end

    end
    
    TrackPassive_pos = TrackPassive_img;
    PassiveRestingState_pos = PassiveRestingState_img;
    TrackRestingState_pos = TrackRestingState_img;

    TrackPassive_neg = TrackPassive_img.*(-1);
    PassiveRestingState_neg = PassiveRestingState_img.*(-1);
    TrackRestingState_neg = TrackRestingState_img.*(-1);

    nifti_file = load_aal;
    offset = load_aal.dat.offset;
    scl_slope = load_aal.dat.scl_slope;
    scl_inter = load_aal.dat.scl_inter;
    
    dtype = 'FLOAT32';
    offset = 0;
    dim = load_aal.dat.dim;

    
    %%% POSITIVE
    
    descrip = 'TrackPassive';
    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-','Passive','-','Helinger','-','R',int2str(irun),'-','pos','.nii');
    input_data = TrackPassive_pos; 
    real_save_image;
    
    descrip = 'PassiveRestingState';
    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-','RestingState','-','Helinger','-','R',int2str(irun),'-','pos','.nii');
    input_data = PassiveRestingState_pos; 
    real_save_image;
    
    descrip = 'TrackRestingState';
    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-','RestingState','-','Helinger','-','R',int2str(irun),'-','pos','.nii');
    input_data = TrackRestingState_pos; 
    real_save_image;
    
    
    %%% NEGATIVE
    
    descrip = 'TrackPassive';
    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-','Passive','-','Helinger','-','R',int2str(irun),'-','neg','.nii');
    input_data = TrackPassive_neg; 
    real_save_image;
    
    descrip = 'PassiveRestingState';
    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-','RestingState','-','Helinger','-','R',int2str(irun),'-','neg','.nii');
    input_data = PassiveRestingState_neg; 
    real_save_image;
    
    descrip = 'TrackRestingState';
    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-','RestingState','-','Helinger','-','R',int2str(irun),'-','neg','.nii');
    input_data = TrackRestingState_neg; 
    real_save_image;
    
end

end

function probability_distribution = getProbabilityDistribution(rho_mat,pval_mat)

pcriterion = 0.01;
idx = find(pval_mat>pcriterion);
rho_mat(idx) = 0;

rho_mat_nozeros = rho_mat;
rho_mat_nozeros(rho_mat==0) = [];

step = 0.001;
interval = -1:step:1;

rho_mat_fisher = (0.5) * log( (1 + rho_mat_nozeros(:)) ./ (1 - rho_mat_nozeros(:)) );

k = ksdensity(rho_mat_fisher(:),interval);

probability_distribution = k .* step;

end

function [interval, probability_distribution] = getProbabilityDistributionZFisher(rho_mat,pval_mat)

pcriterion = 0.01;
idx = find(pval_mat>pcriterion);
rho_mat(idx) = 0;

[col,lin] = size(rho_mat);
identity = eye(col,lin);

rho_mat_nozeros = rho_mat;
rho_mat_nonzeros(logical(identity)) = 0;
rho_mat_nozeros(rho_mat==0) = [];

step = 0.001;
interval = -9:step:9;

rho_mat_fisher = (0.5) * log( (1 + rho_mat_nozeros(:)) ./ (1 - rho_mat_nozeros(:)) );

k = ksdensity(rho_mat_fisher(:),interval);

probability_distribution = k .* step;

end

function probability_distribution = getProbabilityDistributionGroup(rho_mat,pval_mat)

nRhos = length(rho_mat);

all_rhos = [];

for iRho=1:nRhos

    pcriterion = 0.01;
    idx = find(pval_mat(iRho).pval>pcriterion);
    rho_mat(iRho).rho(idx) = 0;

    rho_mat_nozeros = rho_mat(iRho).rho;
    rho_mat_nozeros(rho_mat(iRho).rho==0) = [];

    step = 0.001;
    interval = -1:step:1;

    rho_mat_fisher = (0.5) * log( (1 + rho_mat_nozeros(:)) ./ (1 - rho_mat_nozeros(:)) );

    all_rhos = [all_rhos(:);rho_mat_fisher(:)];
    
end

k = ksdensity(all_rhos(:),interval);

probability_distribution = k .* step;

end

function [interval, probability_distribution] = getProbabilityDistributionGroupZFisher(rho_mat,pval_mat)

nRhos = length(rho_mat);

all_rhos = [];

for iRho=1:nRhos

    pcriterion = 0.01;
    idx = find(pval_mat(iRho).pval>pcriterion);
    rho_mat(iRho).rho(idx) = 0;

    rho_mat_nozeros = rho_mat(iRho).rho;
    rho_mat_nozeros(rho_mat(iRho).rho==0) = [];

    step = 0.001;
    interval = -9:step:9;

    rho_mat_fisher = (0.5) * log( (1 + rho_mat_nozeros(:)) ./ (1 - rho_mat_nozeros(:)) );

    all_rhos = [all_rhos(:);rho_mat_fisher(:)];
    
end

k = ksdensity(all_rhos(:),interval);

probability_distribution = k .* step;

end

function distance = getHellingerDistance(dist1,dist2)

X = dist1;

Y = dist2;
        
XY = X.*Y;
sqrtXY = sqrt(XY);
        
distance = sqrt( 1 - sum( sqrtXY ) );
        
end 

function signal = getSignal(dist)

all_negative = find(dist<0);
all_positive = find(dist>0);

if sum(all_positive) > sum(all_negative)
    
    signal = 1;
    
else
    
    signal = -1;
    
end

end

function plotDistributions(settings,idx_ROI)

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRuns = 4;

%idx_ROI = 1:90;

step = 0.001;
interval = -1:step:1;

TrackPassive = zeros(length(idx_ROI),nRuns);
PassiveRestingState = zeros(length(idx_ROI),nRuns);
TrackRestingState = zeros(length(idx_ROI),nRuns);

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
        
        max_y = max([max(dist1) max(dist2)]);
        
        f = figure;
        
        plot(interval,dist1,'b');
        hold on
        plot(interval,dist2,'r');
        
        xlim([min(interval) max(interval)]);
        ylim([0 max_y]);
        
        legend({'Track','Passive'});
        
        print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-','Passive','-',area_label,'-','distribution','-','R',int2str(irun),'.eps'));
        print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-','Passive','-',area_label,'-','distribution','-','R',int2str(irun),'.jpg'));
        
        close all
        
        disp('Track-RestingState');
        dist1 = Track.TrackDistribution.run(irun).probability_distribution;
        dist2 = RestingState.RestingStateDistribution.run(irun).probability_distribution;
        
        max_y = max([max(dist1) max(dist2)]);
        
        f = figure;
        
        plot(interval,dist1,'b');
        hold on
        plot(interval,dist2,'r');
        
        xlim([min(interval) max(interval)]);
        ylim([0 max_y]);
        
        legend({'Track','RestingState'});
        
        print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-','RestingState','-',area_label,'-','distribution','-','R',int2str(irun),'.eps'));
        print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-','RestingState','-',area_label,'-','distribution','-','R',int2str(irun),'.jpg'));
        
        close all
        
        disp('Passive-RestingState');
        dist1 = Passive.PassiveDistribution.run(irun).probability_distribution;
        dist2 = RestingState.RestingStateDistribution.run(irun).probability_distribution;
        
        max_y = max([max(dist1) max(dist2)]);
        
        f = figure;
        
        plot(interval,dist1,'b');
        hold on
        plot(interval,dist2,'r');
        
        xlim([min(interval) max(interval)]);
        ylim([0 max_y]);
        
        legend({'Passive','RestingState'});
        
        print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-','RestingState','-',area_label,'-','distribution','-','R',int2str(irun),'.eps'));
        print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-','RestingState','-',area_label,'-','distribution','-','R',int2str(irun),'.jpg'));

        close all
        
     end
    
end

end

function compute3DHelingerDistanceAllRuns(settings)

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRuns = 4;

idx_ROI = 1:90;

TrackPassive = zeros(length(idx_ROI),nRuns);
PassiveRestingState = zeros(length(idx_ROI),nRuns);
TrackRestingState = zeros(length(idx_ROI),nRuns);

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
        [shift_pos_tp(iROI,irun), shift_neg_tp(iROI,irun)] = getShift(dist1,dist2);
        
        TrackPassive(iROI,irun) = distance;
        
        disp('Track-RestingState');
        dist1 = Track.TrackDistribution.run(irun).probability_distribution;
        dist2 = RestingState.RestingStateDistribution.run(irun).probability_distribution;
        distance = getHelingerDistance(dist1,dist2);
        [shift_pos_tr(iROI,irun), shift_neg_tr(iROI,irun)] = getShift(dist1,dist2);
        
        TrackRestingState(iROI,irun) = distance;
        
        disp('Passive-RestingState');
        dist1 = Passive.PassiveDistribution.run(irun).probability_distribution;
        dist2 = RestingState.RestingStateDistribution.run(irun).probability_distribution;
        distance = getHelingerDistance(dist1,dist2);
        [shift_pos_pr(iROI,irun), shift_neg_pr(iROI,irun)] = getShift(dist1,dist2);
        
        PassiveRestingState(iROI,irun) = distance;
        
    end
    
end

mn_shift_pos_tp = mean(shift_pos_tp,2);
mn_shift_pos_tr = mean(shift_pos_tr,2);
mn_shift_pos_pr = mean(shift_pos_pr,2);
mn_shift_neg_tp = mean(shift_neg_tp,2);
mn_shift_neg_tr = mean(shift_neg_tr,2);
mn_shift_neg_pr = mean(shift_neg_pr,2);
mn_TrackPassive = mean(TrackPassive,2);
mn_TrackRestingState = mean(TrackRestingState,2);
mn_PassiveRestingState = mean(PassiveRestingState,2);
    
TrackPassive_img = zeros(size(AAL_img));
PassiveRestingState_img = zeros(size(AAL_img));
TrackRestingState_img = zeros(size(AAL_img));

for iROI=1:length(idx_ROI)

    AAL_ID = AAL_ROI(idx_ROI(iROI)).ID;

    idx_voxels = find(AAL_img==AAL_ID);

    nVoxels = length(idx_voxels);

    for iVoxel=1:nVoxels

        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
        
        if mn_shift_pos_tp(iROI)>mn_shift_neg_tp(iROI)
            
            signal = 1;
            
            if mn_shift_pos_tp(iROI) > 0, increase = 1; else increase = 0; end
            
        else
            
            signal = -1;
            
            if mn_shift_neg_tp(iROI) > 0, increase = 1; else increase = 0; end
            
        end

        TrackPassive_img(idxx,idxy,idxz) = increase.*signal.*mn_TrackPassive(iROI);
        
       if mn_shift_pos_pr(iROI)>mn_shift_neg_pr(iROI)
            
            signal = 1;
            
            if mn_shift_pos_pr(iROI) > 0, increase = 1; else increase = 0; end
            
        else
            
            signal = -1;
            
            if mn_shift_neg_pr(iROI) > 0, increase = 1; else increase = 0; end
            
        end
        
        PassiveRestingState_img(idxx,idxy,idxz) = increase.*signal.*mn_PassiveRestingState(iROI);
        
        if mn_shift_pos_tr(iROI)>mn_shift_neg_tr(iROI)
            
            signal = 1;
            
            if mn_shift_pos_tr(iROI) > 0, increase = 1; else increase = 0; end
            
        else
            
            signal = -1;
            
            if mn_shift_neg_tr(iROI) > 0, increase = 1; else increase = 0; end
            
        end
        
        TrackRestingState_img(idxx,idxy,idxz) = increase.*signal.*mn_TrackRestingState(iROI);

    end

end

TrackPassive_pos = TrackPassive_img;
PassiveRestingState_pos = PassiveRestingState_img;
TrackRestingState_pos = TrackRestingState_img;

TrackPassive_neg = TrackPassive_img.*(-1);
PassiveRestingState_neg = PassiveRestingState_img.*(-1);
TrackRestingState_neg = TrackRestingState_img.*(-1);

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;
dim = load_aal.dat.dim;


%%% POSITIVE

descrip = 'TrackPassive';
fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-','Passive','-','Helinger','-','AllRuns','-','pos-shift','.nii');
input_data = TrackPassive_pos; 
real_save_image;

descrip = 'PassiveRestingState';
fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-','RestingState','-','Helinger','-','AllRuns','-','pos-shift','.nii');
input_data = PassiveRestingState_pos; 
real_save_image;

descrip = 'TrackRestingState';
fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-','RestingState','-','Helinger','-','AllRuns','-','pos-shift','.nii');
input_data = TrackRestingState_pos; 
real_save_image;


%%% NEGATIVE

descrip = 'TrackPassive';
fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-','Passive','-','Helinger','-','AllRuns','-','neg-shift','.nii');
input_data = TrackPassive_neg; 
real_save_image;

descrip = 'PassiveRestingState';
fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-','RestingState','-','Helinger','-','AllRuns','-','neg-shift','.nii');
input_data = PassiveRestingState_neg; 
real_save_image;

descrip = 'TrackRestingState';
fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-','RestingState','-','Helinger','-','AllRuns','-','neg-shift','.nii');
input_data = TrackRestingState_neg; 
real_save_image;
       
end

function [shift_pos, shift_neg] = getShift(dist1,dist2)

step = 0.001;
interval = -1:step:1;

mid = round(length(interval)/2);

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

end

function shift = getShiftSignificant(rho1,pval1,rho2,pval2)

NComponent = size(rho1,1);

pcriterion = 0.01;

idx1 = find(pval1>pcriterion);
rho1(idx1) = 0;
rho1(rho1==0) = [];

idx2 = find(pval2>pcriterion);
rho2(idx2) = 0;
rho2(rho2==0) = [];

percent1 = length(rho1)/(NComponent*NComponent);
percent2 = length(rho2)/(NComponent*NComponent);

shift = percent1 - percent2;

end