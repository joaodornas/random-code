function real_FC_voxel_AAL_ROI_Distributions_plot(settings,idx_ROI)

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRuns = 4;

for iROI=1:length(idx_ROI)

    label_ROI = AAL_ROI(idx_ROI(iROI)).Nom_L;

    area_label = strrep(label_ROI,'_','-');       

    disp(area_label);

    TrackDist = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-','distribution','.mat'));
    PassiveDist = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-','distribution','.mat'));
    RestingStateDist = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-','distribution','.mat'));

    step = 0.001;
    interval = -1:step:1;

    f = figure;

    for irun=1:nRuns

       plot(interval,TrackDist.TrackDistribution.run(irun).probability_distribution,'b');
       hold on
       plot(interval,PassiveDist.PassiveDistribution.run(irun).probability_distribution,'r');
       
       max_track(irun) = max(TrackDist.TrackDistribution.run(irun).probability_distribution);
       max_passive(irun) = max(PassiveDist.PassiveDistribution.run(irun).probability_distribution);
       
    end

    xlim([min(interval) max(interval)]);
    ylim([0 max([max_track(:);max_passive(:)])]);

    title(area_label);

    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Distributions','-',area_label,'-','T-P','.jpeg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Distributions','-',area_label,'-','T-P','.eps'));

    f = figure;

    for irun=1:nRuns

       plot(interval,PassiveDist.PassiveDistribution.run(irun).probability_distribution,'r');
       hold on
       plot(interval,RestingStateDist.RestingStateDistribution.run(irun).probability_distribution,'k');

       max_passive(irun) = max(PassiveDist.PassiveDistribution.run(irun).probability_distribution);
       max_rest(irun) = max(RestingStateDist.RestingStateDistribution.run(irun).probability_distribution);

    end

    xlim([min(interval) max(interval)]);
    ylim([0 max([max_passive(:);max_rest(:)])]);
    
    title(area_label);

    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Distributions','-',area_label,'-','P-R','.jpeg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Distributions','-',area_label,'-','P-R','.eps'));

end

end

