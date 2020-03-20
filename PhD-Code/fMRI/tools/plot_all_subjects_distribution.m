
function plot_all_subjects_distribution(P_OnlyTrack_sig)

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

sig_ROI = find(P_OnlyTrack_sig);

nSigROI = length(sig_ROI);

step = 0.001;
interval = -1:step:1;

for iROI=1:nSigROI
    
    label_ROI = AAL_ROI(sig_ROI(iROI)).Nom_L;

    area_label = strrep(label_ROI,'_','-');       

    disp(area_label);

    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-','distribution','.mat'));
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-','distribution','.mat'));
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-','distribution','.mat'));

    f = figure;
    
    plot(interval,TrackDistribution,'b');
    hold on
    plot(interval,PassiveDistribution,'r');
    plot(interval,RestingStateDistribution,'k');
    
    xlim([min(interval) max(interval)]);
    ylim([0 max([TrackDistribution(:);PassiveDistribution(:);RestingStateDistribution(:)])]);
    
    legend({'Track' 'Passive' 'RestingState'});
    title(area_label);
    
    print(f,'-djpeg',strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','T-P-R','-',area_label,'-','distribution','.jpeg'));
    print(f,'-depsc',strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','T-P-R','-',area_label,'-','distribution','.eps'));
    
end

end