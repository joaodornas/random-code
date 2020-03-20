
%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_nodes_frontal = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];

for iclean=1:5
    
    f = figure;
    for iNode=1:length(idx_nodes_frontal)

        subplot(5,4,iNode);

        idx_ROI = AAL_ROI(idx_nodes_frontal(iNode)).ID;
        idx_ROI_Label = AAL_ROI(idx_nodes_frontal(iNode)).Nom_L;
        idx_ROI_Label = strrep(idx_ROI_Label,'_','-');

        idx_voxels = find(AAL_img == idx_ROI);

        cross_idx = ismember(idx_brain,idx_voxels);

        nVoxels = length(find(cross_idx));

        area = zeros(nVoxels,nTR);

        area(1:nVoxels,1:nTR) = condition(cross_idx,1:nTR);

        mean_area = mean(area,1);

        area_clean(1:nVoxels,1:nTR) = clean(iclean).ts(cross_idx,1:nTR);

        mean_area_clean = mean(area_clean,1);

        plot(mean_area,'k');
        hold on
        plot(mean_area_clean,'b');

        xlim([0 nTR]);
        xlabel('TR');
        ylabel('BOLD');
        title(strcat('Threshold-',int2str(iclean),'-',idx_ROI_Label));

        clear area

    end

end
