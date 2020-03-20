
function lowhigh_plot_PCs

%%% SETTINGS

settings_jan_0805; %% MOT4Run2
lowhigh_load_all_data_FSL;
RUN = MOT4Run2;
RUN_label = 'MOT4Run2';

doTheMath(settings,RUN,RUN_label);

%settings_elena_2905; %% RestingStateRun1
%RUN = RestingStateRun1;
%RUN_label = 'RestingStateRun1';

%doTheMath(settings,RUN,RUN_label);

return

function doTheMath(settings,RUN,RUN_label)

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

nTR = size(RUN,4);

disp('Frontal');
idx_nodes_frontal = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];

disp('Occipital');
idx_nodes_occipital = [49 50 51 52 53 54];

disp('Parietal');
idx_nodes_parietal = [59 60 61 62];

disp('Temporal');
idx_nodes_temporal = [81 82 83 84 85 86 87 88 89 90];

idx_nodes = [idx_nodes_frontal, idx_nodes_occipital, idx_nodes_parietal, idx_nodes_temporal];

%for iNode=1:length(idx_nodes)
for iNode=1:1

    AAL_ROI_ID = AAL_ROI(idx_nodes(iNode)).ID;
    AAL_ROI_Label = AAL_ROI(idx_nodes(iNode)).Nom_L;
    AAL_ROI_Label = strrep(AAL_ROI_Label,'_','-');

    idx_voxels = find(AAL_img == AAL_ROI_ID);
    nVoxels = length(idx_voxels);

    disp(AAL_ROI_Label);
    disp(strcat('nVoxels:',int2str(nVoxels)));

    area_RUN = zeros(nVoxels,nTR);

    for iVoxel=1:nVoxels
        
        idx = idx_voxels(iVoxel);

        [idxx(iVoxel), idxy(iVoxel), idxz(iVoxel)] = ind2sub(size(AAL_img),idx);

        voxel = squeeze(RUN(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:));
        
        area_RUN(iVoxel,:) = voxel;
  
    end

    mn_RUN = mean( area_RUN, 2);   
    X0_area_RUN = area_RUN - repmat(mn_RUN,1,nTR);

    [U_area_RUN, S_area_RUN, V_area_RUN] = svd( X0_area_RUN' ); 

    V_area_RUN = V_area_RUN';
    Y0_area_RUN = V_area_RUN*X0_area_RUN;

    %% GET BY AMPLITUDE
%     max_all = max(max(Y0_area_RUN));
%     idx_max = find(Y0_area_RUN == max_all);
%     [idxx, idxy] = ind2sub(size(Y0_area_RUN),idx_max);
    
    disp(strcat('PC:',int2str(idxx)));

    Y0dn_area_RUN = Y0_area_RUN;    
    Y0dn_area_RUN(idxx,:) = zeros(1,nTR);
    
    X0dn_area_RUN = V_area_RUN'*Y0dn_area_RUN;
    Xdn_area_RUN = X0dn_area_RUN + repmat(mn_RUN,1,nTR);
    
    mean_RUN = mean(area_RUN,1);
    meandn_RUN = mean(Xdn_area_RUN,1);
    
    f = figure;
    
    plot(mean_RUN,'b');
    hold on
    plot(meandn_RUN,'r');
    
    legend({'Mean','Mean denoised'});
    xlabel('TR');
    ylabel('BOLD Activity');
  
    min_y = min([min(mean_RUN),min(meandn_RUN)]);
    max_y = max([max(mean_RUN),max(meandn_RUN)]);
    
    ylim([min_y max_y]);
    
    fs = 14;
    
    x = nTR/2;
    y = min_y + (max_y - min_y)/2;
    text(x,y,strcat('PC:',int2str(idxx)),'FontSize',fs);
    
    title(AAL_ROI_Label);
    
    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',RUN_label,'-',AAL_ROI_Label,'.jpg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',RUN_label,'-',AAL_ROI_Label,'.eps'));
    
    close all
    
end

return
    

