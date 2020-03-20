
%% LOAD SETTINGS

settings_jan_0805;

run = 2;

%%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.warped;
file = settings.FSL.files.functional.warped;

settings = lowhigh_concatenate_folders_strings_FSL(settings,get_at_this_preprocessed_step);
     
firstTR = 1;
lastTR = 331;

nTR = lastTR - firstTR + 1;

filename = strcat(settings.mot4.FSL.run(run).fullpath,'\',file,'.nii');

nifti_file = nifti(filename);
dat = nifti_file.dat;
BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),1:nTR) = BOLD(:,:,:,firstTR:lastTR).*dat.scl_slope + dat.scl_inter;
        
%% LOAD MCFLIRT

mc = load('MOT4Run2_mc_abs.mat');
mc_abs = mc.MOT4Run2_mc_abs;

%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

disp('Frontal');
idx_nodes_frontal = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];

disp('Occipital');
idx_nodes_occipital = [49 50 51 52 53 54];

disp('Parietal');
idx_nodes_parietal = [59 60 61 62];

disp('Temporal');
idx_nodes_temporal = [81 82 83 84 85 86 87 88 89 90];

idx_nodes = [idx_nodes_frontal, idx_nodes_occipital, idx_nodes_parietal, idx_nodes_temporal];

%% GET HIGHEST MEAN VOXEL

for j=1:length(idx_nodes)

    idx = idx_nodes(j);
    
    idx_AAL = AAL_ROI(idx).ID;
    idx_voxels = find(AAL_img==idx_AAL); 
    
    label = AAL_ROI(idx).Nom_L;
   
    labels_AAL{j} = strrep(label,'_','-');
    
    area = zeros(length(idx_voxels),nTR);
    
    for iVoxel=1:length(idx_voxels)
        
        [idxx, idxy, idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
        
        area(iVoxel,:) = output_data(idxx,idxy,idxz,:);
        
    end
    
    mean_area = mean(area,2);
    
    [x, I] = sort(mean_area);
    
    idx_max_mean = I(end);
    
    hm_voxel(j,:) = area(idx_max_mean,:);

    clear area
    
end

%% PLOT HIGHEST MEAN + MOTION CORRECTION

for j=1:length(idx_nodes)
        
    f = figure;

    plot(zscore(hm_voxel(j,:)./max(hm_voxel(j,:))),'b');
    hold on
    plot(zscore(mc_abs./max(mc_abs)),'r');

    xlabel('TRs');
    legend({'Highest Mean Voxel (BOLD)', 'Absolute Mean Displacement'});
    
    xlim([0 nTR]);
    
    title(labels_AAL{j});
    
    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',labels_AAL{j},'-HM-voxel-MCFlirt','.jpeg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',labels_AAL{j},'-HM-voxel-MCFlirt','.eps'));
    
    close all
    
end
    

