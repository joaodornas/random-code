
%% LOAD SETTINGS 

settings_jan_0805;
get_at_this_preprocessed_step = 'mcflirt';
settings = lowhigh_concatenate_folders_strings_FSL(settings,get_at_this_preprocessed_step);

run = 2;

%% LOAD DISPLACEMENTS

displaces = load(strcat(settings.mot4.FSL.run(run).fullpath,'\','displacements.mat'));
displacements = displaces.displacements;

firstTR = 1;
lastTR = 331;
nTR = lastTR - firstTR + 1;

prefix_file = 'MOT4-Run2_mcf_';

% costlabel{1} = 'corratio';
% costlabel{2} = 'leastsquares';
% costlabel{3} = 'mutualinfo';
% costlabel{4} = 'normcorr';
% costlabel{5} = 'normmi';
% costlabel{6} = 'woods';

costlabel{1} = 'woods_spline';
costlabel{2} = 'woods_spline_4';
costlabel{3} = 'woods_sinc';

%% LOAD MOTION CORRECTED FUNCTIONS

for iCost=1:length(costlabel)

    file = strcat(prefix_file,costlabel{iCost},'2standard');
    filename = strcat(settings.mot4.FSL.run(run).fullpath,'\',file,'.nii');
    nifti_file = nifti(filename);
    dat = nifti_file.dat;
    BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
    output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),1:nTR) = BOLD(:,:,:,firstTR:lastTR).*dat.scl_slope + dat.scl_inter;
    
    costcorrected(iCost).MOT4Run2 = output_data;
    
    clear output_data

end

%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

grand_area = 'Frontal';
disp(grand_area);
idx_nodes_frontal = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];

idx_nodes = idx_nodes_frontal;

%% GET AND PLOT HIGHEST MEAN VOXEL

for iCost=1:length(costlabel)
    
    output_data = costcorrected(iCost).MOT4Run2;
    
    mc_abs = displacements(1:end,iCost);

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

    f = figure;
    
    for j=1:length(idx_nodes)

        subplot(5,4,j);
        
        plot(zscore(hm_voxel(j,:)./max(hm_voxel(j,:))),'b');
        hold on
        plot(zscore(mc_abs./max(mc_abs)),'r');

        xlabel('TRs');
        %legend({'Highest Mean Voxel (BOLD)', 'Absolute Mean Displacement'});

        xlim([0 nTR]);

        title(strcat(labels_AAL{j},'-',costlabel{iCost}));

        print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',grand_area,'-HM-voxel-MCFlirt-',costlabel{iCost},'.jpeg'));
        print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',grand_area,'-HM-voxel-MCFlirt-',costlabel{iCost},'.eps'));
        print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',grand_area,'-HM-voxel-MCFlirt-',costlabel{iCost},'.pdf'));
        
    end
    
    close all
     
end

