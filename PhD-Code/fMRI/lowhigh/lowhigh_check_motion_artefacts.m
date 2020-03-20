

function lowhigh_check_motion_artefacts

%% SETTINGS

settings_jan_0805;
%settings_elena_2905;

get_at_this_preprocessed_step = settings.FSL.folders.warped;
settings = lowhigh_concatenate_folders_strings_FSL(settings,get_at_this_preprocessed_step);

file{1} = settings.FSL.files.functional.warped;
file{2} = settings.FSL.files.functional.despike;
%file{3} = settings.FSL.files.functional.regressed;

filelabel{1} = 'warped';
filelabel{2} = 'despike';
%filelabel{3} = 'regressed';

Run{1} = 'MOT4-Run1';
Run{2} = 'MOT4-Run2';
Run{3} = 'MOT2-Run1';
Run{4} = 'MOT2-Run2';
Run{5} = 'RestingState-Run1';
Run{6} = 'RestingState-Run2';

Runkind{1} = 'MOT4';
Runkind{2} = 'MOT4';
Runkind{3} = 'MOT2';
Runkind{4} = 'MOT2';
Runkind{5} = 'RestingState';
Runkind{6} = 'RestingState';

Runidx(1) = 1;
Runidx(2) = 2;
Runidx(3) = 1;
Runidx(4) = 2;
Runidx(5) = 1;
Runidx(6) = 2;

nRun = 6;

GM_mask = getGMMask(settings);

%% AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_nodes_frontal = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];
idx_nodes_temporal = [81 82 83 84 85 86 87 88 89 90];

PFROI = 4;
TEMPROI = 3;

for iRun=1:nRun
%for iRun=2:2
    
    [output_data_warped, mask_data_warped, settings] = lowhigh_get_data_FSL(settings,Runkind{iRun},Runidx(iRun),file{1});
    [output_data_despike, mask_data_despike, settings] = lowhigh_get_data_FSL(settings,Runkind{iRun},Runidx(iRun),file{2});
    %[output_data_regressed, mask_data_regressed, settings] = lowhigh_get_data_FSL(settings,Runkind{iRun},Runidx(iRun),file{3});

    DVARS_warped = getDVARS(output_data_warped,GM_mask);
    DVARS_despike = getDVARS(output_data_despike,GM_mask);
    %DVARS_regressed = getDVARS(output_data_regressed,GM_mask);
    
    [x,y,z,FD] = getFD(settings,Runkind{iRun},Runidx(iRun));
    
    fs = 6;
    
    f = figure;
    
    idx = idx_nodes_frontal(PFROI);
    
    idx_AAL = AAL_ROI(idx).ID;
    
    label_AAL = AAL_ROI(idx).Nom_L;
    label_AAL = strrep(label_AAL,'_','-');
    
    idx_voxels = find(AAL_img==idx_AAL);
    
    nVoxels = length(idx_voxels);
    nTR = size(output_data_warped,4);
    
    area_warped = zeros(nVoxels,nTR);
    area_despike = zeros(nVoxels,nTR);
    %area_regressed = zeros(nVoxels,nTR);
    
    for iVoxel=1:nVoxels
        
       [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
       
       area_warped(iVoxel,:) = output_data_warped(idxx,idxy,idxz,:);
       area_despike(iVoxel,:) = output_data_despike(idxx,idxy,idxz,:);
       %area_regressed(iVoxel,:) = output_data_regressed(idxx,idxy,idxz,:);
        
    end
    
    mean_area_warped = zscore(mean(area_warped,1));
    mean_area_despike = zscore(mean(area_despike,1));
    %mean_area_regressed = zscore(mean(area_regressed,1));
 
    subplot(3,2,1);
    
    plot(mean_area_warped,'b');
    hold on
    plot(mean_area_despike,'r');
    %plot(mean_area_regressed,'y');
    
    xlim([0 nTR]);
    ylabel('BOLD - mean, zscored','FontSize',fs);
    xlabel('TRs','FontSize',fs);
    legend({'warped','wavelet'});
    
    title(strcat(Run{iRun},'-',label_AAL));
    
    mean_area_warped = mean(area_warped,1);
    mean_area_despike = mean(area_despike,1);
    %mean_area_regressed = mean(area_regressed,1);
    
    subplot(3,2,2);
    
    plot(mean_area_warped,'b');
    hold on
    plot(mean_area_despike,'r');
    %plot(mean_area_regressed,'y');
    
    xlim([0 nTR]);
    ylabel('BOLD','FontSize',fs);
    xlabel('TRs','FontSize',fs);
    legend({'warped','wavelet'});
    
    title(strcat(Run{iRun},'-',label_AAL));
    
    clear area_warped
    clear area_despike
    clear area_regressed
    
    idx = idx_nodes_temporal(TEMPROI);
    
    idx_AAL = AAL_ROI(idx).ID;
    
    label_AAL = AAL_ROI(idx).Nom_L;
    label_AAL = strrep(label_AAL,'_','-');
    
    idx_voxels = find(AAL_img==idx_AAL);
    
    nVoxels = length(idx_voxels);
    nTR = size(output_data_warped,4);
    
    area_warped = zeros(nVoxels,nTR);
    area_despike = zeros(nVoxels,nTR);
    %area_regressed = zeros(nVoxels,nTR);
    
    for iVoxel=1:nVoxels
        
       [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
       
       area_warped(iVoxel,:) = output_data_warped(idxx,idxy,idxz,:);
       area_despike(iVoxel,:) = output_data_despike(idxx,idxy,idxz,:);
       %area_regressed(iVoxel,:) = output_data_regressed(idxx,idxy,idxz,:);
        
    end
    
    mean_area_warped = zscore(mean(area_warped,1));
    mean_area_despike = zscore(mean(area_despike,1));
    %mean_area_regressed = zscore(mean(area_regressed,1));
    
    subplot(3,2,3);
    
    plot(mean_area_warped,'b');
    hold on
    plot(mean_area_despike,'r');
    %plot(mean_area_regressed,'y');
    
    xlim([0 nTR]);
    ylabel('BOLD - mean, zscored','FontSize',fs);
    xlabel('TRs','FontSize',fs);
    legend({'warped','wavelet'});
    
    title(strcat(Run{iRun},'-',label_AAL));
    
    mean_area_warped = mean(area_warped,1);
    mean_area_despike = mean(area_despike,1);
    %mean_area_regressed = mean(area_regressed,1);
    
    subplot(3,2,4);
    
    plot(mean_area_warped,'b');
    hold on
    plot(mean_area_despike,'r');
    %plot(mean_area_regressed,'y');
    
    xlim([0 nTR]);
    ylabel('BOLD','FontSize',fs);
    xlabel('TRs','FontSize',fs);
    legend({'warped','wavelet'});
    
    title(strcat(Run{iRun},'-',label_AAL));
    
    clear area_warped
    clear area_despike
    clear area_regressed
    
    subplot(3,2,5);
    
    plot(x,'b');
    hold on
    plot(y,'g');
    plot(z,'k');
    legend({'x','y','z'});
    
    xlim([0 nTR]);
    ylabel('mm','FontSize',fs);
    xlabel('TRs','FontSize',fs);
    
    title('Displacement');
    
    subplot(3,2,6);
    
    plot(DVARS_warped,'b');
    hold on
    plot(DVARS_despike,'r');
    %plot(DVARS_regressed,'y');
    legend({'warped','wavelet'});
    
    title('DVARS');
    
    xlim([0 nTR]);
    ylabel('\Delta BOLD','FontSize',fs);
    xlabel('TRs','FontSize',fs);
    
%     subplot(4,1,4);
%     
%     plot(FD,'r');
%     
%     xlim([0 nTR]);
%     ylabel('mm','FontSize',fs);
%     xlabel('TRs','FontSize',fs);
    
    clear output_data_warped
    clear output_data_despike
    clear output_data_regressed
    
    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Run{iRun},'-','AAL-time-series-Displacement-DVARS','.jpg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Run{iRun},'-','AAL-time-series-Displacement-DVARS','.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Run{iRun},'-','AAL-time-series-Displacement-DVARS','.pdf'));
    
    close all;
    
end


end

function DVARS = getDVARS(data,mask)

    nTR = size(data,4);
    
    idx_mask = find(mask);
    
    nVoxels = length(idx_mask);
    
    derivative = zeros(nVoxels,nTR);
    
    for iVoxel=1:nVoxels
            
        [idxx, idxy, idxz] = ind2sub(size(mask),idx_mask(iVoxel));
            
        time_series = data(idxx,idxy,idxz,:);
            
        diff_time_series = diff(time_series);
            
        derivative(iVoxel,2:nTR) = diff_time_series;
        
    end
    
    for iTR=1:nTR
        
        voxels = derivative(1:end,iTR);
        
        sq_voxels = voxels.^2;
        
        mean_voxels = mean(sq_voxels);
        
        sqrt_voxels = sqrt(mean_voxels);
        
        DVARS(iTR) = sqrt_voxels;
        
    end

end

function mask = getGMMask(settings)

    folder = settings.folders.main;
    folder = strcat(folder,'\',settings.folders.experiment);
    folder = strcat(folder,'\',settings.folders.subject);
    folder = strcat(folder,'\',settings.folders.preprocessed);
    FSL = settings.FSL.folders.main;
    
    tmp = strcat(folder,'\',settings.anatomical.folder.main,'\',settings.anatomical.run(1).folder,'\',FSL);
    
    mask_folder = strcat(tmp,'\','segmented.normalized');
    
    file_full_path = strcat(mask_folder,'\',settings.anatomical.run(1).segmentation.mask.GM);
    
    maskdata = nifti(file_full_path);
    
    maskdata.dat.fname = file_full_path;
    
    mask = maskdata.dat(:,:,:);
    
end

function [x,y,z,FD] = getFD(settings,kind,run)

    folder = settings.folders.main;
    folder = strcat(folder,'\',settings.folders.experiment);
    folder = strcat(folder,'\',settings.folders.subject);
    folder = strcat(folder,'\',settings.folders.preprocessed);
    FSL = settings.FSL.folders.main;
    
    all_kind = {'MOT4', 'MOT2', 'RestingState'};
    idx_kind = find(strcmp(kind,all_kind));

switch idx_kind
    
    case 1
        
        tmp = strcat(folder,'\',settings.functional.mot4.folder.main);
        mc_folder = strcat(tmp,'\',settings.functional.mot4.run(run).folder,'\',FSL,'\','Melodic-Fieldmap.ica','\','mc','\','prefiltered_func_data_mcf.mat');
   
    case 2
        
        tmp = strcat(folder,'\',settings.functional.mot2.folder.main);
        mc_folder = strcat(tmp,'\',settings.functional.mot2.run(run).folder,'\',FSL,'\','Melodic-Fieldmap.ica','\','mc','\','prefiltered_func_data_mcf.mat');
        
    case 3
        
        tmp = strcat(folder,'\',settings.functional.restingstate.folder.main);
        mc_folder = strcat(tmp,'\',settings.functional.restingstate.run(run).folder,'\',FSL,'\','Melodic-Fieldmap.ica','\','mc','\','prefiltered_func_data_mcf.mat');
   
end

matrix = load(strcat(mc_folder,'\','mat.mat'));

mc = matrix.mat;

nMAT = length(mc);

for iMAT=1:nMAT
   
    x(iMAT) = mc(iMAT).mc(1,4);
    y(iMAT) = mc(iMAT).mc(2,4);
    z(iMAT) = mc(iMAT).mc(3,4);
    
end

FD = 0;

end
