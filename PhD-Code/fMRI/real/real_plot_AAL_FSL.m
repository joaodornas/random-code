
function real_plot_AAL_FSL

%%% SETTINGS

%settings_subj1_2210;
%settings_subj2_2610;
%settings_subj3_0311;
%settings_subj4_0211;
%settings_subj5_0211;
settings_subj6_2411;

%plotAllROI(settings);

idx_ROI = 3; %% Frontal_Sup_L
plotWarpedWaveletResidualOneROI(settings,idx_ROI);

end

function plotAllROI(settings)

%%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.filtered;
mask = settings.FSL.files.mask.custom;

real_load_all_data_FSL;

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

meanTrack = getConditionMeanROI(Track,AAL_img,AAL_ROI);
meanPassive = getConditionMeanROI(Passive,AAL_img,AAL_ROI);
meanTrials = getConditionMeanROI(Trials,AAL_img,AAL_ROI);
meanRestingState = getConditionMeanROI(RestingState,AAL_img,AAL_ROI);

plotCondition(settings,meanTrack,AAL_ROI,'Track');
plotCondition(settings,meanPassive,AAL_ROI,'Passive');
plotCondition(settings,meanTrials,AAL_ROI,'Trials');
plotCondition(settings,meanRestingState,AAL_ROI,'RestingState');
   
end

function meanCondition = getConditionMeanROI(Condition,AAL_img,AAL_ROI)

nNodes = 90;
for iNode=1:nNodes
    
   idx_region = AAL_ROI(iNode).ID;
   idx_voxels_structures = find(AAL_img == idx_region);
   nVoxels = length(idx_voxels_structures);
   
   idxx = zeros(1,nVoxels);
   idxy = zeros(1,nVoxels);
   idxz = zeros(1,nVoxels);
   
   for iVoxel=1:nVoxels
        
        [idxx(iVoxel), idxy(iVoxel), idxz(iVoxel)] = ind2sub(size(AAL_img),idx_voxels_structures(iVoxel));
    
   end  
   
   nTR = size(Condition(1).run,4);
    
   Run1_voxel_time_series = zeros(nVoxels,nTR);
   Run2_voxel_time_series = zeros(nVoxels,nTR);
   Run3_voxel_time_series = zeros(nVoxels,nTR);
   Run4_voxel_time_series = zeros(nVoxels,nTR);
    
   for iVoxel=1:nVoxels
        
        Run1_voxel_time_series(iVoxel,1:nTR) = Condition(1).run(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        Run2_voxel_time_series(iVoxel,1:nTR) = Condition(2).run(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        Run3_voxel_time_series(iVoxel,1:nTR) = Condition(3).run(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        Run4_voxel_time_series(iVoxel,1:nTR) = Condition(4).run(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
    
   end
    
   Run1_AAL = mean(Run1_voxel_time_series,1);
   Run2_AAL = mean(Run2_voxel_time_series,1);
   Run3_AAL = mean(Run3_voxel_time_series,1);
   Run4_AAL = mean(Run4_voxel_time_series,1);

   meanCondition(1).run(iNode,:) = Run1_AAL(:);
   meanCondition(2).run(iNode,:) = Run2_AAL(:);
   meanCondition(3).run(iNode,:) = Run3_AAL(:);
   meanCondition(4).run(iNode,:) = Run4_AAL(:);
    
end

end

function plotCondition(settings,Condition,AAL_ROI,Condition_Label)
        
    
colors{1} = 'b';
colors{2} = 'r';
colors{3} = 'k';
colors{4} = 'y';

fs = 6;
    
f = figure;

for iROI=1:30
    
    subplot(8,4,iROI);

    label = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    
    for irun=1:4
        
        plot(Condition(irun).run(iROI,:),colors{irun});
        hold on
        
    end
    
   title(label,'FontSize',fs);
   
%    xlabel('TR','FontSize',fs);
%    ylabel('BOLD Activity','FontSize',fs);
   
   set(gca,'FontSize',fs);
   
   if iROI == 30
       
       legend('Run 1','Run 2','Run 3','Run 4');
       
   end
   
   print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','AAL-Mean-1-30-',Condition_Label,'.jpg'));
   print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','AAL-Mean-1-30-',Condition_Label,'.eps'));
        
end

close all

f = figure;

for iROI=31:60
    
    subplot(8,4,iROI-30);
    
    label = strrep(AAL_ROI(iROI).Nom_L,'_','-');

    for irun=1:4
        
        plot(Condition(irun).run(iROI,:),colors{irun});
        hold on
        
    end
    
   title(label,'FontSize',fs);
   
%    xlabel('TR','FontSize',fs);
%    ylabel('BOLD Activity','FontSize',fs);
   
   set(gca,'FontSize',fs);
   
   if iROI == 60
       
       legend('Run 1','Run 2','Run 3','Run 4');
       
   end
   
   print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','AAL-Mean-31-60-',Condition_Label,'.jpg'));
   print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','AAL-Mean-31-60-',Condition_Label,'.eps'));
        
end

close all

f = figure;

for iROI=61:90
    
    subplot(8,4,iROI-60);
    
    label = strrep(AAL_ROI(iROI).Nom_L,'_','-');

    for irun=1:4
        
        plot(Condition(irun).run(iROI,:),colors{irun});
        hold on
        
    end
    
   title(label,'FontSize',fs);
   
%    xlabel('TR','FontSize',fs);
%    ylabel('BOLD Activity','FontSize',fs);
   
   set(gca,'FontSize',fs);
   
   if iROI == 90
       
       legend('Run 1','Run 2','Run 3','Run 4');
       
   end
   
   print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','AAL-Mean-61-90-',Condition_Label,'.jpg'));
   print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','AAL-Mean-61-90-',Condition_Label,'.eps'));
        
end

close all

end

function plotWarpedWaveletResidualOneROI(settings,idx_ROI)

analysis_label = 'WarpedWaveletResidual';

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

area_label = strrep(AAL_ROI(idx_ROI).Nom_L,'_','-');

%%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.filtered;
mask = settings.FSL.files.mask.custom;
real_load_all_data_FSL;

[Track_Filtered,Passive_Filtered,Trials_Filtered,RestingState_Filtered] = getOneROIMeanAllConditions(Track,Passive,Trials,RestingState,AAL_ROI,AAL_img,idx_ROI);

clear Track Passive Trials RestingState

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.wavelet_voxel;
mask = settings.FSL.files.mask.custom;
real_load_all_data_FSL;

[Track_Wavelet,Passive_Wavelet,Trials_Wavelet,RestingState_Wavelet] = getOneROIMeanAllConditions(Track,Passive,Trials,RestingState,AAL_ROI,AAL_img,idx_ROI);

clear Track Passive Trials RestingState

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;
real_load_all_data_FSL;

[Track_Residual,Passive_Residual,Trials_Residual,RestingState_Residual] = getOneROIMeanAllConditions(Track,Passive,Trials,RestingState,AAL_ROI,AAL_img,idx_ROI);

clear Track Passive Trials RestingState

nTR = length(Track_Filtered(1).ROI);

f = figure;
fs = 5;

for irun=1:4
    
    subplot(4,4,irun)

    plot(Track_Filtered(irun).ROI,'b');
    hold on
    plot(Track_Wavelet(irun).ROI,'r');
    plot(Track_Residual(irun).ROI,'k');

    xlim([0 nTR]);
    title(strcat('Track:R',int2str(irun),'-',area_label),'FontSize',fs);
    %xlabel('TR','FontSize',fs);
    %ylabel('BOLD Activity','FontSize',fs);
    set(gca,'FontSize',fs);

end

for irun=1:4
    
    subplot(4,4,irun+4)

    plot(Passive_Filtered(irun).ROI,'b');
    hold on
    plot(Passive_Wavelet(irun).ROI,'r');
    plot(Passive_Residual(irun).ROI,'k');

    xlim([0 nTR]);
    title(strcat('Passive:R',int2str(irun),'-',area_label),'FontSize',fs);
    %xlabel('TR','FontSize',fs);
    %ylabel('BOLD Activity','FontSize',fs);
    set(gca,'FontSize',fs);

end

for irun=1:4
    
    subplot(4,4,irun+8)

    plot(Trials_Filtered(irun).ROI,'b');
    hold on
    plot(Trials_Wavelet(irun).ROI,'r');
    plot(Trials_Residual(irun).ROI,'k');

    xlim([0 nTR]);
    title(strcat('Trials:R',int2str(irun),'-',area_label),'FontSize',fs);
    %xlabel('TR','FontSize',fs);
    %ylabel('BOLD Activity','FontSize',fs);
    set(gca,'FontSize',fs);

end

for irun=1:4
    
    subplot(4,4,irun+12)

    plot(RestingState_Filtered(irun).ROI,'b');
    hold on
    plot(RestingState_Wavelet(irun).ROI,'r');
    plot(RestingState_Residual(irun).ROI,'k');

    xlim([0 nTR]);
    title(strcat('RestingState:R',int2str(irun),'-',area_label),'FontSize',fs);
    %xlabel('TR','FontSize',fs);
    %ylabel('BOLD Activity','FontSize',fs);
    set(gca,'FontSize',fs);

end

legend({'Filtered' 'Wavelet' 'Residual'});

print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-',area_label,'.jpeg'));
print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-',area_label,'.eps'));

end

function [Track_Preprocess,Passive_Preprocess,Trials_Preprocess,RestingState_Preprocess] = getOneROIMeanAllConditions(Track,Passive,Trials,RestingState,AAL_ROI,AAL_img,idx_ROI)

idx_AAL = AAL_ROI(idx_ROI).ID;

idx_voxels = find(AAL_img == idx_AAL);

Track_Preprocess = getMeanROIAllRunsOneCondition(Track,idx_voxels,AAL_img);
Passive_Preprocess = getMeanROIAllRunsOneCondition(Passive,idx_voxels,AAL_img);
Trials_Preprocess = getMeanROIAllRunsOneCondition(Trials,idx_voxels,AAL_img);
RestingState_Preprocess = getMeanROIAllRunsOneCondition(RestingState,idx_voxels,AAL_img);

end

function Runs = getMeanROIAllRunsOneCondition(Condition,idx_voxels,AAL_img)

nVoxels = length(idx_voxels);
nTR = size(Condition(1).run,4);
nRun = length(Condition);

Runs(1:nRun) = struct('voxels',zeros(nVoxels,nTR),'ROI',zeros(1,nTR));
for irun=1:nRun

    for iVoxel=1:nVoxels
        
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
        
        Runs(irun).voxels(iVoxel,:) = squeeze(Condition(irun).run(idxx,idxy,idxz,:));
        
    end
    
    Runs(irun).ROI = mean(Runs(irun).voxels,1);
        
end

end