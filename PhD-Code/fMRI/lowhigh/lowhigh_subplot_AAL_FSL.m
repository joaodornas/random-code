
function lowhigh_subplot_AAL_FSL

%%% SETTINGS

settings_jan_0805;
%settings_elena_2905;

doTheMath(settings);


return

function doTheMath(settings)

%%% LOAD DATA

%file = settings.FSL.files.functional.warped;
%file = settings.FSL.files.functional.despike;
%file = settings.FSL.files.functional.custom.filtered;
%file = settings.FSL.files.functional.custom.wavelet;
%file = settings.FSL.files.functional.custom.residual;
%file = settings.FSL.files.functional.custom.wavelet_voxel;
%file = settings.FSL.files.functional.custom.noise_voxel;
file = settings.FSL.files.functional.custom.residual_voxel;
%datalabel = 'warped';
%datalabel = 'despike';
%datalabel = 'filtered';
%datalabel = 'wavelet';
%datalabel = 'residual';
datalabel = 'residual-voxel';
%datalabel = 'noise-voxel';

get_at_this_preprocessed_step = settings.FSL.folders.custom;
mask = settings.FSL.files.mask.custom;

lowhigh_load_all_data_FSL;

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

%nNodes = 90;

%% AREAS

which_mean = 1; %% AAL mean time series
%which_mean = 2; %% AAL highest mean voxel

MOT4Run = ( MOT4Run1 + MOT4Run2 ) ./ 2;
MOT2Run = ( MOT2Run1 + MOT2Run2 ) ./ 2;
RestingStateRun = ( RestingStateRun1 + RestingStateRun2 ) ./ 2;

arealabel = 'Frontal';
disp(arealabel);
l = 6;
c = 3;
idx_nodes = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];
%doThePlot(MOT4Run1,MOT2Run1,RestingStateRun1,AAL_img,AAL_ROI,settings,idx_nodes,1,arealabel,datalabel,l,c,which_mean);
%doThePlot(MOT4Run2,MOT2Run2,RestingStateRun2,AAL_img,AAL_ROI,settings,idx_nodes,2,arealabel,datalabel,l,c,which_mean);

doThePlot(MOT4Run,MOT2Run,RestingStateRun,AAL_img,AAL_ROI,settings,idx_nodes,12,arealabel,datalabel,l,c,which_mean);

arealabel = 'Occipital';
disp(arealabel);
l = 3;
c = 2;
idx_nodes = [49 50 51 52 53 54];
%doThePlot(MOT4Run1,MOT2Run1,RestingStateRun1,AAL_img,AAL_ROI,settings,idx_nodes,1,arealabel,datalabel,l,c,which_mean);
%doThePlot(MOT4Run2,MOT2Run2,RestingStateRun2,AAL_img,AAL_ROI,settings,idx_nodes,2,arealabel,datalabel,l,c,which_mean);

doThePlot(MOT4Run,MOT2Run,RestingStateRun,AAL_img,AAL_ROI,settings,idx_nodes,12,arealabel,datalabel,l,c,which_mean);

arealabel = 'Parietal';
disp(arealabel);
l = 2;
c = 2;
idx_nodes = [59 60 61 62];
%doThePlot(MOT4Run1,MOT2Run1,RestingStateRun1,AAL_img,AAL_ROI,settings,idx_nodes,1,arealabel,datalabel,l,c,which_mean);
%doThePlot(MOT4Run2,MOT2Run2,RestingStateRun2,AAL_img,AAL_ROI,settings,idx_nodes,2,arealabel,datalabel,l,c,which_mean);

doThePlot(MOT4Run,MOT2Run,RestingStateRun,AAL_img,AAL_ROI,settings,idx_nodes,12,arealabel,datalabel,l,c,which_mean);

arealabel = 'Temporal';
disp(arealabel);
l = 4;
c = 3;
idx_nodes = [81 82 83 84 85 86 87 88 89 90];
%doThePlot(MOT4Run1,MOT2Run1,RestingStateRun1,AAL_img,AAL_ROI,settings,idx_nodes,1,arealabel,datalabel,l,c,which_mean);
%doThePlot(MOT4Run2,MOT2Run2,RestingStateRun2,AAL_img,AAL_ROI,settings,idx_nodes,2,arealabel,datalabel,l,c,which_mean);

doThePlot(MOT4Run,MOT2Run,RestingStateRun,AAL_img,AAL_ROI,settings,idx_nodes,12,arealabel,datalabel,l,c,which_mean);

return

function doThePlot(MOT4Run,MOT2Run,RestingStateRun,AAL_img,AAL_ROI,settings,idx_nodes,run,arealabel,datalabel,l,c,which_mean)

nNodes = length(idx_nodes);

colors{1} = 'b';
colors{2} = 'r';
colors{3} = 'y';

% % c = 4;
% l = ceil(nNodes / c);

f = figure;

fs = 6;

for iNode=1:nNodes
    
   idx_region = AAL_ROI(idx_nodes(iNode)).ID;
   label = AAL_ROI(idx_nodes(iNode)).Nom_L;
   
   label = strrep(label,'_','-');
   
   idx_voxels_structures = find(AAL_img == idx_region);
   nVoxels = length(idx_voxels_structures);
   
   for iVoxel=1:nVoxels
        
        [idxx(iVoxel), idxy(iVoxel), idxz(iVoxel)] = ind2sub(size(AAL_img),idx_voxels_structures(iVoxel));
    
   end  
   
   nTR = size(MOT4Run,4);
    
   MOT4Run_voxel_time_series = zeros(nVoxels,nTR);
   MOT2Run_voxel_time_series = zeros(nVoxels,nTR);
   RestingStateRun_voxel_time_series = zeros(nVoxels,nTR);
    
   for iVoxel=1:nVoxels
        
        MOT4Run_voxel_time_series(iVoxel,1:nTR) = MOT4Run(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        MOT2Run_voxel_time_series(iVoxel,1:nTR) = MOT2Run(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        RestingStateRun_voxel_time_series(iVoxel,1:nTR) = RestingStateRun(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
    
   end
    
   MOT4Run_AAL = mean(MOT4Run_voxel_time_series,1);
   MOT2Run_AAL = mean(MOT2Run_voxel_time_series,1);
   RestingStateRun_AAL = mean(RestingStateRun_voxel_time_series,1);
   
   MOT4Run_AAL_mean_time = mean(MOT4Run_voxel_time_series,2);
   MOT2Run_AAL_mean_time = mean(MOT2Run_voxel_time_series,2);
   RestingStateRun_AAL_mean_time = mean(RestingStateRun_voxel_time_series,2);
   
   [x,I] = sort(MOT4Run_AAL_mean_time);
   idx_hm_mot4 = I(end);
   MOT4Run_AAL_voxel = MOT4Run_voxel_time_series(idx_hm_mot4,:);
   
   [x,I] = sort(MOT2Run_AAL_mean_time);
   idx_hm_mot2 = I(end);
   MOT2Run_AAL_voxel = MOT2Run_voxel_time_series(idx_hm_mot2,:);
   
   [x,I] = sort(RestingStateRun_AAL_mean_time);
   idx_hm_rs = I(end);
   RestingStateRun_AAL_voxel = RestingStateRun_voxel_time_series(idx_hm_rs,:);
   
   subplot(l,c,iNode);
   
   if which_mean == 1
       
        plot(MOT4Run_AAL,colors{1});
        hold on;
        plot(MOT2Run_AAL,colors{2});
        hold on;
        plot(RestingStateRun_AAL,colors{3});
        hold on;
   
   elseif which_mean == 2
       
        plot(MOT4Run_AAL_voxel,colors{1});
        hold on;
        plot(MOT2Run_AAL_voxel,colors{2});
        hold on;
        plot(RestingStateRun_AAL_voxel,colors{3});
        hold on;
       
   end
   
   title(strcat(label,'-','Run-',int2str(run)),'FontSize',fs);
   
   xlabel('TR','FontSize',fs);
   ylabel('BOLD Activity','FontSize',fs);
   %legend('High Attention','Low Attention','Resting State');
   
   set(gca,'FontSize',fs);
   
   xlim([0 nTR]);
   
end

if which_mean == 1
    
    which_mean_label = 'Mean';
    
elseif which_mean == 2
    
    which_mean_label = 'Highest-Mean-Voxel';
    
end

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','AAL-time-series-',which_mean_label,'-Run-',int2str(run),'-',arealabel,'-',datalabel,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','AAL-time-series-',which_mean_label,'-Run-',int2str(run),'-',arealabel,'-',datalabel,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','AAL-time-series-',which_mean_label,'-Run-',int2str(run),'-',arealabel,'-',datalabel,'.pdf'));

close all;

return