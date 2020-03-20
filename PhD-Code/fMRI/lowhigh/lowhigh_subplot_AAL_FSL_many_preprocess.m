
function lowhigh_subplot_AAL_FSL_many_preprocess

%%% SETTINGS

%settings_jan_0805;
settings_elena_2905;

doTheMath(settings);


return

function doTheMath(settings)

%%% LOAD DATA

%file = settings.FSL.files.functional.warped;
%file = settings.FSL.files.functional.despike;
file{1} = settings.FSL.files.functional.custom.filtered;
%file{2} = settings.FSL.files.functional.custom.wavelet;
%file{3} = settings.FSL.files.functional.custom.residual;
%file{2} = settings.FSL.files.functional.custom.wavelet_voxel;
file{2} = settings.FSL.files.functional.custom.residual_voxel;
%file{2} = settings.FSL.files.functional.custom.noise_voxel;

name{1} = 'filtered';
%name{2} = 'wavelet';
%name{3} = 'residual';
%name{2} = 'wavelet-voxel';
name{2} = 'residual-voxel';
%name{2} = 'noise-voxel';

datalabel = 'many';

get_at_this_preprocessed_step = settings.FSL.folders.custom;
mask = settings.FSL.files.mask.custom;

%lowhigh_load_all_data_FSL;

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

%nNodes = 90;

%% AREAS

which_mean = 1; %% AAL mean time series
%which_mean = 2; %% AAL highest mean voxel

arealabel = 'Frontal';
disp(arealabel);
l = 6;
c = 3;
idx_nodes = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];
run = 1;
doThePlot('MOT4',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('MOT2',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('RestingState',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
run = 2;
doThePlot('MOT4',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('MOT2',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('RestingState',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);

arealabel = 'Occipital';
disp(arealabel);
l = 3;
c = 2;
idx_nodes = [49 50 51 52 53 54];
run = 1;
doThePlot('MOT4',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('MOT2',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('RestingState',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
run = 2;
doThePlot('MOT4',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('MOT2',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('RestingState',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);

arealabel = 'Parietal';
disp(arealabel);
l = 2;
c = 2;
idx_nodes = [59 60 61 62];
run = 1;
doThePlot('MOT4',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('MOT2',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('RestingState',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
run = 2;
doThePlot('MOT4',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('MOT2',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('RestingState',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);

arealabel = 'Temporal';
disp(arealabel);
l = 4;
c = 3;
idx_nodes = [81 82 83 84 85 86 87 88 89 90];
run = 1;
doThePlot('MOT4',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('MOT2',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('RestingState',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
run = 2;
doThePlot('MOT4',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('MOT2',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);
doThePlot('RestingState',run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask);

return

function doThePlot(kind,run,file,name,AAL_img,AAL_ROI,settings,idx_nodes,arealabel,datalabel,l,c,which_mean,get_at_this_preprocessed_step,mask)

nNodes = length(idx_nodes);

for ifile=1:length(file)
    
    [Run(ifile).ts, maskdata(ifile).mask, settings] = lowhigh_get_data_FSL(settings,kind,run,file{ifile},mask,get_at_this_preprocessed_step);

end

colors{1} = 'b';
colors{2} = 'r';
colors{3} = 'y';

% c = 4;
% l = ceil(nNodes / c);

f = figure;

fs = 6;

ivalue_mean = 0;
ivalue_hm_mean = 0;

for iNode=1:nNodes
    
   idx_region = AAL_ROI(idx_nodes(iNode)).ID;
   label = AAL_ROI(idx_nodes(iNode)).Nom_L;
   
   label = strrep(label,'_','-');
   
   idx_voxels_structures = find(AAL_img == idx_region);
   nVoxels = length(idx_voxels_structures);
   
   for iVoxel=1:nVoxels
        
        [idxx(iVoxel), idxy(iVoxel), idxz(iVoxel)] = ind2sub(size(AAL_img),idx_voxels_structures(iVoxel));
    
   end  
   
   nTR = size(Run(1).ts,4);
    
   for ifile=1:length(file)
       
    Run_voxel_time_series(ifile).ts = zeros(nVoxels,nTR);
    
   end
    
   for iVoxel=1:nVoxels
      
       for ifile=1:length(file)
           
            Run_voxel_time_series(ifile).ts(iVoxel,1:nTR) = Run(ifile).ts(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
       end
       
   end
    
   for ifile=1:length(file)
       
        ROI(iNode).Run_AAL(ifile).mean = mean(Run_voxel_time_series(ifile).ts,1);
   
        ivalue_mean = ivalue_mean + 1;
        all_min_value_mean(ivalue_mean) = min(ROI(iNode).Run_AAL(ifile).mean);
        all_max_value_mean(ivalue_mean) = max(ROI(iNode).Run_AAL(ifile).mean);
        
   end
   
   for ifile=1:length(file)
       
        Run_AAL(ifile).mean_time = mean(Run_voxel_time_series(ifile).ts,2);
   
        [x,I] = sort(Run_AAL(ifile).mean_time);
        idx_hm = I(end);
        ROI(iNode).Run_AAL_voxel(ifile).ts = Run_voxel_time_series(ifile).ts(idx_hm,:);
   
        ivalue_hm_mean = ivalue_hm_mean + 1;
        all_min_value_hm_mean(ivalue_hm_mean) = min(ROI(iNode).Run_AAL_voxel(ifile).ts);
        all_max_value_hm_mean(ivalue_hm_mean) = max(ROI(iNode).Run_AAL_voxel(ifile).ts);
        
   end
   
end

for iNode=1:nNodes
   
   subplot(l,c,iNode);
   
   if which_mean == 1
       
       for ifile=1:length(file)
           
            my_plot(ifile) = plot(ROI(iNode).Run_AAL(ifile).mean,colors{ifile});
            hold on;
        
       end
       
   elseif which_mean == 2
       
      for ifile=1:length(file)
           
            my_plot(ifile) = plot(ROI(iNode).Run_AAL_voxel(ifile).ts,colors{ifile});
            hold on;
        
       end
       
   end
   
   title(strcat(kind,'-Run-',int2str(run),'-',label),'FontSize',fs);
   
   xlabel('TR','FontSize',fs);
   ylabel('BOLD Activity','FontSize',fs);
   hL = legend(my_plot,name);
   
   newPosition = [0.1 0.1 0.1 0.1];
   newUnits = 'normalized';
   set(hL,'Position',newPosition,'Units',newUnits);
   
   set(gca,'FontSize',fs);
   
   xlim([0 nTR]);
   
%    if which_mean == 1
%        
%         ylim([min(all_min_value_mean) max(all_max_value_mean)]);
%    
%    elseif which_mean == 2
%        
%         ylim([min(all_min_value_hm_mean) max(all_max_value_hm_mean)]);
%        
%    end
   
end

if which_mean == 1
    
    which_mean_label = 'Mean';
    
elseif which_mean == 2
    
    which_mean_label = 'HM-Voxel';
    
end

% print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',kind,'-Run-',int2str(run),'-AAL-time-series-',which_mean_label,'-',arealabel,'-',datalabel,'.jpg'));
% print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',kind,'-Run-',int2str(run),'-AAL-time-series-',which_mean_label,'-',arealabel,'-',datalabel,'.eps'));
% print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',kind,'-Run-',int2str(run),'-AAL-time-series-',which_mean_label,'-',arealabel,'-',datalabel,'.pdf'));

print(f,'-djpeg',strcat('LH','-',settings.folders.subject,'-',kind,'-Run-',int2str(run),'-AAL-',which_mean_label,'-',arealabel,'-',datalabel,'.jpg'));
print(f,'-depsc',strcat('LH','-',settings.folders.subject,'-',kind,'-Run-',int2str(run),'-AAL-',which_mean_label,'-',arealabel,'-',datalabel,'.eps'));
print(f,'-dpdf',strcat('LH','-',settings.folders.subject,'-',kind,'-Run-',int2str(run),'-AAL-',which_mean_label,'-',arealabel,'-',datalabel,'.pdf'));

close all;

return