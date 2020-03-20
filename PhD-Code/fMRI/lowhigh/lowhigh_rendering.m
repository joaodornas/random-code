function lowhigh_rendering

%settings_jan_0502;

%settings_jan_0805;

settings_elena_2905;

doTheRendering(settings);

return

function doTheRendering(settings)

rendfile = getRenderFile(settings);

brt = 1;

outputfolder = getOutputFolder(settings);

zscorefolder = strcat(outputfolder,'\','zscore-whole-brain');
relativefolder = strcat(outputfolder,'\','relative-percentage');
fiveminfolder = strcat(outputfolder,'\','5-min');
stdfolder = strcat(outputfolder,'\','std');

%%% NAMES OF IMAGES

% disp('TMean Zscore Thresholded 0.05');
% disp('MOT4');
% img1_name = strcat('Low-High-',settings.folders.subject,'-MOT4-AllRuns-zscore-whole-brain-TMean-Thresholded-0.05.nii');
% disp('MOT2');
% img2_name = strcat('Low-High-',settings.folders.subject,'-MOT2-AllRuns-zscore-whole-brain-TMean-Thresholded-0.05.nii');
% disp('RestingState');
% img3_name = strcat('Low-High-',settings.folders.subject,'-RestingState-AllRuns-zscore-whole-brain-TMean-Thresholded-0.05.nii');
% 
% plot3imagesRendering(img1_name,img2_name,img3_name,zscorefolder,brt,rendfile);

disp('TMean Zscore Thresholded Positive 0.05');
disp('MOT4');
img1_name = strcat('Low-High-',settings.folders.subject,'-MOT4-AllRuns-zscore-whole-brain-TMean-Thresholded-Positive-0.05.nii');
disp('MOT2');
img2_name = strcat('Low-High-',settings.folders.subject,'-MOT2-AllRuns-zscore-whole-brain-TMean-Thresholded-Positive-0.05.nii');
disp('RestingState');
img3_name = strcat('Low-High-',settings.folders.subject,'-RestingState-AllRuns-zscore-whole-brain-TMean-Thresholded-Positive-0.05.nii');

plot3imagesRendering(img1_name,img2_name,img3_name,zscorefolder,brt,rendfile);

% disp('Diff TMean Zscore Thresholded 0.05');
% disp('MOT4 - RestingState');
% img1_name = strcat('Low-High-',settings.folders.subject,'-MOT4-RestingState-Diff','-','zscore-whole-brain-TMean-0.05','.nii');
% disp('MOT2 - RestingState');
% img2_name = strcat('Low-High-',settings.folders.subject,'-MOT2-RestingState-Diff','-','zscore-whole-brain-TMean-0.05','.nii');
% disp('MOT4 - MOT2');
% img3_name = strcat('Low-High-',settings.folders.subject,'-MOT4-MOT2-Diff','-','zscore-whole-brain-TMean-0.05','.nii');
% 
% plot3imagesRendering(img1_name,img2_name,img3_name,zscorefolder,brt,rendfile);
% 
% disp('Std');
% disp('MOT4');
% img1_name = strcat('Low-High-',settings.folders.subject,'-MOT4-AllRuns','-','std-time','.nii');
% disp('MOT2');
% img2_name = strcat('Low-High-',settings.folders.subject,'-MOT2-AllRuns','-','std-time','.nii');
% disp('RestingState');
% img3_name = strcat('Low-High-',settings.folders.subject,'-RestingState-AllRuns','-','std-time','.nii');
% 
% plot3imagesRendering(img1_name,img2_name,img3_name,stdfolder,brt,rendfile);
% 
% disp('Relative Percentage');
% disp('MOT4 - RestingState');
% img1_name = strcat('Low-High-',settings.folders.subject,'-MOT4-RestingState-Relative-Percentage','.nii');
% disp('MOT2 - RestingState');
% img2_name = strcat('Low-High-',settings.folders.subject,'-MOT2-RestingState-Relative-Percentage','.nii');
% disp('MOT4 - MOT2');
% img3_name = strcat('Low-High-',settings.folders.subject,'-MOT4-MOT2-Relative-Percentage','.nii');
% 
% plot3imagesRendering(img1_name,img2_name,img3_name,relativefolder,brt,rendfile);

% disp('Overall Mean');
% disp('MOT4');
% img1_name = strcat('Low-High-',settings.folders.subject,'-MOT4-Overall-Mean','.nii');
% disp('MOT2');
% img2_name = strcat('Low-High-',settings.folders.subject,'-MOT2-Overall-Mean','.nii');
% disp('RestingState');
% img3_name = strcat('Low-High-',settings.folders.subject,'-RestingState-Overall-Mean','.nii');
% 
% plot3imagesRendering(img1_name,img2_name,img3_name,relativefolder,brt,rendfile);
% 
% disp('5 minutes - Run1');
% disp('MOT4');
% img1_name = strcat('warfdn20_0544-0003-00150-000150-01','.img');
% disp('MOT2');
% img2_name = strcat('warfdn20_0544-0004-00150-000150-01','.img');
% disp('RestingState');
% img3_name = strcat('warfdn20_0544-0005-00150-000150-01','.img');
% 
% plot3imagesRendering(img1_name,img2_name,img3_name,fiveminfolder,brt,rendfile);

return

function plot3imagesRendering(img1_name,img2_name,img3_name,folder,brt,rendfile)

%%% NAMES OF IMAGES CONCATENATED

img1_file = strcat(folder,'\',img1_name);
img2_file = strcat(folder,'\',img2_name);
img3_file = strcat(folder,'\',img3_name);

%%% SPM DATA FILES

img1_dat = fromImageToSPMdat(img1_file);
img2_dat = fromImageToSPMdat(img2_file);
img3_dat = fromImageToSPMdat(img3_file);

%%% MIN AND MAX VALUES

[img1_min_value, img1_max_value] = getMinMaxFromData(img1_dat);
[img2_min_value, img2_max_value] = getMinMaxFromData(img2_dat);
[img3_min_value, img3_max_value] = getMinMaxFromData(img3_dat);

%%% MIN AND MAX VALUES OVERALL

min_value = floor(min([img1_min_value, img2_min_value, img3_min_value]));
max_value = ceil(max([img1_max_value, img2_max_value, img3_max_value]));

min_value = 1.96; %% zscore p-value 0.05
max_value = 2.16;

img1_dat.t(find(img1_dat.t>max_value)) = 0;
img2_dat.t(find(img2_dat.t>max_value)) = 0;
img3_dat.t(find(img3_dat.t>max_value)) = 0;

img1_dat.t(find(img1_dat.t<min_value)) = 0;
img2_dat.t(find(img2_dat.t<min_value)) = 0;
img3_dat.t(find(img3_dat.t<min_value)) = 0;

%%% SET GLOBAL MIN AND MAX

setGlobalMinMax(min_value,max_value);

%%% PLOT 3D RENDERING

plotRender(img1_dat,img2_dat,img3_dat,img1_file,img2_file,img3_file,img1_name,img2_name,img3_name,brt,rendfile);


return

function plotRender(img1_dat,img2_dat,img3_dat,img1_file,img2_file,img3_file,img1_name,img2_name,img3_name,brt,rendfile)

disp('Rendering Image 1');

snapshot_label = img1_name;
setRenderingSnapshot(true,snapshot_label);

spm_render_my(img1_dat,brt,rendfile,img1_file);

%pause;

disp('Rendering Image 2');

snapshot_label = img2_name;
setRenderingSnapshot(true,snapshot_label);

spm_render_my(img2_dat,brt,rendfile,img2_file);

%pause;

disp('Rendering Image 3');

snapshot_label = img3_name;
setRenderingSnapshot(true,snapshot_label);

spm_render_my(img3_dat,brt,rendfile,img3_file);

return

function outputfolder = getOutputFolder(settings)
        
outputfolder = strcat(settings.folders.main,'\',...
    settings.folders.experiment,'\',...
    settings.folders.subject,'\',...
    settings.folders.output);
        
return

function [min_value, max_value] = getMinMaxFromData(dat)

min_value = min(dat.t);
max_value = max(dat.t);

return


