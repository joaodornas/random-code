
function lowhigh_zscore_whole_brain_diff_FSL

%settings_jan_0502;

%settings_jan_0805;

settings_elena_2905;

doTheMath(settings);

return

function doTheMath(settings)

lowhigh_load_all_data_FSL;

scale = 9;

pvalue(1) = 0.05;
pvalue(2) = 0.001;

disp('get Common Mask');

masks(1).mask = mask_MOT4Run1;
masks(2).mask = mask_MOT4Run2;
masks(3).mask = mask_MOT2Run1;
masks(4).mask = mask_MOT2Run2;
masks(5).mask = mask_RestingStateRun1;
masks(6).mask =  mask_RestingStateRun2;

common_mask = get_common_mask(masks);

saveCommonMask(common_mask,settings);

disp('get ZScore');

[MOT4_Tmean_WB_ZScore, MOT4_Tstd_WB_ZScore, MOT4_Tstd_WB_ZScore_twice] = getZScore(MOT4Run1,MOT4Run2,common_mask,common_mask,scale);

[MOT2_Tmean_WB_ZScore, MOT2_Tstd_WB_ZScore, MOT2_Tstd_WB_ZScore_twice] = getZScore(MOT2Run1,MOT2Run2,common_mask,common_mask,scale);

[RestingState_Tmean_WB_ZScore, RestingState_Tstd_WB_ZScore, RestingState_Tstd_WB_ZScore_twice] = getZScore(RestingStateRun1,RestingStateRun2,common_mask,common_mask,scale);

disp('save ZScore');

kind = 'TMean';
saveZScore(MOT4_Tmean_WB_ZScore,MOT2_Tmean_WB_ZScore,RestingState_Tmean_WB_ZScore,settings,scale,kind);
kind = 'Tstd';
saveZScore(MOT4_Tstd_WB_ZScore,MOT2_Tstd_WB_ZScore,RestingState_Tstd_WB_ZScore,settings,scale,kind);
kind = 'Tstd-twice';
saveZScore(MOT4_Tstd_WB_ZScore_twice,MOT2_Tstd_WB_ZScore_twice,RestingState_Tstd_WB_ZScore_twice,settings,scale,kind);

disp('get Differences');

common_idx = find(common_mask < 2);

[MOT4RestingStateDiff_positive, MOT4RestingStateDiff_negative, MOT4RestingStateDiff] = getDiff(MOT4_Tmean_WB_ZScore,RestingState_Tmean_WB_ZScore,common_mask,common_idx);

[MOT2RestingStateDiff_positive, MOT2RestingStateDiff_negative, MOT2RestingStateDiff] = getDiff(MOT2_Tmean_WB_ZScore,RestingState_Tmean_WB_ZScore,common_mask,common_idx);

[MOT4MOT2Diff_positive, MOT4MOT2Diff_negative, MOT4MOT2Diff] = getDiff(MOT4_Tmean_WB_ZScore,MOT2_Tmean_WB_ZScore,common_mask,common_idx);

kind = 'no threshold';
common_idx = 0;
plotHistogram(MOT4_Tmean_WB_ZScore,MOT2_Tmean_WB_ZScore,RestingState_Tmean_WB_ZScore,MOT4RestingStateDiff,MOT2RestingStateDiff,MOT4MOT2Diff,0,kind,common_mask,common_idx,settings,scale);

disp('Threshold Images');

for ip=1:1
    
    disp(strcat('p-value:',num2str(pvalue(ip))));
    
    disp('get Std Threshold');

    [MOT4_Std_Thresholded, MOT4_Std_Thresholded_positive, MOT4_Std_Thresholded_negative, idx_threshold_MOT4_Std_positive, idx_threshold_MOT4_Std_negative, idx_threshold_MOT4_Std] = getThresholdedImages(MOT4_Tstd_WB_ZScore_twice,pvalue(ip),0);
    [MOT2_Std_Thresholded, MOT2_Std_Thresholded_positive, MOT2_Std_Thresholded_negative, idx_threshold_MOT2_Std_positive, idx_threshold_MOT2_Std_negative, idx_threshold_MOT2_Std] = getThresholdedImages(MOT2_Tstd_WB_ZScore_twice,pvalue(ip),0);
    [RestingState_Std_Thresholded, RestingState_Std_Thresholded_positive, RestingState_Std_Thresholded_negative, idx_threshold_RestingState_Std_positive, idx_threshold_RestingState_Std_negative, idx_threshold_RestingState_Std] = getThresholdedImages(RestingState_Tstd_WB_ZScore_twice,pvalue(ip),0);
    
    disp('save Std Threshold');

    kind = 'Thresholded-';
    saveStdThreshold(MOT4_Std_Thresholded,MOT2_Std_Thresholded,RestingState_Std_Thresholded,pvalue(1),settings,kind);
    kind = 'Thresholded-Positive-';
    saveStdThreshold(MOT4_Std_Thresholded_positive,MOT2_Std_Thresholded_positive,RestingState_Std_Thresholded_positive,pvalue(1),settings,kind);
    kind = 'Thresholded-Negative-';
    saveStdThreshold(MOT4_Std_Thresholded_negative,MOT2_Std_Thresholded_negative,RestingState_Std_Thresholded_negative,pvalue(1),settings,kind);

    MOT4_Std_Thresholded_Negative_Inverted = MOT4_Std_Thresholded_negative .* (-1);
    MOT2_Std_Thresholded_Negative_Inverted = MOT2_Std_Thresholded_negative .* (-1);
    RestingState_Std_Thresholded_Negative_Inverted = RestingState_Std_Thresholded_negative .* (-1);
    
    kind = 'Thresholded-Negative-Inverted-';
    saveStdThreshold(MOT4_Std_Thresholded_Negative_Inverted,MOT2_Std_Thresholded_Negative_Inverted,RestingState_Std_Thresholded_Negative_Inverted,pvalue(1),settings,kind);
    
    disp('get TMean Threshold');
    
    [MOT4_Threshold, MOT4_Threshold_Positive, MOT4_Threshold_Negative, idx_threshold_MOT4_Positive, idx_threshold_MOT4_Negative, idx_threshold_MOT4] = getThresholdedImages(MOT4_Tmean_WB_ZScore,pvalue(ip),scale);

    MOT4_Threshold_Negative_Inverted = MOT4_Threshold_Negative .* (-1);
    
    [MOT2_Threshold, MOT2_Threshold_Positive, MOT2_Threshold_Negative, idx_threshold_MOT2_Positive, idx_threshold_MOT2_Negative, idx_threshold_MOT2] = getThresholdedImages(MOT2_Tmean_WB_ZScore,pvalue(ip),scale);

    MOT2_Threshold_Negative_Inverted = MOT2_Threshold_Negative .* (-1);
    
    [RestingState_Threshold, RestingState_Threshold_Positive, RestingState_Threshold_Negative, idx_threshold_RestingState_Positive, idx_threshold_RestingState_Negative, idx_threshold_RestingState] = getThresholdedImages(RestingState_Tmean_WB_ZScore,pvalue(ip),scale);

    RestingState_Threshold_Negative_Inverted = RestingState_Threshold_Negative .* (-1);
    
    common_idx = getCommonIdx(idx_threshold_MOT4, idx_threshold_MOT2, idx_threshold_RestingState);
    
    disp('save TMean Threshold');
    
    kind = 'MOT4';
    saveThresholdImages(MOT4_Threshold_Positive,MOT4_Threshold_Negative,MOT4_Threshold_Negative_Inverted,MOT4_Threshold,kind,pvalue(ip),settings,scale);
    
    kind = 'MOT2';
    saveThresholdImages(MOT2_Threshold_Positive,MOT2_Threshold_Negative,MOT2_Threshold_Negative_Inverted,MOT2_Threshold,kind,pvalue(ip),settings,scale);
    
    kind = 'RestingState';
    saveThresholdImages(RestingState_Threshold_Positive,RestingState_Threshold_Negative,RestingState_Threshold_Negative_Inverted,RestingState_Threshold,kind,pvalue(ip),settings,scale);
    
    disp('get Difference');
    
    [MOT4RestingStateDiff_positive, MOT4RestingStateDiff_negative, MOT4RestingStateDiff] = getDiff(MOT4_Tmean_WB_ZScore,RestingState_Tmean_WB_ZScore,common_mask,common_idx);

    [MOT2RestingStateDiff_positive, MOT2RestingStateDiff_negative, MOT2RestingStateDiff] = getDiff(MOT2_Tmean_WB_ZScore,RestingState_Tmean_WB_ZScore,common_mask,common_idx);

    [MOT4MOT2Diff_positive, MOT4MOT2Diff_negative, MOT4MOT2Diff] = getDiff(MOT4_Tmean_WB_ZScore,MOT2_Tmean_WB_ZScore,common_mask,common_idx);
    
    disp('save Difference');
    
    MOT4RestingStateDiff_negative_inverted = MOT4RestingStateDiff_negative .* (-1);
    MOT2RestingStateDiff_negative_inverted = MOT2RestingStateDiff_negative .* (-1);
    MOT4MOT2Diff_negative_inverted = MOT4MOT2Diff_negative .* (-1);
    
    saveDiff(MOT4RestingStateDiff_positive, MOT4RestingStateDiff_negative, MOT4RestingStateDiff_negative_inverted, MOT4RestingStateDiff,MOT2RestingStateDiff_positive, MOT2RestingStateDiff_negative, MOT2RestingStateDiff_negative_inverted, MOT2RestingStateDiff, MOT4MOT2Diff_positive, MOT4MOT2Diff_negative, MOT4MOT2Diff_negative_inverted, MOT4MOT2Diff,pvalue(ip),settings,scale);
    
    disp('Plot Histogram');
    
    kind = 'threshold';
    plotHistogram(MOT4_Threshold,MOT2_Threshold,RestingState_Threshold,MOT4RestingStateDiff,MOT2RestingStateDiff,MOT4MOT2Diff,pvalue(ip),kind,common_mask,common_idx,settings,scale);
    
end

return

function [Tmean_WB_ZScore, Tstd_WB_ZScore, Tstd_WB_ZScore_twice] = getZScore(Run1,Run2,mask_Run1,mask_Run2,scale)

[Run1ZStint, Run1ZStdouble] = zscore_whole_brain(Run1,mask_Run1,scale);
[Run2ZStint, Run2ZStdouble] = zscore_whole_brain(Run2,mask_Run2,scale);

Run1ZStmean = mean(Run1ZStint,4,'native');
Run2ZStmean = mean(Run2ZStint,4,'native');

Run1ZStstd = std(Run1,0,4);
Run2ZStstd = std(Run2,0,4);

masks(1).mask = mask_Run1;
masks(2).mask = mask_Run2;

common_mask = get_common_mask(masks);
idx_all = find(common_mask < 2);
idx_all(find(common_mask == 1)) = [];

Run1ZStmean(idx_all) = 0;
Run2ZStmean(idx_all) = 0;

Run1ZStstd(idx_all) = 0;
Run2ZStstd(idx_all) = 0;

ZStAll = int32(zeros(size(Run1,1),size(Run1,2),size(Run1,3),2));
ZStstdAll = double(zeros(size(Run1,1),size(Run1,2),size(Run1,3),2));

ZStAll(:,:,:,1) = Run1ZStmean(:,:,:);
ZStAll(:,:,:,2) = Run2ZStmean(:,:,:);

ZStstdAll(:,:,:,1) = Run1ZStstd(:,:,:);
ZStstdAll(:,:,:,2) = Run2ZStstd(:,:,:);

Tmean_WB_ZScore = mean(ZStAll,4,'native');

Tstd_WB_ZScore = mean(ZStstdAll,4);

ZStstdAllmean = mean(ZStstdAll,4);

get_std_data_vec = reshape(ZStstdAllmean,[1,size(ZStstdAllmean,1)*size(ZStstdAllmean,2)*size(ZStstdAllmean,3)]);
mask_vec = reshape(common_mask,[1,size(common_mask,1)*size(common_mask,2)*size(common_mask,3)]);
idx_voxels_brain = find(mask_vec == 1);
idx_voxels_space = find(mask_vec == 0);

voxels_data = get_std_data_vec(idx_voxels_brain);
idx_non_zero = find(voxels_data ~= 0);
voxels_data_zscored = zeros(1,length(voxels_data));
voxels_data_zscored(idx_non_zero) = zscore(voxels_data(idx_non_zero));
new_get_std_data_vec = zeros(1,size(common_mask,1)*size(common_mask,2)*size(common_mask,3));
new_get_std_data_vec(idx_voxels_brain) = voxels_data_zscored;

Tstd_WB_ZScore_twice = reshape(new_get_std_data_vec,[size(common_mask,1),size(common_mask,2),size(common_mask,3)]);
    
return

function [img_thresholded, img_thresholded_positive, img_thresholded_negative, idx_positive, idx_negative, idx_all] = getThresholdedImages(zscored_image,pvalue,scale)

%     scale = 4;

%     if pvalue == 0.05
%         
%         threshold_positive = 1.96;
%         threshold_negative = -1.96;
%         
%     elseif pvalue == 0.001
%         
%         threshold_positive = 3.2905;
%         threshold_negative = -3.2905;
%         
%     end

    threshold_positive = abs(norminv(pvalue/2,0,1));
    threshold_negative = norminv(pvalue/2,0,1);
    
    threshold_positive = threshold_positive*10^(scale);
    threshold_negative = threshold_negative*10^(scale);
    
    img_thresholded_positive = zscored_image;
    idx_pos = find(img_thresholded_positive<threshold_positive);
    img_thresholded_positive(idx_pos) = 0;

    img_thresholded_negative = zscored_image;
    idx_neg = find(img_thresholded_negative>threshold_negative);
    img_thresholded_negative(idx_neg) = 0;
    
    img_thresholded = zeros(size(zscored_image,1),size(zscored_image,2),size(zscored_image,3));
    img_thresholded = img_thresholded_positive + img_thresholded_negative;

    idx_positive = find(img_thresholded_positive);
    idx_negative = find(img_thresholded_negative);
    if size(idx_positive,1) ~= 1, idx_positive = idx_positive'; end
    if size(idx_negative,1) ~= 1, idx_negative = idx_negative'; end
    idx_all = [idx_positive, idx_negative];
    
return

function saveZScore(MOT4,MOT2,RestingState,settings,scale,kind)

run = 1;

dim = [settings.mot4.FSL.run(run).dim(1), settings.mot4.FSL.run(run).dim(2),settings.mot4.FSL.run(run).dim(3)];

nifti_file = settings.mot4.FSL.run(run).nifti_file;
dtype = settings.mot4.FSL.run(run).dtype;
offset = settings.mot4.FSL.run(run).offset;
scl_slope = settings.mot4.FSL.run(run).scl_slope;
scl_inter = settings.mot4.FSL.run(run).scl_inter;
    
dtype = 'FLOAT32';
offset = 0;

descrip = 'mean Z-Score';
 
fname = strcat('Low-High-',settings.folders.subject,'-MOT4-AllRuns-','zscore-whole-brain-',kind,'.nii');

input_data = MOT4; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;     

fname = strcat('Low-High-',settings.folders.subject,'-MOT2-AllRuns-','zscore-whole-brain-',kind,'.nii');

input_data = MOT2; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;   

fname = strcat('Low-High-',settings.folders.subject,'-RestingState-AllRuns-','zscore-whole-brain-',kind,'.nii');

input_data = RestingState; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;   

return

function saveThresholdImages(img_thresholded_positive,img_thresholded_negative,img_thresholded_negative_inverted,img_thresholded,kind,pvalue,settings,scale)

run = 1;

dim = [settings.mot4.FSL.run(run).dim(1), settings.mot4.FSL.run(run).dim(2),settings.mot4.FSL.run(run).dim(3)];

nifti_file = settings.mot4.FSL.run(run).nifti_file;
dtype = settings.mot4.FSL.run(run).dtype;
offset = settings.mot4.FSL.run(run).offset;
scl_slope = settings.mot4.FSL.run(run).scl_slope;
scl_inter = settings.mot4.FSL.run(run).scl_inter;
    
dtype = 'FLOAT32';
offset = 0;

descrip = 'Z-Score Thresholded';
 
fname = strcat('Low-High-',settings.folders.subject,'-',kind,'-AllRuns-','zscore-whole-brain-TMean-Thresholded-Positive-',num2str(pvalue),'.nii');

input_data = img_thresholded_positive; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;     

fname = strcat('Low-High-',settings.folders.subject,'-',kind,'-AllRuns-','zscore-whole-brain-TMean-Thresholded-Negative-',num2str(pvalue),'.nii');

input_data = img_thresholded_negative;
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;   

fname = strcat('Low-High-',settings.folders.subject,'-',kind,'-AllRuns-','zscore-whole-brain-TMean-Thresholded-Negative-Inverted-',num2str(pvalue),'.nii');

input_data = img_thresholded_negative_inverted; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image; 

fname = strcat('Low-High-',settings.folders.subject,'-',kind,'-AllRuns-','zscore-whole-brain-TMean-Thresholded-',num2str(pvalue),'.nii');

input_data = img_thresholded;  
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;   

return

function saveCommonMask(common_mask,settings)

run = 1;

dim = [settings.mot4.FSL.run(run).dim(1), settings.mot4.FSL.run(run).dim(2),settings.mot4.FSL.run(run).dim(3)];

nifti_file = settings.mot4.FSL.run(run).nifti_file;
dtype = settings.mot4.FSL.run(run).dtype;
offset = settings.mot4.FSL.run(run).offset;
scl_slope = settings.mot4.FSL.run(run).scl_slope;
scl_inter = settings.mot4.FSL.run(run).scl_inter;
    
dtype = 'FLOAT32';
offset = 0;

descrip = 'Common Mask';
 
fname = strcat('Low-High-',settings.folders.subject,'-','CommonMask-AllRuns','.nii');

input_data = common_mask; 

lowhigh_save_image;     

return

function common_idx = getCommonIdx(MOT4_idx, MOT2_idx, RestingState_idx)

MOT4MOT2_idx = MOT4_idx(ismember(MOT4_idx,MOT2_idx));

MOT4MOT2RestingState_idx = MOT4MOT2_idx(ismember(MOT4MOT2_idx,RestingState_idx));

common_idx = MOT4MOT2RestingState_idx;

return

function [difference_image_positive, difference_image_negative, difference_image] = getDiff(image_1,image_2,mask,idx_threshold)

    idx_mask = find(mask);
    
    difference_image =  int32(zeros(size(image_1,1),size(image_1,2),size(image_1,3)));
    
%     %%% remove voxels out of mask
%     
%     all_idx = find(mask<2);
%     
%     all_idx(idx_mask) = [];
%     
%     image_1(all_idx) = 0;
%     
%     image_2(all_idx) = 0;
%     
%     %%% remove non statistically significant voxels
%     
%     all_idx = find(mask<2);
%     
%     all_idx(idx_threshold) = [];
%     
%     image_1(all_idx) = 0;
%     
%     image_2(all_idx) = 0;
    
    %%% get only relevant indices
    
    idx_mask = idx_mask(ismember(idx_mask,idx_threshold));
    
    difference_image(idx_mask) = image_1(idx_mask) - image_2(idx_mask);
    
    difference_image_positive = difference_image;
    idx_pos = find(difference_image_positive<0);
    difference_image_positive(idx_pos) = 0;
    
    difference_image_negative = difference_image;
    idx_neg = find(difference_image_negative>0);
    difference_image_negative(idx_neg) = 0;
    
return

function saveDiff(MOT4RestingStateDiff_positive, MOT4RestingStateDiff_negative, MOT4RestingStateDiff_negative_inverted, MOT4RestingStateDiff,MOT2RestingStateDiff_positive, MOT2RestingStateDiff_negative, MOT2RestingStateDiff_negative_inverted, MOT2RestingStateDiff, MOT4MOT2Diff_positive, MOT4MOT2Diff_negative, MOT4MOT2Diff_negative_inverted, MOT4MOT2Diff,pvalue,settings,scale)

run = 1;
nifti_file = settings.restingstate.FSL.run(run).nifti_file;
dtype = settings.restingstate.FSL.run(run).dtype;
offset = settings.restingstate.FSL.run(run).offset;
scl_slope = settings.restingstate.FSL.run(run).scl_slope;
scl_inter = settings.restingstate.FSL.run(run).scl_inter;
    
dim = [settings.restingstate.FSL.run(run).dim(1), settings.restingstate.FSL.run(run).dim(2),settings.restingstate.FSL.run(run).dim(3)];

dtype = 'FLOAT32';
offset = 0;

descrip = 'mean ZScore Diff';

fname = strcat('Low-High-',settings.folders.subject,'-MOT4-RestingState-Diff','-','zscore-whole-brain-TMean-',num2str(pvalue),'.nii');

input_data = MOT4RestingStateDiff; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT2-RestingState-Diff','-','zscore-whole-brain-TMean-',num2str(pvalue),'.nii');

input_data = MOT2RestingStateDiff; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT4-MOT2-Diff','-','zscore-whole-brain-TMean-',num2str(pvalue),'.nii');

input_data = MOT4MOT2Diff; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT4-RestingState-Diff-Positive','-','zscore-whole-brain-TMean-',num2str(pvalue),'.nii');

input_data = MOT4RestingStateDiff_positive; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT2-RestingState-Diff-Positive','-','zscore-whole-brain-TMean-',num2str(pvalue),'.nii');

input_data = MOT2RestingStateDiff_positive; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT4-MOT2-Diff-Positive','-','zscore-whole-brain-TMean-',num2str(pvalue),'.nii');

input_data = MOT4MOT2Diff_positive; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT4-RestingState-Diff-Negative','-','zscore-whole-brain-TMean-',num2str(pvalue),'.nii');

input_data = MOT4RestingStateDiff_negative; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT2-RestingState-Diff-Negative','-','zscore-whole-brain-TMean-',num2str(pvalue),'.nii');

input_data = MOT2RestingStateDiff_negative; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT4-MOT2-Diff-Negative','-','zscore-whole-brain-TMean-',num2str(pvalue),'.nii');

input_data = MOT4MOT2Diff_negative; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT4-RestingState-Diff-Negative-Inverted','-','zscore-whole-brain-TMean-',num2str(pvalue),'.nii');

input_data = MOT4RestingStateDiff_negative_inverted; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT2-RestingState-Diff-Negative-Inverted','-','zscore-whole-brain-TMean-',num2str(pvalue),'.nii');

input_data = MOT2RestingStateDiff_negative_inverted; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT4-MOT2-Diff-Negative-Inverted','-','zscore-whole-brain-TMean-',num2str(pvalue),'.nii');

input_data = MOT4MOT2Diff_negative_inverted; 
input_data = double(input_data)*10^(-scale);

lowhigh_save_image;

return

function plotHistogram(MOT4_Threshold,MOT2_Threshold,RestingState_Threshold,MOT4RestingStateDiff,MOT2RestingStateDiff,MOT4MOT2Diff,pvalue,kind,common_mask,common_idx,settings,scale)

if strcmpi('threshold',kind)
    
    %mask_vec = reshape(common_mask,[1,size(common_mask,1)*size(common_mask,2)*size(common_mask,3)]);
    %all_idx = find(common_mask<2);
    %all_idx(common_idx) = [];
    idx_brain = common_idx;
    
else
    
    %mask_vec = reshape(common_mask,[1,size(common_mask,1)*size(common_mask,2)*size(common_mask,3)]);
    idx_brain = find(common_mask ~= 0);

end

%MOT4_threshold_vec = reshape(MOT4_Threshold,[1,size(MOT4_Threshold,1)*size(MOT4_Threshold,2)*size(MOT4_Threshold,3)]);
% MOT4_threshold_vec(idx_space) = [];
MOT4_threshold_vec = double(MOT4_Threshold(idx_brain))*10^(-scale);
%
%MOT2_threshold_vec = reshape(MOT2_Threshold,[1,size(MOT2_Threshold,1)*size(MOT2_Threshold,2)*size(MOT2_Threshold,3)]);
% MOT2_threshold_vec(idx_space) = [];
MOT2_threshold_vec = double(MOT2_Threshold(idx_brain))*10^(-scale);
% 
%RestingState_threshold_vec = reshape(RestingState_Threshold,[1,size(RestingState_Threshold,1)*size(RestingState_Threshold,2)*size(RestingState_Threshold,3)]);
% RestingState_threshold_vec(idx_space) = [];
RestingState_threshold_vec = double(RestingState_Threshold(idx_brain))*10^(-scale);
% 
%MOT4RestingStateDiff_vec = reshape(MOT4RestingStateDiff,[1,size(MOT4RestingStateDiff,1)*size(MOT4RestingStateDiff,2)*size(MOT4RestingStateDiff,3)]);
% MOT4RestingStateDiff_vec(idx_space) = [];
MOT4RestingStateDiff_vec = double(MOT4RestingStateDiff(idx_brain))*10^(-scale);
% 
%MOT2RestingStateDiff_vec = reshape(MOT2RestingStateDiff,[1,size(MOT2RestingStateDiff,1)*size(MOT2RestingStateDiff,2)*size(MOT2RestingStateDiff,3)]);
% MOT2RestingStateDiff_vec(idx_space) = [];
MOT2RestingStateDiff_vec = double(MOT2RestingStateDiff(idx_brain))*10^(-scale);
% 
%MOT4MOT2Diff_vec = reshape(MOT4MOT2Diff,[1,size(MOT4MOT2Diff,1)*size(MOT4MOT2Diff,2)*size(MOT4MOT2Diff,3)]);
% MOT4MOT2Diff_vec(idx_space) = [];
MOT4MOT2Diff_vec = double(MOT4MOT2Diff(idx_brain))*10^(-scale);

MOT4_RestingState_covariance = getCovariance(MOT4_threshold_vec,RestingState_threshold_vec);
MOT2_RestingState_covariance = getCovariance(MOT2_threshold_vec,RestingState_threshold_vec);
MOT4_MOT2_covariance = getCovariance(MOT4_threshold_vec,MOT2_threshold_vec);

max_threshold = ceil(max([max(MOT4_threshold_vec),max(MOT2_threshold_vec),max(RestingState_threshold_vec)]));
min_threshold = floor(min([min(MOT4_threshold_vec),min(MOT2_threshold_vec),min(RestingState_threshold_vec)]));

max_diff = ceil(max([max(MOT4RestingStateDiff_vec),max(MOT2RestingStateDiff_vec),max(MOT4MOT2Diff_vec)]));
min_diff = floor(min([min(MOT4RestingStateDiff_vec),min(MOT2RestingStateDiff_vec),min(MOT4MOT2Diff_vec)]));

nPoints = 1000;

bins_centers_threshold = linspace(min_threshold,max_threshold,nPoints);
bins_centers_diff = linspace(min_diff,max_diff,nPoints);

f = figure;
hold on;
fs = 9;

%%% PLOT Z-SCORE THRESHOLDED

[nquality,xquality] = hist(MOT4_threshold_vec,bins_centers_threshold);
fquality = nquality / sum(nquality);
dquality_MOT4 = fquality;
xquality_MOT4 = xquality;
%dquality = fquality / (xquality(2) - xquality(1));
mean_x_MOT4_threshold_vec = double(mean(MOT4_Threshold(idx_brain),'native'))*10^(-scale);
X_MOT4 = mean_x_MOT4_threshold_vec;
Y_MOT4 = max(dquality_MOT4) / 2;

[nquality,xquality] = hist(MOT2_threshold_vec,bins_centers_threshold);
fquality = nquality / sum(nquality);
dquality_MOT2 = fquality;
xquality_MOT2 = xquality;
%dquality = fquality / (xquality(2) - xquality(1));
mean_x_MOT2_threshold_vec = double(mean(MOT2_Threshold(idx_brain),'native'))*10^(-scale);
X_MOT2 = mean_x_MOT2_threshold_vec;
Y_MOT2 = max(dquality_MOT2) / 2;

[nquality,xquality] = hist(RestingState_threshold_vec,bins_centers_threshold);
fquality = nquality / sum(nquality);
dquality_RestingState = fquality;
xquality_RestingState = xquality;
%dquality = fquality / (xquality(2) - xquality(1));
mean_x_RestingState_threshold_vec = double(mean(RestingState_Threshold(idx_brain),'native'))*10^(-scale);
X_RestingState = mean_x_RestingState_threshold_vec;
Y_RestingState = max(dquality_RestingState) / 2;

min_dquality = min([dquality_MOT4, dquality_MOT2, dquality_RestingState]);
max_dquality = max([dquality_MOT4, dquality_MOT2, dquality_RestingState]);

subplot(2,3,1);
bar( xquality_MOT4, dquality_MOT4  ); 
xlabel('ZScore','FontSize',fs);
ylabel('density','FontSize',fs);
text(min_threshold,Y_MOT4,strcat('mean:',num2str(mean_x_MOT4_threshold_vec)));
title('MOT 4');
xlim([min_threshold max_threshold]);
ylim([min_dquality max_dquality]);

subplot(2,3,2);
bar( xquality_MOT2, dquality_MOT2  ); 
xlabel('ZScore','FontSize',fs);
ylabel('density','FontSize',fs);
text(min_threshold,Y_MOT2,strcat('mean:',num2str(mean_x_MOT2_threshold_vec)));
title('MOT 2');
xlim([min_threshold max_threshold]);
ylim([min_dquality max_dquality]);

subplot(2,3,3);
bar( xquality_RestingState, dquality_RestingState  ); 
xlabel('ZScore','FontSize',fs);
ylabel('density','FontSize',fs);
text(min_threshold,Y_RestingState,strcat('mean:',num2str(mean_x_RestingState_threshold_vec)));
title('Resting State');
xlim([min_threshold max_threshold]);
ylim([min_dquality max_dquality]);

%%% PLOT DIFFERENCE

[nquality,xquality] = hist(MOT4RestingStateDiff_vec,bins_centers_diff);
fquality = nquality / sum(nquality);
dquality_MOT4RestingState = fquality;
xquality_MOT4RestingState = xquality;
%dquality = fquality / (xquality(2) - xquality(1));
mean_x_MOT4RestingState = double(mean(MOT4RestingStateDiff(idx_brain),'native'))*10^(-scale);
X_MOT4RestingState = mean_x_MOT4RestingState;
Y_MOT4RestingState = max(dquality_MOT4RestingState) / 2;

[nquality,xquality] = hist(MOT2RestingStateDiff_vec,bins_centers_diff);
fquality = nquality / sum(nquality);
dquality_MOT2RestingState = fquality;
xquality_MOT2RestingState = xquality;
%dquality = fquality / (xquality(2) - xquality(1));
mean_x_MOT2RestingState = double(mean(MOT2RestingStateDiff(idx_brain),'native'))*10^(-scale);
X_MOT2RestingState = mean_x_MOT2RestingState;
Y_MOT2RestingState = max(dquality_MOT2RestingState) / 2;

[nquality,xquality] = hist(MOT4MOT2Diff_vec,bins_centers_diff);
fquality = nquality / sum(nquality);
dquality_MOT4MOT2 = fquality;
xquality_MOT4MOT2 = xquality;
%dquality = fquality / (xquality(2) - xquality(1));
mean_x_MOT4MOT2 = double(mean(MOT4MOT2Diff(idx_brain),'native'))*10^(-scale);
X_MOT4MOT2 = mean_x_MOT4MOT2;
Y_MOT4MOT2 = max(dquality_MOT4MOT2) / 2;

min_dquality = min([dquality_MOT4RestingState, dquality_MOT2RestingState, dquality_MOT4MOT2]);
max_dquality = max([dquality_MOT4RestingState, dquality_MOT2RestingState, dquality_MOT4MOT2]);

subplot(2,3,4);
bar( xquality_MOT4RestingState, dquality_MOT4RestingState  ); 
xlabel('ZScore Difference','FontSize',fs);
ylabel('density','FontSize',fs);
title('MOT 4 - Resting State');
text(min_diff,Y_MOT4RestingState,strcat('mean:',num2str(mean_x_MOT4RestingState)));
text(min_diff,Y_MOT4RestingState/2,strcat('cov:',num2str(MOT4_RestingState_covariance)));
xlim([min_diff max_diff]);
ylim([min_dquality max_dquality]);

subplot(2,3,5);
bar( xquality_MOT2RestingState, dquality_MOT2RestingState  ); 
xlabel('ZScore Difference','FontSize',fs);
ylabel('density','FontSize',fs);
title('MOT 2 - Resting State');
text(min_diff,Y_MOT2RestingState,strcat('mean:',num2str(mean_x_MOT2RestingState)));
text(min_diff,Y_MOT2RestingState/2,strcat('cov:',num2str(MOT2_RestingState_covariance)));
xlim([min_diff max_diff]);
ylim([min_dquality max_dquality]);

subplot(2,3,6);
bar( xquality_MOT4MOT2, dquality_MOT4MOT2  ); 
xlabel('ZScore Difference','FontSize',fs);
ylabel('density','FontSize',fs);
title('MOT 4 - MOT 2');
text(min_diff,Y_MOT4MOT2,strcat('mean:',num2str(mean_x_MOT4MOT2)));
text(min_diff,Y_MOT4MOT2/2,strcat('cov:',num2str(MOT4_MOT2_covariance)));
xlim([min_diff max_diff]);
ylim([min_dquality max_dquality]);

if strcmpi('threshold',kind)
    
    suptitle(strcat('p-value:',num2str(pvalue)));

else
    
    suptitle('No Threshold');
    
end

print(f,'-djpeg',strcat('Low-High-',settings.folders.subject,'-','ZScore-Diff-Histogram-AllRuns-p-value-',num2str(pvalue),'.jpeg'));
print(f,'-depsc',strcat('Low-High-',settings.folders.subject,'-','ZScore-Diff-Histogram-AllRuns-p-value-',num2str(pvalue),'.eps'));

return

function covariance = getCovariance(X,Y)

diagonal = diag(flip(cov(X,Y)));
covariance = diagonal(1);

return

function saveStdThreshold(MOT4_Std_Thresholded,MOT2_Std_Thresholded,RestingState_Std_Thresholded,pvalue,settings,kind)

run = 1;
nifti_file = settings.restingstate.FSL.run(run).nifti_file;
dtype = settings.restingstate.FSL.run(run).dtype;
offset = settings.restingstate.FSL.run(run).offset;
scl_slope = settings.restingstate.FSL.run(run).scl_slope;
scl_inter = settings.restingstate.FSL.run(run).scl_inter;
    
dim = [settings.restingstate.FSL.run(run).dim(1), settings.restingstate.FSL.run(run).dim(2),settings.restingstate.FSL.run(run).dim(3)];

dtype = 'FLOAT32';
offset = 0;

descrip = 'std ZScore Thresholded';

fname = strcat('Low-High-',settings.folders.subject,'-MOT4-AllRuns-','zscore-whole-brain-Tstd-',kind,num2str(pvalue),'.nii');

input_data = MOT4_Std_Thresholded; 

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-MOT2-AllRuns-','zscore-whole-brain-Tstd-',kind,num2str(pvalue),'.nii');

input_data = MOT2_Std_Thresholded; 

lowhigh_save_image;

fname = strcat('Low-High-',settings.folders.subject,'-RestingState-AllRuns-','zscore-whole-brain-Tstd-',kind,num2str(pvalue),'.nii');

input_data = RestingState_Std_Thresholded; 

lowhigh_save_image;

return