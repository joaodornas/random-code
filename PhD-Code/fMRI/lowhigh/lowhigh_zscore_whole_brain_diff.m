
function lowhigh_zscore_whole_brain_diff

settings_jan_0502;
doTheMath(settings);

return

function doTheMath(settings)

k = 3; %% polynomial order
frame_size = 120 - 1;

%% REALIGNED
% get_at_this_preprocessed_step = settings.folders.realigned.name;
% prefix_for_the_preprocessed_step = settings.folders.realigned.prefix;

%% WARPED
get_at_this_preprocessed_step = settings.folders.normalization.name;
prefix_for_the_preprocessed_step = settings.folders.normalization.prefix;

lowhigh_load_all_data;

descrip = 'mean Z-Score';

flag_zscore = 1;
flag_std = 0;

nBins = 21;

disp('MOT4 - zscore');

nTR = size(MOT4Run1,4);

MOT4Run1sgolay = sgolayfilt(MOT4Run1,k,frame_size,[],4);
MOT4Run2sgolay = sgolayfilt(MOT4Run2,k,frame_size,[],4);

MOT4Run1removed = MOT4Run1 - MOT4Run1sgolay;
MOT4Run2removed = MOT4Run2 - MOT4Run2sgolay;
 
MOT4Run1ZSt = zscore_whole_brain(MOT4Run1removed);
MOT4Run2ZSt = zscore_whole_brain(MOT4Run2removed);

MOT4Run1ZStmean = mean(MOT4Run1ZSt,4);
MOT4Run2ZStmean = mean(MOT4Run2ZSt,4);

MOT4All(:,:,:,1) = MOT4Run1ZStmean(:,:,:);
MOT4All(:,:,:,2) = MOT4Run2ZStmean(:,:,:);

MOT4mean = mean(MOT4All,4);

for run=1:2;
    
    nifti_file = settings.mot4.run(run).nifti_file;
    dtype = settings.mot4.run(run).dtype;
    offset = settings.mot4.run(run).offset;
    scl_slope = settings.mot4.run(run).scl_slope;
    scl_inter = settings.mot4.run(run).scl_inter;
    
    dtype = 'FLOAT32';
    offset = 0;
    
    if run == 1; 
        
        fname = strcat('Low-High-',settings.subject,'-MOT4-Run-',int2str(run),'-','zscore-whole-brain-t','.nii');
        dim = [settings.mot4.run(run).dim, nTR];
        input_data = MOT4Run1ZSt; 
        
        lowhigh_save_image;
        
        fname = strcat('Low-High-',settings.subject,'-MOT4-Run-',int2str(run),'-','zscore-whole-brain-mean-t','.nii');
        dim = settings.mot4.run(run).dim;
        input_data = MOT4Run1ZStmean; 
        
        lowhigh_save_image;
        
%         fname = strcat('Low-High-',settings.subject,'-MOT4-Run-',int2str(run),'-','zscore-whole-brain-std-t','.nii');
%         dim = settings.mot4.run(run).dim;
%         input_data = MOT4Run1ZStstd;
%         
%         lowhigh_save_image;
    
    end
    
    if run == 2; 
        
        fname = strcat('Low-High-',settings.subject,'-MOT4-Run-',int2str(run),'-','zscore-whole-brain-t','.nii');
        dim = [settings.mot4.run(run).dim, nTR];
        input_data = MOT4Run2ZSt; 
        
        lowhigh_save_image;
        
        fname = strcat('Low-High-',settings.subject,'-MOT4-Run-',int2str(run),'-','zscore-whole-brain-mean-t','.nii');
        dim = settings.mot4.run(run).dim;
        input_data = MOT4Run2ZStmean; 
        
        lowhigh_save_image;

%         fname = strcat('Low-High-',settings.subject,'-MOT4-Run-',int2str(run),'-','zscore-whole-brain-std-t','.nii');
%         dim = settings.mot4.run(run).dim;
%         input_data = MOT4Run2ZStstd;
%         
%         lowhigh_save_image;

    end

end

disp('MOT2 - zscore');

nTR = size(MOT2Run1,4);

MOT2Run1sgolay = sgolayfilt(MOT2Run1,k,frame_size,[],4);
MOT2Run2sgolay = sgolayfilt(MOT2Run2,k,frame_size,[],4);

MOT2Run1removed = MOT2Run1 - MOT2Run1sgolay;
MOT2Run2removed = MOT2Run2 - MOT2Run2sgolay;

MOT2Run1ZSt = zscore_whole_brain(MOT2Run1removed);
MOT2Run2ZSt = zscore_whole_brain(MOT2Run2removed);

MOT2Run1ZStmean = mean(MOT2Run1ZSt,4);
MOT2Run2ZStmean = mean(MOT2Run2ZSt,4);

MOT2All(:,:,:,1) = MOT2Run1ZStmean(:,:,:);
MOT2All(:,:,:,2) = MOT2Run2ZStmean(:,:,:);

MOT2mean = mean(MOT2All,4);

for run=1:2;
    
    nifti_file = settings.mot2.run(run).nifti_file;
    dtype = settings.mot2.run(run).dtype;
    offset = settings.mot2.run(run).offset;
    scl_slope = settings.mot2.run(run).scl_slope;
    scl_inter = settings.mot2.run(run).scl_inter;
    
    dtype = 'FLOAT32';
    offset = 0;
    
    if run == 1; 
        
        fname = strcat('Low-High-',settings.subject,'-MOT2-Run-',int2str(run),'-','zscore-whole-brain-t','.nii');
        dim = [settings.mot2.run(run).dim, nTR];
        input_data = MOT2Run1ZSt; 
        
        lowhigh_save_image;
        
        fname = strcat('Low-High-',settings.subject,'-MOT2-Run-',int2str(run),'-','zscore-whole-brain-mean-t','.nii');
        dim = settings.mot2.run(run).dim;
        input_data = MOT2Run1ZStmean; 
        
        lowhigh_save_image;
        
%         fname = strcat('Low-High-',settings.subject,'-MOT2-Run-',int2str(run),'-','zscore-whole-brain-std-t','.nii');
%         dim = settings.mot2.run(run).dim;
%         input_data = MOT2Run1ZStstd;
%         
%         lowhigh_save_image;
    
    end
    
    if run == 2; 
        
        fname = strcat('Low-High-',settings.subject,'-MOT2-Run-',int2str(run),'-','zscore-whole-brain-t','.nii');
        dim = [settings.mot2.run(run).dim, nTR];
        input_data = MOT2Run2ZSt; 
        
        lowhigh_save_image;
        
        fname = strcat('Low-High-',settings.subject,'-MOT2-Run-',int2str(run),'-','zscore-whole-brain-mean-t','.nii');
        dim = settings.mot2.run(run).dim;
        input_data = MOT2Run2ZStmean; 
        
        lowhigh_save_image;
        
%         fname = strcat('Low-High-',settings.subject,'-MOT2-Run-',int2str(run),'-','zscore-whole-brain-std-t','.nii');
%         dim = settings.mot2.run(run).dim;
%         input_data = MOT2Run2ZStstd;
%         
%         lowhigh_save_image;
    
    end

end

disp('RestingState - zscore');

nTR = size(RestingStateRun1,4);

RestingStateRun1sgolay = sgolayfilt(RestingStateRun1,k,frame_size,[],4);
RestingStateRun2sgolay = sgolayfilt(RestingStateRun2,k,frame_size,[],4);

RestingStateRun1removed = RestingStateRun1 - RestingStateRun1sgolay;
RestingStateRun2removed = RestingStateRun2 - RestingStateRun2sgolay;

RestingStateRun1ZSt = zscore_whole_brain(RestingStateRun1removed);
RestingStateRun2ZSt = zscore_whole_brain(RestingStateRun2removed);

RestingStateRun1ZStmean = mean(RestingStateRun1ZSt,4);
RestingStateRun2ZStmean = mean(RestingStateRun2ZSt,4);

RestingStateAll(:,:,:,1) = RestingStateRun1ZStmean(:,:,:);
RestingStateAll(:,:,:,2) = RestingStateRun2ZStmean(:,:,:);

RestingStatemean = mean(RestingStateAll,4);

disp('Save data');

for run=1:2;
    
    nifti_file = settings.restingstate.run(run).nifti_file;
    dtype = settings.restingstate.run(run).dtype;
    offset = settings.restingstate.run(run).offset;
    scl_slope = settings.restingstate.run(run).scl_slope;
    scl_inter = settings.restingstate.run(run).scl_inter;
    
    dtype = 'FLOAT32';
    offset = 0;
    
    if run == 1; 
        
        fname = strcat('Low-High-',settings.subject,'-RestingState-Run-',int2str(run),'-','zscore-whole-brain-t','.nii');
        dim = [settings.restingstate.run(run).dim, nTR];
        input_data = RestingStateRun1ZSt; 
        
        lowhigh_save_image;
        
        fname = strcat('Low-High-',settings.subject,'-RestingState-Run-',int2str(run),'-','zscore-whole-brain-mean-t','.nii');
        dim = settings.restingstate.run(run).dim;
        input_data = RestingStateRun1ZStmean; 
        
        lowhigh_save_image;
        
%         fname = strcat('Low-High-',settings.subject,'-RestingState-Run-',int2str(run),'-','zscore-whole-brain-std-t','.nii');
%         dim = settings.restingstate.run(run).dim;
%         input_data = RestingStateRun1ZStstd;
%         
%         lowhigh_save_image;
    
    end
    
    if run == 2; 
        
        fname = strcat('Low-High-',settings.subject,'-RestingState-Run-',int2str(run),'-','zscore-whole-brain-t','.nii');
        dim = [settings.restingstate.run(run).dim, nTR];
        input_data = RestingStateRun2ZSt; 
        
        lowhigh_save_image;
        
        fname = strcat('Low-High-',settings.subject,'-RestingState-Run-',int2str(run),'-','zscore-whole-brain-mean-t','.nii');
        dim = settings.restingstate.run(run).dim;
        input_data = RestingStateRun2ZStmean; 
        
        lowhigh_save_image;
        
%         fname = strcat('Low-High-',settings.subject,'-RestingState-Run-',int2str(run),'-','zscore-whole-brain-std-t','.nii');
%         dim = settings.restingstate.run(run).dim;
%         input_data = RestingStateRun2ZStstd;
%         
%         lowhigh_save_image;
    
    end

end

disp('Difference');

MOT4RestingStateDiff = MOT4mean - RestingStatemean;
MOT2RestingStateDiff = MOT2mean - RestingStatemean;
MOT4MOT2Diff = MOT4mean - MOT2mean;

disp('Save Data');

run = 1;
nifti_file = settings.restingstate.run(run).nifti_file;
dtype = settings.restingstate.run(run).dtype;
offset = settings.restingstate.run(run).offset;
scl_slope = settings.restingstate.run(run).scl_slope;
scl_inter = settings.restingstate.run(run).scl_inter;
    
dtype = 'FLOAT32';
offset = 0;

fname = strcat('Low-High-',settings.subject,'-MOT4-RestingState-Diff','-','zscore-whole-brain-t','.nii');
dim = settings.restingstate.run(run).dim;
input_data = MOT4RestingStateDiff; 

lowhigh_save_image;

fname = strcat('Low-High-',settings.subject,'-MOT2-RestingState-Diff','-','zscore-whole-brain-t','.nii');
dim = settings.restingstate.run(run).dim;
input_data = MOT2RestingStateDiff; 

lowhigh_save_image;

fname = strcat('Low-High-',settings.subject,'-MOT4-MOT2-Diff','-','zscore-whole-brain-t','.nii');
dim = settings.restingstate.run(run).dim;
input_data = MOT4MOT2Diff; 

lowhigh_save_image;

return

