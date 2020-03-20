function lowhigh_ttest_ROI

settings_jan_0502;
doTheMath(settings);

end

function doTheMath(settings)

%% WARPED
get_at_this_preprocessed_step = settings.folders.smooth.name;
prefix_for_the_preprocessed_step = settings.folders.smooth.prefix;

lowhigh_load_all_data;

descrip = 'p-value';

summarization = 'mean';

k = 3; %% polynomial order
frame_size = 120 - 1;

alpha = 0.05;

% disp('low-drift');
% 
% MOT4Run1sgolay = sgolayfilt(MOT4Run1,k,frame_size,[],4);
% MOT4Run2sgolay = sgolayfilt(MOT4Run2,k,frame_size,[],4);
% 
% MOT4Run1 = MOT4Run1 - MOT4Run1sgolay;
% MOT4Run2 = MOT4Run2 - MOT4Run2sgolay;

% 
% MOT2Run1sgolay = sgolayfilt(MOT2Run1,k,frame_size,[],4);
% MOT2Run2sgolay = sgolayfilt(MOT2Run2,k,frame_size,[],4);
% 
% MOT2Run1 = MOT2Run1 - MOT2Run1sgolay;
% MOT2Run2 = MOT2Run2 - MOT2Run2sgolay;
% 
% RestingStateRun1sgolay = sgolayfilt(RestingStateRun1,k,frame_size,[],4);
% RestingStateRun2sgolay = sgolayfilt(RestingStateRun2,k,frame_size,[],4);
% 
% RestingStateRun1 = RestingStateRun1 - RestingStateRun1sgolay;
% RestingStateRun2 = RestingStateRun2 - RestingStateRun2sgolay;

disp('MOT4 - ROIs');

[ MOT4Run1_ROIs, labels ] = get_ROI_aal( MOT4Run1, summarization );
[ MOT4Run2_ROIs, labels ] = get_ROI_aal( MOT4Run2, summarization );


disp('MOT2 - ROIs');

[ MOT2Run1_ROIs, labels ] = get_ROI_aal( MOT2Run1, summarization );
[ MOT2Run2_ROIs, labels ] = get_ROI_aal( MOT2Run2, summarization );

disp('RestingState - ROIs');

[ RestingStateRun1_ROIs, labels ] = get_ROI_aal( RestingStateRun1, summarization );
[ RestingStateRun2_ROIs, labels ] = get_ROI_aal( RestingStateRun2, summarization );

disp('t-test');

[H_MOT4Run1_RestingStateRun1, P_MOT4Run1_RestingStateRun1, label_MOT4Run1_RestingStateRun1] = tcontrast(MOT4Run1_ROIs,RestingStateRun1_ROIs,labels,alpha);
[H_MOT4Run1_RestingStateRun2, P_MOT4Run1_RestingStateRun2, label_MOT4Run1_RestingStateRun2] = tcontrast(MOT4Run1_ROIs,RestingStateRun2_ROIs,labels,alpha);
[H_MOT4Run2_RestingStateRun1, P_MOT4Run2_RestingStateRun1, label_MOT4Run2_RestingStateRun1] = tcontrast(MOT4Run2_ROIs,RestingStateRun1_ROIs,labels,alpha);
[H_MOT4Run2_RestingStateRun2, P_MOT4Run2_RestingStateRun2, label_MOT4Run2_RestingStateRun2] = tcontrast(MOT4Run2_ROIs,RestingStateRun2_ROIs,labels,alpha);

[H_MOT2Run1_RestingStateRun1, P_MOT2Run1_RestingStateRun1, label_MOT2Run1_RestingStateRun1] = tcontrast(MOT2Run1_ROIs,RestingStateRun1_ROIs,labels,alpha);
[H_MOT2Run1_RestingStateRun2, P_MOT2Run1_RestingStateRun2, label_MOT2Run1_RestingStateRun2] = tcontrast(MOT2Run1_ROIs,RestingStateRun2_ROIs,labels,alpha);
[H_MOT2Run2_RestingStateRun1, P_MOT2Run2_RestingStateRun1, label_MOT2Run2_RestingStateRun1] = tcontrast(MOT2Run2_ROIs,RestingStateRun1_ROIs,labels,alpha);
[H_MOT2Run2_RestingStateRun2, P_MOT2Run2_RestingStateRun2, label_MOT2Run2_RestingStateRun2] = tcontrast(MOT2Run2_ROIs,RestingStateRun2_ROIs,labels,alpha);

[H_MOT4Run1_MOT2Run1, P_MOT4Run1_MOT2Run1, label_MOT4Run1_MOT2Run1] = tcontrast(MOT4Run1_ROIs,MOT2Run1_ROIs,labels,alpha);
[H_MOT4Run1_MOT2Run2, P_MOT4Run1_MOT2Run2, label_MOT4Run1_MOT2Run2] = tcontrast(MOT4Run1_ROIs,MOT2Run2_ROIs,labels,alpha);
[H_MOT4Run2_MOT2Run1, P_MOT4Run2_MOT2Run1, label_MOT4Run2_MOT2Run1] = tcontrast(MOT4Run2_ROIs,MOT2Run1_ROIs,labels,alpha);
[H_MOT4Run2_MOT2Run2, P_MOT4Run2_MOT2Run2, label_MOT4Run2_MOT2Run2] = tcontrast(MOT4Run2_ROIs,MOT2Run2_ROIs,labels,alpha);

[H_MOT4Run1_MOT4Run2, P_MOT4Run1_MOT4Run2, label_MOT4Run1_MOT4Run2] = tcontrast(MOT4Run1_ROIs,MOT4Run2_ROIs,labels,alpha);
[H_MOT2Run1_MOT2Run2, P_MOT2Run1_MOT2Run2, label_MOT2Run1_MOT2Run2] = tcontrast(MOT2Run1_ROIs,MOT2Run2_ROIs,labels,alpha);
[H_RestingStateRun1_RestingStateRun2, P_RestingStateRun1_RestingStateRun2, label_RestingStateRun1_RestingStateRun2] = tcontrast(RestingStateRun1_ROIs,RestingStateRun2_ROIs,labels,alpha);

MOT4_ROIs = [MOT4Run1_ROIs;MOT4Run2_ROIs];
MOT2_ROIs = [MOT2Run1_ROIs;MOT2Run2_ROIs];
RestingState_ROIs = [RestingStateRun1_ROIs;RestingStateRun2_ROIs];

[H_MOT4_MOT2,P_MOT4_MOT2,label_MOT4_MOT2] = tcontrast(MOT4_ROIs,MOT2_ROIs,labels,alpha);
[H_MOT4_RestingState,P_MOT4_RestingState,label_MOT4_RestingState] = tcontrast(MOT4_ROIs,RestingState_ROIs,labels,alpha);
[H_MOT2_RestingState,P_MOT2_RestingState,label_MOT2_RestingState] = tcontrast(MOT2_ROIs,RestingState_ROIs,labels,alpha);

name = strcat('Low-High-',settings.subject,'-','ttest');

disp('Save Everything');

save(name,'P_MOT4Run1_RestingStateRun1','P_MOT4Run1_RestingStateRun2',...
    'P_MOT4Run2_RestingStateRun1','P_MOT4Run2_RestingStateRun2',...
    'P_MOT2Run1_RestingStateRun1','P_MOT2Run1_RestingStateRun2',...
    'P_MOT2Run2_RestingStateRun1','P_MOT2Run2_RestingStateRun2',...
    'P_MOT4Run1_MOT2Run1','P_MOT4Run1_MOT2Run2',...
    'P_MOT4Run2_MOT2Run1','P_MOT4Run2_MOT2Run2',...
    'P_MOT4Run1_MOT4Run2','P_MOT4Run1_MOT4Run2','P_RestingStateRun1_RestingStateRun2',...
    'P_MOT4_MOT2','P_MOT4_RestingState','P_MOT2_RestingState',...
    'label_MOT4Run1_RestingStateRun1','label_MOT4Run1_RestingStateRun2',...
    'label_MOT4Run2_RestingStateRun1','label_MOT4Run2_RestingStateRun2',...
    'label_MOT2Run1_RestingStateRun1','label_MOT2Run1_RestingStateRun2',...
    'label_MOT2Run2_RestingStateRun1','label_MOT2Run2_RestingStateRun2',...
    'label_MOT4Run1_MOT2Run1','label_MOT4Run1_MOT2Run2',...
    'label_MOT4Run2_MOT2Run1','label_MOT4Run2_MOT2Run2',...
    'label_MOT4Run1_MOT4Run2','label_MOT2Run1_MOT2Run2','label_RestingStateRun1_RestingStateRun2',...
    'label_MOT4_MOT2','label_MOT4_RestingState','label_MOT2_RestingState',...
    'labels');

run = 1;

nifti_file = settings.mot4.run(run).nifti_file;
dtype = settings.mot4.run(run).dtype;
offset = settings.mot4.run(run).offset;
scl_slope = settings.mot4.run(run).scl_slope;
scl_inter = settings.mot4.run(run).scl_inter;

dtype = 'FLOAT32';
offset = 0;

fname = strcat('Low-High-',settings.subject,'-MOT4-Run1-','RestingState-Run1','-','t-test-p-value','.nii');
ROIs = P_MOT4Run1_RestingStateRun1;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);

fname = strcat('Low-High-',settings.subject,'-MOT4-Run1-','RestingState-Run2','-','t-test-p-value','.nii');
ROIs = P_MOT4Run1_RestingStateRun2;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);

fname = strcat('Low-High-',settings.subject,'-MOT4-Run2-','RestingState-Run1','-','t-test-p-value','.nii');
ROIs = P_MOT4Run2_RestingStateRun1;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);

fname = strcat('Low-High-',settings.subject,'-MOT4-Run2-','RestingState-Run2','-','t-test-p-value','.nii');
ROIs = P_MOT4Run2_RestingStateRun2;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);


fname = strcat('Low-High-',settings.subject,'-MOT2-Run1-','RestingState-Run1','-','t-test-p-value','.nii');
ROIs = P_MOT2Run1_RestingStateRun1;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);

fname = strcat('Low-High-',settings.subject,'-MOT2-Run1-','RestingState-Run2','-','t-test-p-value','.nii');
ROIs = P_MOT2Run1_RestingStateRun2;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);

fname = strcat('Low-High-',settings.subject,'-MOT2-Run2-','RestingState-Run1','-','t-test-p-value','.nii');
ROIs = P_MOT2Run2_RestingStateRun1;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);

fname = strcat('Low-High-',settings.subject,'-MOT2-Run2-','RestingState-Run2','-','t-test-p-value','.nii');
ROIs = P_MOT2Run2_RestingStateRun2;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);


fname = strcat('Low-High-',settings.subject,'-MOT4-Run1-','MOT2-Run1','-','t-test-p-value','.nii');
ROIs = P_MOT4Run1_MOT2Run1;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);

fname = strcat('Low-High-',settings.subject,'-MOT4-Run1-','MOT2-Run2','-','t-test-p-value','.nii');
ROIs = P_MOT4Run1_MOT2Run2;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);

fname = strcat('Low-High-',settings.subject,'-MOT4-Run2-','MOT2-Run1','-','t-test-p-value','.nii');
ROIs = P_MOT4Run2_MOT2Run1;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);

fname = strcat('Low-High-',settings.subject,'-MOT4-Run2-','MOT2-Run2','-','t-test-p-value','.nii');
ROIs = P_MOT4Run2_MOT2Run2;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);


fname = strcat('Low-High-',settings.subject,'-MOT4-Run1-','MOT4-Run2','-','t-test-p-value','.nii');
ROIs = P_MOT4Run1_MOT4Run2;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);

fname = strcat('Low-High-',settings.subject,'-MOT2-Run1-','MOT2-Run2','-','t-test-p-value','.nii');
ROIs = P_MOT2Run1_MOT2Run2;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);

fname = strcat('Low-High-',settings.subject,'-RestingState-Run1-','RestingState-Run2','-','t-test-p-value','.nii');
ROIs = P_RestingStateRun1_RestingStateRun2;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);
  
fname = strcat('Low-High-',settings.subject,'-MOT4-','MOT2','-','t-test-p-value','.nii');
ROIs = P_MOT4_MOT2;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);

fname = strcat('Low-High-',settings.subject,'-MOT4-','RestingState','-','t-test-p-value','.nii');
ROIs = P_MOT4_RestingState;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);

fname = strcat('Low-High-',settings.subject,'-MOT2-','RestingState','-','t-test-p-value','.nii');
ROIs = P_MOT2_RestingState;
save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip);

end

function [H, P, threshold_labels] = tcontrast(ROIs1, ROIs2, labels, alpha)

contrast = ROIs1 - ROIs2;

[H, P] = ttest(contrast);

idx = find(P < alpha);

threshold_labels = labels(idx);

end

