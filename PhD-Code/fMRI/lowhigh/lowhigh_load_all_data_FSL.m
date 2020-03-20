

%file = settings.FSL.files.functional.warped;
%file = settings.FSL.files.functional.masked;
%file = settings.FSL.files.functional.despike;

%file = settings.FSL.files.functional.custom.residual_voxel;
%mask = settings.FSL.files.mask.custom;

%mask = settings.FSL.files.mask.warped;
%mask = settings.FSL.files.mask.custom;

%get_at_this_preprocessed_step = settings.FSL.folders.warped;
%get_at_this_preprocessed_step = settings.FSL.folders.custom;

disp(file);

kind = 'MOT4';
run = 1;
[MOT4Run1, mask_MOT4Run1, settings] = lowhigh_get_data_FSL(settings,kind,run,file,mask,get_at_this_preprocessed_step);
run = 2;
[MOT4Run2, mask_MOT4Run2, settings] = lowhigh_get_data_FSL(settings,kind,run,file,mask,get_at_this_preprocessed_step);

kind = 'MOT2';
run = 1;
[MOT2Run1, mask_MOT2Run1, settings] = lowhigh_get_data_FSL(settings,kind,run,file,mask,get_at_this_preprocessed_step);
run = 2;
[MOT2Run2, mask_MOT2Run2, settings] = lowhigh_get_data_FSL(settings,kind,run,file,mask,get_at_this_preprocessed_step);

kind = 'RestingState';
run = 1;
[RestingStateRun1, mask_RestingStateRun1, settings] = lowhigh_get_data_FSL(settings,kind,run,file,mask,get_at_this_preprocessed_step);
run = 2;
[RestingStateRun2, mask_RestingStateRun2, settings] = lowhigh_get_data_FSL(settings,kind,run,file,mask,get_at_this_preprocessed_step);
