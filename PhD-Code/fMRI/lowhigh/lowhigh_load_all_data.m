
kind = 'MOT4';
run = 1;
[MOT4Run1, settings] = lowhigh_get_data(settings,kind,run,get_at_this_preprocessed_step,prefix_for_the_preprocessed_step);
run = 2;
[MOT4Run2, settings] = lowhigh_get_data(settings,kind,run,get_at_this_preprocessed_step,prefix_for_the_preprocessed_step);

kind = 'MOT2';
run = 1;
[MOT2Run1, settings] = lowhigh_get_data(settings,kind,run,get_at_this_preprocessed_step,prefix_for_the_preprocessed_step);
run = 2;
[MOT2Run2, settings] = lowhigh_get_data(settings,kind,run,get_at_this_preprocessed_step,prefix_for_the_preprocessed_step);

kind = 'RestingState';
run = 1;
[RestingStateRun1, settings] = lowhigh_get_data(settings,kind,run,get_at_this_preprocessed_step,prefix_for_the_preprocessed_step);
run = 2;
[RestingStateRun2, settings] = lowhigh_get_data(settings,kind,run,get_at_this_preprocessed_step,prefix_for_the_preprocessed_step);
