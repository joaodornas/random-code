
function settings = real_concatenate_folders_strings_FSL(settings,get_at_this_preprocessed_step)

    %%% FSL
    
    folder = settings.folders.main;
    folder = strcat(folder,'\',settings.folders.experiment);
    folder = strcat(folder,'\',settings.folders.subject);
    folder = strcat(folder,'\',settings.folders.preprocessed);
    FSL = settings.FSL.folders.main;
    
    tmp = strcat(folder,'\',settings.functional.track.folder.main);
    settings.track.FSL.run(1).fullpath = strcat(tmp,'\',settings.functional.track.run(1).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.track.FSL.run(2).fullpath = strcat(tmp,'\',settings.functional.track.run(2).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.track.FSL.run(3).fullpath = strcat(tmp,'\',settings.functional.track.run(3).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.track.FSL.run(4).fullpath = strcat(tmp,'\',settings.functional.track.run(4).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    
    tmp = strcat(folder,'\',settings.functional.passive.folder.main);
    settings.passive.FSL.run(1).fullpath = strcat(tmp,'\',settings.functional.passive.run(1).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.passive.FSL.run(2).fullpath = strcat(tmp,'\',settings.functional.passive.run(2).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.passive.FSL.run(3).fullpath = strcat(tmp,'\',settings.functional.passive.run(3).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.passive.FSL.run(4).fullpath = strcat(tmp,'\',settings.functional.passive.run(4).folder,'\',FSL,'\',get_at_this_preprocessed_step);

    tmp = strcat(folder,'\',settings.functional.trials.folder.main);
    settings.trials.FSL.run(1).fullpath = strcat(tmp,'\',settings.functional.trials.run(1).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.trials.FSL.run(2).fullpath = strcat(tmp,'\',settings.functional.trials.run(2).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.trials.FSL.run(3).fullpath = strcat(tmp,'\',settings.functional.trials.run(3).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.trials.FSL.run(4).fullpath = strcat(tmp,'\',settings.functional.trials.run(4).folder,'\',FSL,'\',get_at_this_preprocessed_step);

    tmp = strcat(folder,'\',settings.functional.restingstate.folder.main);
    settings.restingstate.FSL.run(1).fullpath = strcat(tmp,'\',settings.functional.restingstate.run(1).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.restingstate.FSL.run(2).fullpath = strcat(tmp,'\',settings.functional.restingstate.run(2).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.restingstate.FSL.run(3).fullpath = strcat(tmp,'\',settings.functional.restingstate.run(3).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.restingstate.FSL.run(4).fullpath = strcat(tmp,'\',settings.functional.restingstate.run(4).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    
end