
function settings = real_concatenate_folders_strings_FSL_SUBJ5_Again(settings,get_at_this_preprocessed_step)

    %%% FSL
    
    folder = settings.folders.main;
    folder = strcat(folder,'\',settings.folders.experiment);
    folder = strcat(folder,'\',settings.folders.subject);
    folder = strcat(folder,'\',settings.folders.preprocessed);
    FSL = settings.FSL.folders.main;
    
    tmp = strcat(folder,'\',settings.functional.restingstate.folder.main);
    settings.restingstate.FSL.run(1).fullpath = strcat(tmp,'\',settings.functional.restingstate.run(1).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.restingstate.FSL.run(2).fullpath = strcat(tmp,'\',settings.functional.restingstate.run(2).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.restingstate.FSL.run(3).fullpath = strcat(tmp,'\',settings.functional.restingstate.run(3).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    settings.restingstate.FSL.run(4).fullpath = strcat(tmp,'\',settings.functional.restingstate.run(4).folder,'\',FSL,'\',get_at_this_preprocessed_step);
    
end