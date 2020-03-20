
function settings = real_concatenate_folders_strings_FSL_Martin(settings,get_at_this_preprocessed_step)

    %%% FSL
    
    folder = settings.folders.main;
    folder = strcat(folder,'\',settings.folders.experiment);
    % folder = strcat(folder,'\',settings.folders.subject);
    folder = strcat(folder,'\',settings.folders.preprocessed);
    FSL = settings.FSL.folders.main;
    
    tmp = strcat(folder,'\',settings.functional.restingstate.folder.main);
    settings.restingstate.FSL.run(1).fullpath = strcat(tmp,'\','Run-1','\',FSL,'\',get_at_this_preprocessed_step);
    settings.restingstate.FSL.run(2).fullpath = strcat(tmp,'\','Run-2','\',FSL,'\',get_at_this_preprocessed_step);
    settings.restingstate.FSL.run(3).fullpath = strcat(tmp,'\','Run-3','\',FSL,'\',get_at_this_preprocessed_step);
    settings.restingstate.FSL.run(4).fullpath = strcat(tmp,'\','Run-4','\',FSL,'\',get_at_this_preprocessed_step);
    
 
end