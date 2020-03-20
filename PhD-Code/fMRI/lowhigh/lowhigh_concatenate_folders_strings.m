
function settings = lowhigh_concatenate_folders_strings(settings,get_at_this_preprocessed_step)


    %%% SPM
    
    folder = settings.folders.main;
    folder = strcat(folder,'\',settings.folders.experiment);
    folder = strcat(folder,'\',settings.folders.subject);
    settings.output.fullpath = strcat(folder,'\',settings.folders.output);

    folder = strcat(folder,'\',settings.folders.preprocessed);

    tmp = strcat(folder,'\',settings.anatomical.folder.main);
    settings.anatomica.run(1).fullpath = strcat(tmp,'\',settings.anatomical.run(1).folder);
    
    tmp = strcat(folder,'\',settings.functional.mot4.folder.main);
    settings.mot4.run(1).fullpath = strcat(tmp,'\',settings.functional.mot4.run(1).folder,'\',get_at_this_preprocessed_step);
    settings.mot4.run(2).fullpath = strcat(tmp,'\',settings.functional.mot4.run(2).folder,'\',get_at_this_preprocessed_step);
    
    tmp = strcat(folder,'\',settings.functional.mot2.folder.main);
    settings.mot2.run(1).fullpath = strcat(tmp,'\',settings.functional.mot2.run(1).folder,'\',get_at_this_preprocessed_step);
    settings.mot2.run(2).fullpath = strcat(tmp,'\',settings.functional.mot2.run(2).folder,'\',get_at_this_preprocessed_step);

    tmp = strcat(folder,'\',settings.functional.restingstate.folder.main);
    settings.restingstate.run(1).fullpath = strcat(tmp,'\',settings.functional.restingstate.run(1).folder,'\',get_at_this_preprocessed_step);
    settings.restingstate.run(2).fullpath = strcat(tmp,'\',settings.functional.restingstate.run(2).folder,'\',get_at_this_preprocessed_step);
    
end