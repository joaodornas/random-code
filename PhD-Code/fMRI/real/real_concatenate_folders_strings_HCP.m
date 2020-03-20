
function settings = real_concatenate_folders_strings_HCP(settings)

    %%% FSL
    
    folder = settings.folders.main;
    folder = strcat(folder,'\',settings.folders.experiment);
    folder = strcat(folder,'\',settings.folders.preprocessed);
    
    settings.restingstate.run(1).fullpath = strcat(folder,'\',settings.folders.subject(1).path,'\',settings.folders.subject(1).results(1).path,'\',settings.functional.restingstate.subject(1).run(1).file);
    settings.restingstate.run(2).fullpath = strcat(folder,'\',settings.folders.subject(1).path,'\',settings.folders.subject(1).results(2).path,'\',settings.functional.restingstate.subject(1).run(2).file)
    settings.restingstate.run(3).fullpath = strcat(folder,'\',settings.folders.subject(2).path,'\',settings.folders.subject(2).results(1).path,'\',settings.functional.restingstate.subject(2).run(1).file)
    settings.restingstate.run(4).fullpath = strcat(folder,'\',settings.folders.subject(2).path,'\',settings.folders.subject(2).results(2).path,'\',settings.functional.restingstate.subject(2).run(2).file)
    
 
end