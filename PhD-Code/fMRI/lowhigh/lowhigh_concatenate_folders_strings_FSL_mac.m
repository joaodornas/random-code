
function settings = lowhigh_concatenate_folders_strings_FSL_mac(settings,get_at_this_preprocessed_step)

    %%% FSL
    
    folder = settings.folders.main;
    folder = strcat(folder,'/',settings.folders.experiment);
    folder = strcat(folder,'/',settings.folders.subject);
    folder = strcat(folder,'/',settings.folders.preprocessed);
    FSL = settings.FSL.folders.main;
    
    tmp = strcat(folder,'/',settings.functional.mot4.folder.main);
    settings.mot4.FSL.run(1).fullpath = strcat(tmp,'/',settings.functional.mot4.run(1).folder,'/',FSL,'/',get_at_this_preprocessed_step);
    settings.mot4.FSL.run(2).fullpath = strcat(tmp,'/',settings.functional.mot4.run(2).folder,'/',FSL,'/',get_at_this_preprocessed_step);
    
    tmp = strcat(folder,'/',settings.functional.mot2.folder.main);
    settings.mot2.FSL.run(1).fullpath = strcat(tmp,'/',settings.functional.mot2.run(1).folder,'/',FSL,'/',get_at_this_preprocessed_step);
    settings.mot2.FSL.run(2).fullpath = strcat(tmp,'/',settings.functional.mot2.run(2).folder,'/',FSL,'/',get_at_this_preprocessed_step);

    tmp = strcat(folder,'/',settings.functional.restingstate.folder.main);
    settings.restingstate.FSL.run(1).fullpath = strcat(tmp,'/',settings.functional.restingstate.run(1).folder,'/',FSL,'/',get_at_this_preprocessed_step);
    settings.restingstate.FSL.run(2).fullpath = strcat(tmp,'/',settings.functional.restingstate.run(2).folder,'/',FSL,'/',get_at_this_preprocessed_step);
    
end