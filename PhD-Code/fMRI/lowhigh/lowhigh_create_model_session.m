
settings_jan_0502;

get_at_this_preprocessed_step = settings.folders.smooth.name;
prefix_for_the_preprocessed_step = settings.folders.smooth.prefix;

settings = lowhigh_concatenate_folders_strings(settings,get_at_this_preprocessed_step);

for run=1:2
    
    new_iTR = 0;
    
    iTR = 0;
    
    for i=settings.functional.mot4.firstTR:settings.functional.mot4.lastTR

        iTR = iTR + 1;
        
        new_iTR = new_iTR + 1;

        get_iTR_number;

        prefix = settings.functional.mot4.run(run).prefix;

        originalfilehdr = strcat(settings.mot4.run(run).fullpath,'\',prefix_for_the_preprocessed_step,prefix,iTR_number,'-','0',iTR_number,'-','01','.hdr');
        originalfileimg = strcat(settings.mot4.run(run).fullpath,'\',prefix_for_the_preprocessed_step,prefix,iTR_number,'-','0',iTR_number,'-','01','.img');

        destinationfilehdr = strcat('MOT4-RestingState-Run-',int2str(run),'-',new_iTR_number,'.hdr');
        destinationfileimg = strcat('MOT4-RestingState-Run-',int2str(run),'-',new_iTR_number,'.img');
        
        copyfile(originalfilehdr,destinationfilehdr,'f');
        copyfile(originalfileimg,destinationfileimg,'f');

    end
    
    iTR = 0;
    
    for i=settings.functional.restingstate.firstTR:settings.functional.restingstate.lastTR

        iTR = iTR + 1;
        
        new_iTR = new_iTR + 1;

        get_iTR_number;

        prefix = settings.functional.restingstate.run(run).prefix;
        
        originalfilehdr = strcat(settings.restingstate.run(run).fullpath,'\',prefix_for_the_preprocessed_step,prefix,iTR_number,'-','0',iTR_number,'-','01','.hdr');
        originalfileimg = strcat(settings.restingstate.run(run).fullpath,'\',prefix_for_the_preprocessed_step,prefix,iTR_number,'-','0',iTR_number,'-','01','.img');

        destinationfilehdr = strcat('MOT4-RestingState-Run-',int2str(run),'-',new_iTR_number,'.hdr');
        destinationfileimg = strcat('MOT4-RestingState-Run-',int2str(run),'-',new_iTR_number,'.img');
        
        copyfile(originalfilehdr,destinationfilehdr,'f');
        copyfile(originalfileimg,destinationfileimg,'f');
        
    end

end