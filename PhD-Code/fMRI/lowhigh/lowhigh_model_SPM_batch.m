function lowhigh_model_SPM_batch(matlabbatch)

%settings_jan_0805;
settings_elena_2905;

model = 'MOT4 > RestingState';
updateModel2(settings,model,matlabbatch);

model = 'MOT2 > RestingState';
updateModel2(settings,model,matlabbatch);

model = 'MOT4 > MOT2';
updateModel2(settings,model,matlabbatch);

% model = 'MOT4 > MOT2 > RestingState';
% updateModel3(settings,model,matlabbatch);

end

function updateModel2(settings,model,matlabbatch)

%get_at_this_preprocessed_step = settings.folders.smooth.name;
%prefix_for_the_preprocessed_step = settings.folders.smooth.prefix;

get_at_this_preprocessed_step = settings.FSL.folders.custom;

file = settings.FSL.files.functional.custom.residual_voxel;
volume_prefix = 'residual';
volume_folder = 'volumes';

settings = lowhigh_concatenate_folders_strings_FSL(settings,get_at_this_preprocessed_step);

iTR = 0;
iTR_mot4 = 0;
iTR_RestingState = 0;

for iRun=1:10
    
    if iRun <= 5, run = 1; else run = 2; end
    
    if strcmpi(model,'MOT4 > RestingState')
        
        condition1(run).fullpath = settings.mot4.FSL.run(run).fullpath;    
        condition2(run).fullpath = settings.restingstate.FSL.run(run).fullpath;
    
        name = 'MOT4-RestingState';
        
    elseif strcmpi(model,'MOT2 > RestingState')
        
        condition1(run).fullpath = settings.mot2.FSL.run(run).fullpath;
        condition2(run).fullpath = settings.restingstate.FSL.run(run).fullpath;
        
        name = 'MOT2-RestingState';
        
    elseif strcmpi(model,'MOT4 > MOT2')
        
        condition1(run).fullpath = settings.mot4.FSL.run(run).fullpath;
        condition2(run).fullpath = settings.mot2.FSL.run(run).fullpath;
        
        name = 'MOT4-MOT2';
        
    end
    
    if iTR_mot4 == 300, iTR_mot4 = 0; end
    if iTR_RestingState == 300, iTR_RestingState = 0; end
    
    for iV=1:120
       
        new_iTR = 1;
        
        if (iV/12) <= 1 || (iV/12) > 2 && (iV/12) <= 3 || (iV/12) > 4 && (iV/12) <= 5 || (iV/12) > 6 && (iV/12) <= 7 || (iV/12) > 8 && (iV/12) <= 9
        
            iTR_mot4 = iTR_mot4 + 1;
            
            iTR_FSL = iTR_mot4 - 1;
            
            get_iTR_number;

            file = strcat(condition1(run).fullpath,'\',volume_folder,'\',volume_prefix,iTR_number_FSL,'.nii');
 
        end
        
        if (iV/12) > 1 && (iV/12) <= 2 || (iV/12) > 3 && (iV/12) <= 4 || (iV/12) > 5 && (iV/12) <= 6 || (iV/12) > 7 && (iV/12) <= 8 || (iV/12) > 9
        
            iTR_RestingState = iTR_RestingState + 1;
            
            iTR_FSL = iTR_RestingState - 1;
            
            get_iTR_number;

            file = strcat(condition2(run).fullpath,'\',volume_folder,'\',volume_prefix,iTR_number_FSL,'.nii');
        end
        
        matlabbatch{1,1}.spm.stats.fmri_spec.sess(iRun).scans{iV,1} = strcat(file,',1');
        
    end
    
end

save(strcat('SPM-',name,'.mat'),'matlabbatch');

end

function updateModel3(settings,model,matlabbatch)

%get_at_this_preprocessed_step = settings.folders.smooth.name;
%prefix_for_the_preprocessed_step = settings.folders.smooth.prefix;

get_at_this_preprocessed_step = settings.FSL.folders.custom;

file = settings.FSL.files.functional.custom.residual_voxel;
volume_prefix = 'residual';
volume_folder = 'volumes';

settings = lowhigh_concatenate_folders_strings_FSL(settings,get_at_this_preprocessed_step);

iTR = 0;
iTR_mot4 = 0 + 30;
iTR_mot2 = 0 + 30;
iTR_RestingState = 0;

for iRun=1:10
    
    if iRun <= 5, run = 1; else run = 2; end
    
    if strcmpi(model,'MOT4 > MOT2 > RestingState')
        
        condition1(run).fullpath = settings.mot4.FSL.run(run).fullpath;
        condition2(run).fullpath = settings.mot2.FSL.run(run).fullpath;
        condition3(run).fullpath = settings.restingstate.FSL.run(run).fullpath;
        
        name = 'MOT4-MOT2-RestingState';
        
    end
    
    if iTR_mot4 == 300, iTR_mot4 = 0; end
    if iTR_mot2 == 300, iTR_mot2 = 0; end
    if iTR_RestingState == 300, iTR_RestingState = 0; end
    
    for iV=1:180
       
        new_iTR = 1;
        
        if (iV/12) <= 1 || (iV/12) > 2 && (iV/12) <= 3 || (iV/12) > 4 && (iV/12) <= 5 || (iV/12) > 6 && (iV/12) <= 7 || (iV/12) > 8 && (iV/12) <= 9
        
            iTR_mot4 = iTR_mot4 + 1;
            
            iTR_FSL = iTR_mot4 - 1;
            
            get_iTR_number;
            
            prefix = condition1(run).prefix;
            file = strcat(condition1(run).fullpath,'\',volume_folder,'\',volume_prefix,iTR_number_FSL,'.nii');
 
        end
        
        if (iV/12) > 1 && (iV/12) <= 2 || (iV/12) > 3 && (iV/12) <= 4 || (iV/12) > 5 && (iV/12) <= 6 || (iV/12) > 7 && (iV/12) <= 8 || (iV/12) > 9 && (iV/12) <= 10
        
            iTR_mot2 = iTR_mot2 + 1;
            
            iTR_FSL = iTR_mot2 - 1;
            
            get_iTR_number;
            
            prefix = condition2(run).prefix;
            file = strcat(condition2(run).fullpath,'\',volume_folder,'\',volume_prefix,iTR_number_FSL,'.nii');
        end
        
        if (iV/12) > 10 && (iV/12) <= 11 || (iV/12) > 12 && (iV/12) <= 13 || (iV/12) > 14 && (iV/12) <= 15 || (iV/12) > 16 && (iV/12) <= 17 || (iV/12) > 18 
        
            iTR_RestingState = iTR_RestingState + 1;
            
            iTR_FSL = iTR_RestingState - 1;
            
            get_iTR_number;
            
            prefix = condition3(run).prefix;
            file = strcat(condition3(run).fullpath,'\',volume_folder,'\',volume_prefix,iTR_number_FSL,'.nii');
      
        end
        
        matlabbatch{1,1}.spm.stats.fmri_spec.sess(iRun).scans{iV,1} = strcat(file,',1');
        
    end
    
end

save(strcat('SPM-',name,'.mat'),'matlabbatch');

end






