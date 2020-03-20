
function [output_data, settings] = lowhigh_get_data(settings,kind,run,get_at_this_preprocessed_step,prefix_for_the_preprocessed_step)

settings = lowhigh_concatenate_folders_strings(settings,get_at_this_preprocessed_step);

all_kind = {'MOT4', 'MOT2', 'RestingState'};
idx_kind = find(strcmp(kind,all_kind));

new_iTR = 1;

switch idx_kind

    case 1
        
        iTR = 0;
        
        disp('Getting MOT4');
        for i=settings.functional.mot4.firstTR:settings.functional.mot4.lastTR

            iTR = iTR + 1;
            
            get_iTR_number;

            prefix = settings.functional.mot4.run(run).prefix;
            file = strcat(settings.mot4.run(run).fullpath,'\',prefix_for_the_preprocessed_step,prefix,iTR_number,'-','0',iTR_number,'-','01','.hdr');
            
            nifti_file = nifti(file);
            dat = nifti_file.dat;
            BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
            output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),iTR) = BOLD(:,:,:).*dat.scl_slope + dat.scl_inter;

        end
        
        settings.mot4.run(run).dim = dat.dim;
        settings.mot4.run(run).dtype = dat.dtype;
        settings.mot4.run(run).offset = dat.offset;
        settings.mot4.run(run).scl_slope = dat.scl_slope;
        settings.mot4.run(run).scl_inter = dat.scl_inter;

        settings.mot4.run(run).nifti_file = nifti_file;

    case 2
        
        iTR = 0;
        
        disp('Getting MOT2');
        for i=settings.functional.mot2.firstTR:settings.functional.mot2.lastTR

            iTR = iTR + 1;
             
            get_iTR_number;

            prefix = settings.functional.mot2.run(run).prefix;
            file = strcat(settings.mot2.run(run).fullpath,'\',prefix_for_the_preprocessed_step,prefix,iTR_number,'-','0',iTR_number,'-','01','.hdr');

            nifti_file = nifti(file);
            dat = nifti_file.dat;
            BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
            output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),iTR) = BOLD(:,:,:).*dat.scl_slope + dat.scl_inter;

        end
        
        settings.mot2.run(run).dim = dat.dim;
        settings.mot2.run(run).dtype = dat.dtype;
        settings.mot2.run(run).offset = dat.offset;
        settings.mot2.run(run).scl_slope = dat.scl_slope;
        settings.mot2.run(run).scl_inter = dat.scl_inter;

        settings.mot2.run(run).nifti_file = nifti_file;
    
    case 3
        
        iTR = 0;

        disp('Getting RestingState');
        for i=settings.functional.restingstate.firstTR:settings.functional.restingstate.lastTR

            iTR = iTR + 1;
             
            get_iTR_number;

            prefix = settings.functional.restingstate.run(run).prefix;
            file = strcat(settings.restingstate.run(run).fullpath,'\',prefix_for_the_preprocessed_step,prefix,iTR_number,'-','0',iTR_number,'-','01','.hdr');

            nifti_file = nifti(file);
            dat = nifti_file.dat;
            BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
            output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),iTR) = BOLD(:,:,:).*dat.scl_slope + dat.scl_inter;

        end
        
        settings.restingstate.run(run).dim = dat.dim;
        settings.restingstate.run(run).dtype = dat.dtype;
        settings.restingstate.run(run).offset = dat.offset;
        settings.restingstate.run(run).scl_slope = dat.scl_slope;
        settings.restingstate.run(run).scl_inter = dat.scl_inter;

        settings.restingstate.run(run).nifti_file = nifti_file;
    
end

end



