

function [output_data, mask_data, settings] = lowhigh_get_data_FSL_mac(settings,kind,run,file,mask,get_at_this_preprocessed_step)

settings = lowhigh_concatenate_folders_strings_FSL_mac(settings,get_at_this_preprocessed_step);

all_kind = {'MOT4', 'MOT2', 'RestingState'};
idx_kind = find(strcmp(kind,all_kind));

switch idx_kind

    case 1
        
        disp('Getting MOT4');
        
        nTR = settings.functional.mot4.lastTR - settings.functional.mot4.firstTR + 1;
        
        firstTR = settings.functional.mot4.firstTR;
        lastTR = settings.functional.mot4.lastTR;
        
        filename = strcat(settings.mot4.FSL.run(run).fullpath,'/',file,'.nii');
        
        nifti_file = nifti(filename);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),1:nTR) = BOLD(:,:,:,firstTR:lastTR).*dat.scl_slope + dat.scl_inter;
        
        settings.mot4.FSL.run(run).dim = dat.dim;
        settings.mot4.FSL.run(run).dtype = dat.dtype;
        settings.mot4.FSL.run(run).offset = dat.offset;
        settings.mot4.FSL.run(run).scl_slope = dat.scl_slope;
        settings.mot4.FSL.run(run).scl_inter = dat.scl_inter;

        settings.mot4.FSL.run(run).nifti_file = nifti_file;
        
        maskname = strcat(settings.mot4.FSL.run(run).fullpath,'/',mask,'.nii');
        
        nifti_file = nifti(maskname);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        mask_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3)) = BOLD(:,:,:).*dat.scl_slope + dat.scl_inter;
        
        
    case 2
        
        disp('Getting MOT2');
        
        nTR = settings.functional.mot2.lastTR - settings.functional.mot2.firstTR + 1;
        
        firstTR = settings.functional.mot2.firstTR;
        lastTR = settings.functional.mot2.lastTR;
        
        filename = strcat(settings.mot2.FSL.run(run).fullpath,'/',file,'.nii');
        
        nifti_file = nifti(filename);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),1:nTR) = BOLD(:,:,:,firstTR:lastTR).*dat.scl_slope + dat.scl_inter;

        settings.mot2.FSL.run(run).dim = dat.dim;
        settings.mot2.FSL.run(run).dtype = dat.dtype;
        settings.mot2.FSL.run(run).offset = dat.offset;
        settings.mot2.FSL.run(run).scl_slope = dat.scl_slope;
        settings.mot2.FSL.run(run).scl_inter = dat.scl_inter;

        settings.mot2.FSL.run(run).nifti_file = nifti_file;
        
        maskname = strcat(settings.mot2.FSL.run(run).fullpath,'/',mask,'.nii');
        
        nifti_file = nifti(maskname);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        mask_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3)) = BOLD(:,:,:).*dat.scl_slope + dat.scl_inter;
        
    case 3

        disp('Getting RestingState');
        
        nTR = settings.functional.restingstate.lastTR - settings.functional.restingstate.firstTR + 1;
        
        firstTR = settings.functional.restingstate.firstTR;
        lastTR = settings.functional.restingstate.lastTR;
        
        filename = strcat(settings.restingstate.FSL.run(run).fullpath,'/',file,'.nii');
        
        nifti_file = nifti(filename);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),1:nTR) = BOLD(:,:,:,firstTR:lastTR).*dat.scl_slope + dat.scl_inter;

        settings.restingstate.FSL.run(run).dim = dat.dim;
        settings.restingstate.FSL.run(run).dtype = dat.dtype;
        settings.restingstate.FSL.run(run).offset = dat.offset;
        settings.restingstate.FSL.run(run).scl_slope = dat.scl_slope;
        settings.restingstate.FSL.run(run).scl_inter = dat.scl_inter;

        settings.restingstate.FSL.run(run).nifti_file = nifti_file;
        
        maskname = strcat(settings.restingstate.FSL.run(run).fullpath,'/',mask,'.nii');
        
        nifti_file = nifti(maskname);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        mask_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3)) = BOLD(:,:,:).*dat.scl_slope + dat.scl_inter;
          
end

end