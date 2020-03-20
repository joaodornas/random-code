

function [output_data, mask_data, settings] = real_get_data_FSL_MAC(settings,kind,run,file,mask,get_at_this_preprocessed_step,mm)

settings = real_concatenate_folders_strings_FSL_MAC(settings,get_at_this_preprocessed_step);

all_kind = {'Track', 'Passive', 'Trials', 'RestingState'};
idx_kind = find(strcmp(kind,all_kind));

switch idx_kind

    case 1
        
        disp('Getting Track');
        
        nTR = settings.functional.track.lastTR - settings.functional.track.firstTR + 1;
        
        firstTR = settings.functional.track.firstTR;
        lastTR = settings.functional.track.lastTR;
        
        filename = strcat(settings.track.FSL.run(run).fullpath,'/',file,'-',int2str(mm),'mm','.nii');
        
        nifti_file = nifti(filename);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),1:nTR) = BOLD(:,:,:,firstTR:lastTR).*dat.scl_slope + dat.scl_inter;
        
        settings.track.FSL.run(run).dim = dat.dim;
        settings.track.FSL.run(run).dtype = dat.dtype;
        settings.track.FSL.run(run).offset = dat.offset;
        settings.track.FSL.run(run).scl_slope = dat.scl_slope;
        settings.track.FSL.run(run).scl_inter = dat.scl_inter;

        settings.track.FSL.run(run).nifti_file = nifti_file;
        
        maskname = strcat(settings.track.FSL.run(run).fullpath,'/',mask,'.nii');
        
        nifti_file = nifti(maskname);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        mask_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3)) = BOLD(:,:,:).*dat.scl_slope + dat.scl_inter;
        
    case 2
        
        disp('Getting Passive');
        
        nTR = settings.functional.passive.lastTR - settings.functional.passive.firstTR + 1;
        
        firstTR = settings.functional.passive.firstTR;
        lastTR = settings.functional.passive.lastTR;
        
        filename = strcat(settings.passive.FSL.run(run).fullpath,'/',file,'.nii');
        
        nifti_file = nifti(filename);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),1:nTR) = BOLD(:,:,:,firstTR:lastTR).*dat.scl_slope + dat.scl_inter;
        
        settings.passive.FSL.run(run).dim = dat.dim;
        settings.passive.FSL.run(run).dtype = dat.dtype;
        settings.passive.FSL.run(run).offset = dat.offset;
        settings.passive.FSL.run(run).scl_slope = dat.scl_slope;
        settings.passive.FSL.run(run).scl_inter = dat.scl_inter;

        settings.passive.FSL.run(run).nifti_file = nifti_file;
        
        maskname = strcat(settings.passive.FSL.run(run).fullpath,'/',mask,'.nii');
        
        nifti_file = nifti(maskname);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        mask_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3)) = BOLD(:,:,:).*dat.scl_slope + dat.scl_inter;
        
    case 3
        
        disp('Getting Trials');
        
        nTR = settings.functional.trials.lastTR - settings.functional.trials.firstTR + 1;
        
        firstTR = settings.functional.trials.firstTR;
        lastTR = settings.functional.trials.lastTR;
        
        filename = strcat(settings.trials.FSL.run(run).fullpath,'/',file,'.nii');
        
        nifti_file = nifti(filename);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),1:nTR) = BOLD(:,:,:,firstTR:lastTR).*dat.scl_slope + dat.scl_inter;
        
        settings.trials.FSL.run(run).dim = dat.dim;
        settings.trials.FSL.run(run).dtype = dat.dtype;
        settings.trials.FSL.run(run).offset = dat.offset;
        settings.trials.FSL.run(run).scl_slope = dat.scl_slope;
        settings.trials.FSL.run(run).scl_inter = dat.scl_inter;

        settings.trials.FSL.run(run).nifti_file = nifti_file;
        
        maskname = strcat(settings.trials.FSL.run(run).fullpath,'/',mask,'.nii');
        
        nifti_file = nifti(maskname);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        mask_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3)) = BOLD(:,:,:).*dat.scl_slope + dat.scl_inter;

    case 4
        
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