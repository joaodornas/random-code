

function [output_data, settings] = real_get_data_FSL_Petra(settings,run,file,get_at_this_preprocessed_step)

settings = real_concatenate_folders_strings_FSL_Petra(settings,get_at_this_preprocessed_step);
        
disp('Getting RestingState');

nTR = settings.functional.restingstate.lastTR - settings.functional.restingstate.firstTR + 1;

firstTR = settings.functional.restingstate.firstTR;
lastTR = settings.functional.restingstate.lastTR;

filename = strcat(settings.restingstate.FSL.run(run).fullpath,'\',file,'.nii');

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

end