

function [output_data, settings] = real_get_data_HCP(settings,run)

settings = real_concatenate_folders_strings_HCP(settings);
        
disp('Getting RestingState');

nTR = settings.functional.restingstate.lastTR - settings.functional.restingstate.firstTR + 1;

firstTR = settings.functional.restingstate.firstTR;
lastTR = settings.functional.restingstate.lastTR;

filename = strcat(settings.restingstate.run(run).fullpath,'.nii');

nifti_file = nifti(filename);
dat = nifti_file.dat;
BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),1:nTR) = BOLD(:,:,:,firstTR:lastTR).*dat.scl_slope + dat.scl_inter;

settings.restingstate.run(run).dim = dat.dim;
settings.restingstate.run(run).dtype = dat.dtype;
settings.restingstate.run(run).offset = dat.offset;
settings.restingstate.run(run).scl_slope = dat.scl_slope;
settings.restingstate.run(run).scl_inter = dat.scl_inter;

settings.restingstate.run(run).nifti_file = nifti_file;

end