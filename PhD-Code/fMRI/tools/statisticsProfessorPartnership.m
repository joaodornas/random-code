function statisticsProfessorPartnership

getCopyOfData;

end


function getCopyOfData

nRuns = 4;
nSubjects = 8;
nTR = 150;

all_settings = getAllSettings;

iiRun_track = 0;
iiRun_passive = 0;
iiRun_rest = 0;

for iSubject=1:nSubjects

    settings = all_settings(iSubject).settings;

    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
    
    settings = real_concatenate_folders_strings_FSL(settings,get_at_this_preprocessed_step);

    all_kind = {'Track', 'Passive', 'RestingState'};
    
    for run=1:nRuns
        
        iiRun_track = iiRun_track + 1;

        nTR = settings.functional.track.lastTR - settings.functional.track.firstTR + 1;
        
        firstTR = settings.functional.track.firstTR;
        lastTR = settings.functional.track.lastTR;
        
        filename = strcat(settings.track.FSL.run(run).fullpath,'\',file,'.nii');
        
        nifti_file = nifti(filename);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),1:nTR) = BOLD(:,:,:,firstTR:lastTR).*dat.scl_slope + dat.scl_inter;
        
        dim = [91 109 91 150];
        dtype = dat.dtype;
        offset = dat.offset;
        scl_slope = dat.scl_slope;
        scl_inter = dat.scl_inter;
        
        input_data = output_data;
        
        descrip = 'Track';
        
        fname = strcat('Track','-',int2str(iiRun_track),'.nii');
        
        real_save_image;
        
    end  
        
    for run=1:nRuns
        
        iiRun_passive = iiRun_passive + 1;

        nTR = settings.functional.passive.lastTR - settings.functional.passive.firstTR + 1;
        
        firstTR = settings.functional.passive.firstTR;
        lastTR = settings.functional.passive.lastTR;
        
        filename = strcat(settings.passive.FSL.run(run).fullpath,'\',file,'.nii');
        
        nifti_file = nifti(filename);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),1:nTR) = BOLD(:,:,:,firstTR:lastTR).*dat.scl_slope + dat.scl_inter;
        
        dim = [91 109 91 150];
        dtype = dat.dtype;
        offset = dat.offset;
        scl_slope = dat.scl_slope;
        scl_inter = dat.scl_inter;
        
        input_data = output_data;
        
        descrip = 'Passive';
        
        fname = strcat('Passive','-',int2str(iiRun_passive),'.nii');
        
        real_save_image;
        
    end
    
    for run=1:nRuns
        
        iiRun_rest = iiRun_rest + 1;

        nTR = settings.functional.restingstate.lastTR - settings.functional.restingstate.firstTR + 1;
        
        firstTR = settings.functional.restingstate.firstTR;
        lastTR = settings.functional.restingstate.lastTR;
        
        filename = strcat(settings.restingstate.FSL.run(run).fullpath,'\',file,'.nii');
        
        nifti_file = nifti(filename);
        dat = nifti_file.dat;
        BOLD = file_array(dat.fname,dat.dim,dat.dtype,dat.offset,dat.scl_slope,dat.scl_inter);
        output_data(1:dat.dim(1),1:dat.dim(2),1:dat.dim(3),1:nTR) = BOLD(:,:,:,firstTR:lastTR).*dat.scl_slope + dat.scl_inter;
        
        dim = [91 109 91 150];
        dtype = dat.dtype;
        offset = dat.offset;
        scl_slope = dat.scl_slope;
        scl_inter = dat.scl_inter;
        
        input_data = output_data;
        
        descrip = 'Rest';
        
        fname = strcat('Rest','-',int2str(iiRun_rest),'.nii');
        
        real_save_image;
        
    end
    
end
    
end

