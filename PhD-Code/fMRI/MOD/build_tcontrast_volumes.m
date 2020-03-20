
idx_k_MOD_Color_mv1 = [146, 1100, 2135, 2385, 2896, 3220, 3764, 4511, 5426, 5974];

idx_k_MOD_Color_mv2 = [106, 1170, 2279, 2472, 3992, 4988, 5386, 5608, 5940, 6589];

idx_frame_MOD_Color_mv1 = idx_k_MOD_Color_mv1.*3;
idx_frame_MOD_Color_mv2 = idx_k_MOD_Color_mv2.*3;

FrameRate = 60;
TR = 2;

idx_TR_MOD_Color_mv1 = round( idx_frame_MOD_Color_mv1 ./ (FrameRate*TR) ) + 1;
idx_TR_MOD_Color_mv2 = round( idx_frame_MOD_Color_mv2 ./ (FrameRate*TR) ) + 1;

idx_TR_MOD_Color_mv1(idx_TR_MOD_Color_mv1<15) = [];
idx_TR_MOD_Color_mv2(idx_TR_MOD_Color_mv2<15) = [];

idx_TR_MOD_Color_mv1(idx_TR_MOD_Color_mv1>150) = [];
idx_TR_MOD_Color_mv2(idx_TR_MOD_Color_mv2>150) = [];

idx_TR_MOD_Color_mv1 = idx_TR_MOD_Color_mv1 - 15;
idx_TR_MOD_Color_mv2 = idx_TR_MOD_Color_mv2 - 15;

TR_interval = 15;
nRuns = 32;

all_settings = getAllSettings;

iirun_track = 0;
iirun_passive = 0;
for iSet=1:length(all_settings)
    
    settings = all_settings(iSet).settings;
    
    %% LOAD DATA

    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;

    kind = 'Track';
    for irun=1:4
        iirun_track = iirun_track + 1;
        [Track(iirun_track).run, out_mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end

end

%%% LOAD AAL

load_example = nifti('filtered_func_data.nii');
load_example.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-1-22-10-2015\preprocessed\T2-Stimulus-MOT-4-Balls-Track\Run-1-1-1\FSL\Melodic-Fieldmap.ica\',load_example.dat.fname);

for iMOD=2:length(idx_TR_MOD_Color_mv1)
   
    MOD_volume = zeros(91,109,91,(nRuns/2)*2*TR_interval);
    
    iiRun = 0;
    for iRun=1:2:nRuns
        
       iiRun = iiRun + 1;
        
       MOD_volume(:,:,:,((iiRun-1)*(2*TR_interval)+1):((iiRun)*(2*TR_interval))) = Track(iRun).run(:,:,:,(idx_TR_MOD_Color_mv1(iMOD)-TR_interval+1):(idx_TR_MOD_Color_mv1(iMOD)+TR_interval));
        
    end
    
    nifti_file = load_example;
    offset = load_example.dat.offset;
    scl_slope = load_example.dat.scl_slope;
    scl_inter = load_example.dat.scl_inter;

    dtype = 'FLOAT32';
    offset = 0;

    dim = [91,109,91,(nRuns/2)*2*TR_interval];

    descrip = 'MOD';

    fname = strcat('MOD-Color-mv1-k-',int2str(idx_k_MOD_Color_mv1(iMOD)),'.nii');
    input_data = MOD_volume; 
    real_save_image;

end

for iMOD=2:length(idx_TR_MOD_Color_mv2)
   
    MOD_volume = zeros(91,109,91,(nRuns/2)*2*TR_interval);
    
    iiRun = 0;
    for iRun=2:2:nRuns
        
       iiRun = iiRun + 1;
        
       MOD_volume(:,:,:,((iiRun-1)*(2*TR_interval)+1):((iiRun)*(2*TR_interval))) = Track(iRun).run(:,:,:,(idx_TR_MOD_Color_mv1(iMOD)-TR_interval+1):(idx_TR_MOD_Color_mv1(iMOD)+TR_interval));
        
    end
    
    nifti_file = load_example;
    offset = load_example.dat.offset;
    scl_slope = load_example.dat.scl_slope;
    scl_inter = load_example.dat.scl_inter;

    dtype = 'FLOAT32';
    offset = 0;

    dim = [91,109,91,(nRuns/2)*2*TR_interval];

    descrip = 'MOD';

    fname = strcat('MOD-Color-mv2-k-',int2str(idx_k_MOD_Color_mv2(iMOD)),'.nii');
    input_data = MOD_volume; 
    real_save_image;
    
end

