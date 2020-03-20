function lowhigh_rendering_std

%settings_jan_0502;

%settings_jan_0805;

settings_elena_2905;

doTheRendering(settings);

return

function doTheRendering(settings)

rendfile = 'K:\Dropbox (Uni Magdeburg)\_TOOLBOX\spm8\canonical\cortex_20484.surf.gii';

brt = 1;

%%% NAMES OF IMAGES - ALL RUNS

%get_at_this_preprocessed_step = settings.FSL.folders.warped;
get_at_this_preprocessed_step = settings.FSL.folders.custom;
%file = settings.FSL.files.functional.warped;
%file = settings.FSL.files.functional.despike;
%file = settings.FSL.files.functional.custom.residual;
%file = settings.FSL.files.functional.custom.wavelet_voxel;
file = settings.FSL.files.functional.custom.residual_voxel;

%analysis = 'std';
analysis = 'mean';
%kindoffile = 'warped';
%kindoffile = 'despike-25';
%kindoffile = 'residual';
%kindoffile = 'wavelet-voxel';
kindoffile = 'residual-voxel';

file = strcat(file,'-mean');
%file = strcat(file,'-std');

settings = lowhigh_concatenate_folders_strings_FSL(settings,get_at_this_preprocessed_step);

run = 1;
filename{1} = strcat(settings.mot4.FSL.run(run).fullpath,'\',file,'.nii');
filename{2} = strcat(settings.mot2.FSL.run(run).fullpath,'\',file,'.nii');
filename{3} = strcat(settings.restingstate.FSL.run(run).fullpath,'\',file,'.nii');

run = 2;
filename{4} = strcat(settings.mot4.FSL.run(run).fullpath,'\',file,'.nii');
filename{5} = strcat(settings.mot2.FSL.run(run).fullpath,'\',file,'.nii');
filename{6} = strcat(settings.restingstate.FSL.run(run).fullpath,'\',file,'.nii');

names{1} = 'High-Attention-Run1';
names{2} = 'Low-Attention-Run1';
names{3} = 'Resting-State-Run1';
names{4} = 'High-Attention-Run2';
names{5} = 'Low-Attention-Run2';
names{6} = 'Resting-State-Run2';

%%% SPM DATA FILES

for ifile=1:length(filename)

    img(ifile).dat = fromImageToSPMdat(filename(ifile));

end

%%% MIN AND MAX VALUES

for ifile=1:length(filename)
    
    [img(ifile).min_value, img(ifile).max_value] = getMinMaxFromData(img(ifile).dat);

end

%%% MIN AND MAX VALUES OVERALL

min_value = floor(min([img(1).min_value, img(2).min_value, img(3).min_value, img(4).min_value, img(5).min_value, img(6).min_value]));
max_value = ceil(max([img(1).max_value, img(2).max_value, img(3).max_value, img(4).max_value, img(5).max_value, img(6).max_value]));

%min_value = 1.96; %% zscore p-value 0.05
%max_value = 2.16;

%min_value = 0;
%max_value = 600/2;
%max_value = 120;
parcellation = false;

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;
idx_brain = find(AAL_img);
idx_space = find(AAL_img == 0);

for ifile=1:length(filename)

    if parcellation
    
        img(ifile).dat.t(idx_space) = 0;
        
    else
        
        img(ifile).dat.t(find(img(ifile).dat.t<min_value)) = 0;
        img(ifile).dat.t(find(img(ifile).dat.t>max_value)) = 0;

    end
    
end

%%% SET GLOBAL MIN AND MAX

setGlobalMinMax(min_value,max_value);

%%% PLOT 3D RENDERING

disp('Rendering Image 1');

for ifile=1:length(filename)

    snapshot_label = strcat(settings.folders.experiment,'-',settings.folders.subject,'-',names{ifile},'-',analysis,'-',kindoffile);
    setRenderingSnapshot(true,snapshot_label);

    spm_render_my(img(ifile).dat,brt,rendfile,filename{ifile});

end

return

function [min_value, max_value] = getMinMaxFromData(dat)

min_value = min(dat.t);
max_value = max(dat.t);

return


