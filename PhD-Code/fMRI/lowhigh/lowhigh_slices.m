
%%% SETTINGS

settings_jan_0805;

%settings_elena_2905;

%%% FILES

get_at_this_preprocessed_step = settings.FSL.folders.custom;

%file = settings.FSL.files.functional.custom.residual;
%file = settings.FSL.files.functional.custom.residual_voxel;
file = settings.FSL.files.functional.custom.filtered;
file = strcat(file,'-std');

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


%%% AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_AAL_space = find(AAL_img == 0);

%%% CANONICAL T1 TEMPLATE

MNI = nifti('avg152T1.nii');
MNI.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\spm8\canonical\',MNI.dat.fname);

MNI_img = MNI.dat(:,:,:);

%%% PLOT SLICES

%for ifile=1:length(filename)
for ifile=1
   
    filedata = nifti(filename{ifile});
    img = filedata.dat(:,:,:);
    
    %img(idx_AAL_space) = 0;
    
    %img(img>50) = 0;
    
%     figure;
%     plot_brain2d(img,4,4,1);
%     colorbar;
%     figure;
%     plot_brain2d(img,4,4,2);
%     colorbar;
%     figure;
%     plot_brain2d(img,4,4,3);
%     colorbar;

    plot_my_slices(img,MNI_img);

end