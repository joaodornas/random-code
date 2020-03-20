
function gen_global_lFCD_dataset

gen_runs;
% gen_runs_numbers;

% gen_masks;

% gen_nVoxels;

% gen_global_Voxels;

end

function gen_runs

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

%%% ALL SETTINGS DATA

all_settings = getAllSettings;
iirun_track = 0;
iirun_passive = 0;
iirun_rest = 0;
iirun = 0;
for iSet=1:length(all_settings)
    
    settings = all_settings(iSet).settings;
    
    %%% LOAD DATA

    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;

    real_load_all_data_FSL;
    
    for irun=1:4
       
        iirun = iirun + 1;
        iirun_track = iirun_track + 1;
        nifti_file = load_aal;
        offset = load_aal.dat.offset;
        scl_slope = load_aal.dat.scl_slope;
        scl_inter = load_aal.dat.scl_inter;
        dtype = 'FLOAT32';
        offset = 0;
        dim = [91 109 91 150];
        descrip = 'run';
        fname = strcat('run-track','-',int2str(iirun_track),'.nii');
        input_data = Track(irun).run; 
        real_save_image;
        %extract_time_series(Track(irun).run,AAL_img,iirun_track,'track');
        
        iirun_passive = iirun_passive + 1;
        nifti_file = load_aal;
        offset = load_aal.dat.offset;
        scl_slope = load_aal.dat.scl_slope;
        scl_inter = load_aal.dat.scl_inter;
        dtype = 'FLOAT32';
        offset = 0;
        dim = [91 109 91 150];
        descrip = 'run';
        fname = strcat('run-passive','-',int2str(iirun_passive),'.nii');
        input_data = Passive(irun).run; 
        real_save_image;
        %extract_time_series(Passive(irun).run,AAL_img,iirun_passive,'passive');
        
        iirun_rest = iirun_rest + 1;
        nifti_file = load_aal;
        offset = load_aal.dat.offset;
        scl_slope = load_aal.dat.scl_slope;
        scl_inter = load_aal.dat.scl_inter;
        dtype = 'FLOAT32';
        offset = 0;
        dim = [91 109 91 150];
        descrip = 'run';
        fname = strcat('run-rest','-',int2str(iirun_rest),'.nii');
        input_data = RestingState(irun).run; 
        real_save_image;
        %extract_time_series(RestingState(irun).run,AAL_img,iirun_rest,'rest');

    end

end

end

function gen_runs_numbers

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

%%% ALL SETTINGS DATA

all_settings = getAllSettings;
iirun_track = 0;
iirun_passive = 0;
iirun_rest = 0;
iirun_jump = 32;
for iSet=1:length(all_settings)
    
    settings = all_settings(iSet).settings;
    
    %%% LOAD DATA

    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;

    real_load_all_data_FSL;
    
    for irun=1:4
       
        iirun_track = iirun_track + 1;
        nifti_file = load_aal;
        offset = load_aal.dat.offset;
        scl_slope = load_aal.dat.scl_slope;
        scl_inter = load_aal.dat.scl_inter;
        dtype = 'FLOAT32';
        offset = 0;
        dim = [91 109 91 150];
        descrip = 'run';
        fname = strcat('run','-',int2str(iirun_track),'.nii');
        input_data = Track(irun).run; 
        real_save_image;
        %extract_time_series(Track(irun).run,AAL_img,iirun_track,'track');
        
        iirun_passive = iirun_passive + 1;
        nifti_file = load_aal;
        offset = load_aal.dat.offset;
        scl_slope = load_aal.dat.scl_slope;
        scl_inter = load_aal.dat.scl_inter;
        dtype = 'FLOAT32';
        offset = 0;
        dim = [91 109 91 150];
        descrip = 'run';
        fname = strcat('run','-',int2str(iirun_passive + iirun_jump),'.nii');
        input_data = Passive(irun).run; 
        real_save_image;
        %extract_time_series(Passive(irun).run,AAL_img,iirun_passive,'passive');
        
        iirun_rest = iirun_rest + 1;
        nifti_file = load_aal;
        offset = load_aal.dat.offset;
        scl_slope = load_aal.dat.scl_slope;
        scl_inter = load_aal.dat.scl_inter;
        dtype = 'FLOAT32';
        offset = 0;
        dim = [91 109 91 150];
        descrip = 'run';
        fname = strcat('run','-',int2str(iirun_rest + iirun_jump*2),'.nii');
        input_data = RestingState(irun).run; 
        real_save_image;
        %extract_time_series(RestingState(irun).run,AAL_img,iirun_rest,'rest');

    end

end

end

function extract_time_series(run,AAL_img,idx_run,condition)

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

mkdir(strcat('run-',condition,'-',int2str(idx_run)));

for iROI=1:nNodes
    
   idx_ROI = AAL_ROI(iROI).ID;
    
   idx_voxels = find(AAL_img==idx_ROI);
   
   for iVoxel=1:length(idx_voxels)
      
       [x,y,z] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
       
       voxel = squeeze(run(x,y,z,:));
       
       dlmwrite(strcat('run-',condition,'-',int2str(idx_run),'/','ROI-',int2str(iROI),'-','voxel-',int2str(iVoxel),'.txt'),voxel);
       
   end
    
end

end

function gen_masks

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

for iROI=1:nNodes
    
   idx_ROI = AAL_ROI(iROI).ID;
    
   idx_voxels = find(AAL_img==idx_ROI);
   
   AAL_ROI_outside = AAL_img;
   AAL_ROI_inside = AAL_img;
   
   AAL_ROI_outside(idx_voxels) = 0;
   AAL_ROI_outside(find(AAL_ROI_outside)) = 1;
   
   AAL_ROI_inside = zeros(size(AAL_img));
   
   for iVoxel=1:length(idx_voxels)
      
       [x,y,z] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

       AAL_ROI_inside(x,y,z) = iVoxel;
       
   end
   
   nifti_file = load_aal;
   offset = load_aal.dat.offset;
   scl_slope = load_aal.dat.scl_slope;
   scl_inter = load_aal.dat.scl_inter;
   dtype = 'FLOAT32';
   offset = 0;
   dim = [91 109 91];
   descrip = 'run';
   fname = strcat('ROI','-',int2str(iROI),'-inside','.nii');
   input_data = AAL_ROI_inside; 
   real_save_image;
   
   nifti_file = load_aal;
   offset = load_aal.dat.offset;
   scl_slope = load_aal.dat.scl_slope;
   scl_inter = load_aal.dat.scl_inter;
   dtype = 'FLOAT32';
   offset = 0;
   dim = [91 109 91];
   descrip = 'run';
   fname = strcat('ROI','-',int2str(iROI),'-outside','.nii');
   input_data = AAL_ROI_outside; 
   real_save_image;
        
   end
   
end

function gen_nVoxels

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

for iROI=1:nNodes
    
   idx_ROI = AAL_ROI(iROI).ID;
    
   idx_voxels = find(AAL_img==idx_ROI);
   
   nVoxels = length(idx_voxels);
   
   dlmwrite(strcat('ROI-',int2str(iROI),'.txt'),nVoxels);
   
end


end

function gen_global_Voxels


nRuns = 4;
nSubjects = 8;
nConditions = 3;

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

iGlobalVoxel = 0;

condition_label{1} = 'track';
condition_label{2} = 'passive';
condition_label{3} = 'rest';

for iCondition=1:nConditions
    
    iiRun = 0;
    
   for iSubject=1:nSubjects
      
      for iRun=1:nRuns
          
          iiRun = iiRun + 1;
          
         for iROI=1:nNodes
             
              idx_ROI = AAL_ROI(iROI).ID;
    
              idx_voxels = find(AAL_img==idx_ROI);
   
              nVoxels = length(idx_voxels);
              
              for iVoxel=1:nVoxels
             
                 iGlobalVoxel = iGlobalVoxel + 1;
                 
                 run{iGlobalVoxel} = strcat('run','-',condition_label{iCondition},'-',int2str(iiRun));
                 
                 roi(iGlobalVoxel) = iROI;
                 
                 idx_all_voxels(iGlobalVoxel) = iVoxel;
                  
              end
   
         end
          
          
      end
       
   end
    
end

% fid = fopen(strcat('iVoxel','-run','.txt'),'wt');
% fprintf(fid, run);
% fclose(fid);

dlmwrite(strcat('iVoxel','-ROI','.txt'),roi);
dlmwrite(strcat('iVoxel','-iVoxel','.txt'),idx_all_voxels);

end

function gen_global_Voxels_v2

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

iiVoxel = 0;
 for iROI=1:nNodes
             
    idx_ROI = AAL_ROI(iROI).ID;

    idx_voxels = find(AAL_img==idx_ROI);

    nVoxels = length(idx_voxels);

    for iVoxel=1:nVoxels
        
        iiVoxel = iiVoxel + 1;
        
        roi(iiVoxel) = iROI;
        voxel(iiVoxel) = iVoxel;

    end
              
 end

% dlmwrite(strcat('iVoxel','-ROI','.txt'),roi,'newline','pc');
% dlmwrite(strcat('iVoxel','-iVoxel','.txt'),voxel,'newline','pc');

fid = fopen(strcat('iVoxel','-ROI','.txt'),'wt');
fprintf(fid,'%d\n',roi);
fclose(fid);

fid = fopen(strcat('iVoxel','-iVoxel','.txt'),'wt');
fprintf(fid,'%d\n',voxel);
fclose(fid);

end


