
%%% GET SUBJECTS SETTINGS

all_settings = getAllSettings;
nSubjects = length(all_settings);

%%% LOAD AAL INFO

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nAllVoxels = length(find(AAL_img));

nRuns = 4;

nROI = 90;

neg_vol = zeros(size(AAL_img));
pos_vol = zeros(size(AAL_img));

for iROI=1:nROI
   
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(strcat(int2str(iROI),':',label_ROI{iROI},':',int2str(length(idx_voxels)),':voxels',':',datestr(now)));
    
    nVoxels = length(idx_voxels);
    
    for iVoxel=1:nVoxels
       
        if mod(iVoxel,1000) == 0
            
            disp(strcat('iVoxel:',int2str(iVoxel),':',datestr(now)));
        
        end
        
        it = 0;
        
        for isubject=1:1
            
            settings = all_settings(isubject).settings;
            
            %%% TRACK
            load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-density-Track-',area_label,'.mat'));
            
            for irun=1:1
            
                it = it + 1;
                samples_track_pos(it) = FC_Density.run(irun).voxels_pos_density(iVoxel);
                samples_track_neg(it) = FC_Density.run(irun).voxels_neg_density(iVoxel);
            
            end
            
            clear FC_Density
            
        end
        
        mean_track_pos = mean(samples_track_pos);
        mean_track_neg = mean(samples_track_neg);
        
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
        
        neg_vol(idxx,idxy,idxz) = mean_track_neg;
        pos_vol(idxx,idxy,idxz) = mean_track_pos;
        
    end
    
end

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = load_aal.dat.dim;

descrip = 'paint';

fname = strcat('local-mean-pos-vol','.nii');
input_data = pos_vol; 
real_save_image;

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = load_aal.dat.dim;

descrip = 'paint';

fname = strcat('local-mean-neg-vol','.nii');
input_data = neg_vol; 
real_save_image;

