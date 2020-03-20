function real_FC_voxel_AAL_ROI_densities_change_LUT


filename{1} = 'LHR-Correlation-Contrast-Passive-All-z-clu-both-fl-2-mm';
filename{2} = 'LHR-Correlation-Contrast-Passive-Neg-z-clu-both-fl-2-mm';
filename{3} = 'LHR-Correlation-Contrast-Passive-Pos-z-clu-both-fl-2-mm';
filename{4} = 'LHR-Correlation-Contrast-Track-All-z-clu-both-fl-2-mm';
filename{5} = 'LHR-Correlation-Contrast-Track-Neg-z-clu-both-fl-2-mm';
filename{6} = 'LHR-Correlation-Contrast-Track-Pos-z-clu-both-fl-2-mm';

for iFile=1:length(filename)
    
    changeLUT(filename{iFile});
    
end

end

function changeLUT(filename)

load_img = nifti(strcat(filename,'.nii'));
folder = '/Volumes/dropbox/_DATA/LOW-HIGH-ATTENTION/all-subjects/Real/FC_Voxels_AAL_ROI/densities/rendering/both-flipped/clusters-dlh=1.5-both-fl-2mm-FS_spm_CanonicalBrain_norecon-unzip';
load_img.dat.fname = strcat(folder,'/',load_img.dat.fname);

my_img = load_img.dat(:,:,:);

idx_negative = find(my_img<0);
idx_positive = find(my_img>0);

new_img = uint8(zeros(size(my_img)));

nVoxelsPositive = length(idx_positive);
nVoxelsNegative = length(idx_negative);

LUTP = 23:42;
LUTN = 22:-1:1;

max_z = 4.5*10;
min_z = 1.6*10;

iiVoxel = 0;

for iVoxel=1:nVoxelsPositive
   
    [idxx,idxy,idxz] = ind2sub(size(my_img),idx_positive(iVoxel));
    
    voxel = floor(abs(my_img(idxx,idxy,idxz)) * 10);
    
    if (voxel >= min_z) & (voxel <= max_z)
        
        iiVoxel = iiVoxel + 1;
    
        new_voxel = round(((voxel-min_z)/(max_z-min_z))*(max(LUTP)-min(LUTP))) + min(LUTP);
    
        all_positives(iiVoxel) = new_voxel;
    
        new_img(idxx,idxy,idxz) = uint8(new_voxel);
    
    else
        
        new_img(idxx,idxy,idxz) = 0;
        
    end
    
end

iiVoxel = 0;

for iVoxel=1:nVoxelsNegative
   
    [idxx,idxy,idxz] = ind2sub(size(my_img),idx_negative(iVoxel));
    
    voxel = floor(abs(my_img(idxx,idxy,idxz)) * 10);
    
    if (voxel >= min_z) & (voxel <= max_z)
        
        iiVoxel = iiVoxel + 1;
    
        new_voxel = round(((voxel-min_z)/(max_z-min_z))*(max(LUTN)-min(LUTN))) + min(LUTN);
    
        all_negatives(iiVoxel) = new_voxel;
    
        new_img(idxx,idxy,idxz) = uint8(new_voxel);
    
    else
        
        new_img(idxx,idxy,idxz) = 0;
        
    end
    
end

nifti_file = load_img;
offset = load_img.dat.offset;
scl_slope = load_img.dat.scl_slope;
scl_inter = load_img.dat.scl_inter;

dtype = 'uint8-le';
offset = 0;

dim = load_img.dat.dim;

descrip = 'results';

fname = strcat(strcat(filename,'-LUT','.nii'));
input_data = new_img; 
real_save_image;

end
