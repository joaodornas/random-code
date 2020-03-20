
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_labels = load('ROI_MNI_V4_List.mat');
ROI = load_labels.ROI;

nROI = length(ROI);

for iROI=1:nROI
    
    IDs(iROI) = ROI(1,iROI).ID;
    
end

transformed_IDs = floor(IDs./100);

[ID_areas, nareas] = count_unique(transformed_IDs);

nAreas = length(ID_areas);

col = colormap(jet);

mi = 1;
ma = nAreas;

idx_areas = 1:nAreas;

%C = squeeze(ind2rgb(floor(((idx_areas(:)-mi)/(ma-mi))*size(col,1)),col));
C = idx_areas;

nVoxels = size(AAL_img,1)*size(AAL_img,2)*size(AAL_img,3);

for iVoxel=1:nVoxels
   
    ID_voxel = AAL_img(iVoxel);
    
    if ID_voxel ~= 0
        
        idx_ID_area_voxel = find(ismember(ID_areas,floor(ID_voxel/100)));
    
        AAL_img(iVoxel) = C(idx_ID_area_voxel);
        
    end
        
end

fname = strcat('AAL-RGB','.nii');
scl_slope = 1;
scl_inter = 0; 
dim = [size(AAL_img,1),size(AAL_img,2),size(AAL_img,3)];
dtype = 'FLOAT32';
offset = 0;
descrip = 'AAL - RGB';
nifti_file.mat = load_aal.mat;
nifti_file.mat_intent = load_aal.mat_intent;
nifti_file.mat0 = load_aal.mat0;
nifti_file.mat0_intent = load_aal.mat0_intent;

input_data = AAL_img; 

lowhigh_save_image;