function save_ROI_aal_img( ROIs, fname, nifti_file, dtype, offset, scl_slope, scl_inter, descrip)

AAL_img = nifti('ROI_MNI_V4.nii');
AAL_dat = AAL_img.dat;
AAL_dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',AAL_dat.fname);
AAL_matrix = AAL_dat(:,:,:);

AAL_labels = load('ROI_MNI_V4_List.mat');

nROI = length(AAL_labels.ROI);

whole_brain = zeros(size(AAL_matrix));

for iROI=1:nROI
   
    RGB = AAL_labels.ROI(iROI).ID;
    
    idx = find(AAL_matrix == RGB);
    
    whole_brain(idx) = ROIs(iROI);
    
end

dim = [size(whole_brain)];

input_data = whole_brain; 
        
lowhigh_save_image;

end

