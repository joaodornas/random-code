function ROIname = getROIName(X,Y,Z)

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('/Users/joaodornas/Dropbox/_Research/_toolBox/aal_for_SPM8/',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

[x_s, y_s, z_s] = size(AAL_img);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

MNI_x_center = 45;
MNI_y_center = 63;
MNI_z_center = 36;

size_voxels_mm = 2;

x_center = MNI_x_center - round(X/size_voxels_mm);
y_center = MNI_y_center + round(Y/size_voxels_mm);
z_center = MNI_z_center + round(Z/size_voxels_mm);

idx_ROI = AAL_img(x_center + 1,y_center + 1,z_center + 1);

ID = vertcat(AAL_ROI.ID);

idx_ROI_AAL = find(ID == idx_ROI);

%disp(strcat('ROI:',string(AAL_ROI(idx_ROI_AAL).Nom_L)));

ROIname = string(AAL_ROI(idx_ROI_AAL).Nom_L);

end

