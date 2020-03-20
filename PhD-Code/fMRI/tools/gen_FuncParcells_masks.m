
function gen_FuncParcells_masks

filepath = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\ParcellationVolume';
volume = 'LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation';

load_aal = nifti(strcat(volume,'.nii'));
load_aal.dat.fname = strcat(filepath,'\',load_aal.dat.fname);

FuncParcells = load_aal.dat(:,:,:);

nClusters = 758;

for iCluster=1:nClusters
    
   idx_voxels = find(FuncParcells==iCluster);
   
   AAL_ROI_outside = FuncParcells;
   AAL_ROI_inside = FuncParcells;
   
   AAL_ROI_outside(idx_voxels) = 0;
   AAL_ROI_outside(find(AAL_ROI_outside)) = 1;
   
   AAL_ROI_inside = zeros(size(FuncParcells));
   
   for iVoxel=1:length(idx_voxels)
      
       [x,y,z] = ind2sub(size(FuncParcells),idx_voxels(iVoxel));

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
   fname = strcat('Cluster','-',int2str(iCluster),'-inside','.nii');
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
   fname = strcat('Cluster','-',int2str(iCluster),'-outside','.nii');
   input_data = AAL_ROI_outside; 
   real_save_image;
        
   end
   
end
