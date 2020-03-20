function netPerMD = getNetPerMD(network_label)

nTotalClusters = 758;

seeds = getFunctionalSeeds_v6(network_label);

func_vol = nifti('LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
func_vol.dat.fname = 'Z:\_DATA\Parcellation\758-Cluster\LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii';
func_clu = func_vol.dat(:,:,:);

MNI_x_center = 45;
MNI_y_center = 63;
MNI_z_center = 36;

size_voxels_mm = 2;
ROI_radius_mm = 6 - size_voxels_mm;
ROI_radius_voxels = ROI_radius_mm/size_voxels_mm;

netPerMD = zeros(1,nTotalClusters);

nSeeds = length(seeds.ROI);

for iSeed=1:nSeeds

    x_center = MNI_x_center + ceil(seeds.ROI(iSeed).x/size_voxels_mm) + 1;
    y_center = MNI_y_center + ceil(seeds.ROI(iSeed).y/size_voxels_mm) + 1;
    z_center = MNI_z_center + ceil(seeds.ROI(iSeed).z/size_voxels_mm) + 1;
    
    idx_clu = func_clu(x_center,y_center,z_center);
    
    if idx_clu == 0
        
        idx_clu = func_clu(x_center + 1,y_center,z_center);
        
    end
    
    if idx_clu == 0
        
        idx_clu = func_clu(x_center - 1,y_center,z_center);
        
    end
    
    if idx_clu == 0
        
        idx_clu = func_clu(x_center,y_center + 1,z_center);
        
    end
    
    if idx_clu == 0
        
        idx_clu = func_clu(x_center,y_center - 1,z_center);
        
    end
    
    if idx_clu == 0
        
        idx_clu = func_clu(x_center,y_center,z_center + 1);
        
    end
    
    if idx_clu == 0
        
        idx_clu = func_clu(x_center,y_center,z_center - 1);
        
    end
    
    if idx_clu ~= 0; netPerMD(idx_clu) = 1; end

end


end

