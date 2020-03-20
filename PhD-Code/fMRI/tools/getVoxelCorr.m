function rho_mat = getVoxelCorr(all_voxels)


nVoxels = size(all_voxels,2);

rho_mat = zeros(nVoxels,nVoxels);
pval_mat = zeros(nVoxels,nVoxels);

for iVoxel=1:nVoxels
    
    voxel_time_series = all_voxels(:,iVoxel);

    for iiVoxel=iVoxel:nVoxels
        
        next_voxel_time_series = all_voxels(:,iiVoxel);
        
        rho = corr2(voxel_time_series,next_voxel_time_series);
        
        rho_mat(iVoxel,iiVoxel) = rho;
        
        rho_mat(iiVoxel,iVoxel) = rho;
        
    end
    
end

end

