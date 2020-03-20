function processed_voxel  = preprocessVoxel(voxel)

    voxel_detrended = detrend(voxel);

    processed_voxel = voxel_detrended - mean(voxel_detrended);

end
