
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
setenv('FSL_DIR','/usr/local/fsl');

%system('sudo csh /usr/local/fsl/etc/fslconf/fsl.csh');

nClusters = 758;

for iCluster=1:nClusters
    
    system(sprintf('/usr/local/fsl/bin/fslmaths FuncClusterParcels2nodif.nii.gz -thr %s -uthr %s cluster-%s.nii.gz',int2str(iCluster),int2str(iCluster),int2str(iCluster)));
    
end

