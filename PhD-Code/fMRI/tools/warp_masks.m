
path_ref = '/Volumes/dropbox/_DATA/LOW-HIGH-ATTENTION/SUBJECT-1-22-10-2015/preprocessed/B0-DTI/4.BET';
path_mat = '/Volumes/dropbox/_DATA/LOW-HIGH-ATTENTION/SUBJECT-1-22-10-2015/preprocessed/B0-DTI/7.FLIRT';

nClusters = 758;

for iCluster=1:nClusters
    
    %system(sprintf('/usr/local/fsl/bin/applywarp -i Cluster-%s-inside.nii -r %s/nodif_brain.nii.gz -o Cluster-%s-inside2nodif -w %s/FuncClusterParcels2nodif.nii.gz',int2str(iCluster),path_nodif,int2str(iCluster),path_warp));
    %system(sprintf('/usr/local/fsl/bin/applywarp -i Cluster-%s-outside.nii -r %s/nodif_brain.nii.gz -o Cluster-%s-outside2nodif -w %s/FuncClusterParcels2nodif.nii.gz',int2str(iCluster),path_nodif,int2str(iCluster),path_warp));
    
    
    system(sprintf('/usr/local/fsl/bin/flirt -in /Volumes/dropbox/__tmp/Func_masks/Cluster-%s-inside.nii -applyxfm -init %s/std2nodif.mat -out /Volumes/dropbox/__tmp/Cluster-%s-inside2nodif -paddingsize 0.0 -interp nearestneighbour -ref %s/nodif_brain.nii.gz',int2str(iCluster),path_mat,int2str(iCluster),path_ref));
    system(sprintf('/usr/local/fsl/bin/flirt -in /Volumes/dropbox/__tmp/Func_masks/Cluster-%s-outside.nii -applyxfm -init %s/std2nodif.mat -out /Volumes/dropbox/__tmp/Cluster-%s-outside2nodif -paddingsize 0.0 -interp nearestneighbour -ref %s/nodif_brain.nii.gz',int2str(iCluster),path_mat,int2str(iCluster),path_ref));
    
end