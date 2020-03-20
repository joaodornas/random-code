clear all

filepath = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\ParcellationVolume';
volume = 'LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation';

data = nifti(strcat(volume,'.nii'));
data.dat.fname = strcat(filepath,'\',data.dat.fname);

FuncParcells = data.dat(:,:,:);

nClusters = 758;

clusters_stats = zeros(nClusters,7);

iiVoxel = 0;
for iCluster=1:nClusters
    
   idx_voxels = find(FuncParcells == iCluster);
   
   this_cluster = zeros(1,7);
   
   for iVoxel=1:length(idx_voxels)
       
       iiVoxel = iiVoxel + 1;
       
        [x,y,z] = ind2sub(size(FuncParcells),idx_voxels(iVoxel));
        
        count = 0;
        
        if FuncParcells(x+1,y,z) == iCluster; count = count + 1; end
        if FuncParcells(x,y+1,z) == iCluster; count = count + 1; end
        if FuncParcells(x,y,z+1) == iCluster; count = count + 1; end
        if FuncParcells(x-1,y,z) == iCluster; count = count + 1; end
        if FuncParcells(x,y-1,z) == iCluster; count = count + 1; end
        if FuncParcells(x,y,z-1) == iCluster; count = count + 1; end
       
        this_cluster(count+1) = this_cluster(count+1) + 1;
        
   end
   
   clusters_stats(iCluster,:) = this_cluster(:);
   
   clear this_cluster
   
end

stat_sum = sum(clusters_stats,1)./iiVoxel;

bar(stat_sum,'k');
xlabel('# of voxels as neighbors');
ylabel('% of voxels');