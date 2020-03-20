clear all

cluster_file = '758-Clusters.nii';
cluster_folder = '/Volumes/dropbox/_DATA/Parcellation/758-Cluster';

load_cluster = nifti(strcat(cluster_folder,'/',cluster_file));
load_cluster.dat.fname = strcat(cluster_folder,'/',cluster_file);

cluster_img = load_cluster.dat(:,:,:);

new_cluster_img = zeros(size(cluster_img));

load('/Volumes/dropbox/_DATA/LOW-HIGH-ATTENTION/all-subjects/Real/FC_Voxels_AAL_ROI/FC-Voxels-AAL-ROI-corr-KMeans/FC-Voxels-AAL-ROI-corr-KMeans-Info.mat')

nROI = 90;

iiCluster = 0;

for iROI=1:nROI
    
    nClusters = length(ROI(iROI).clusters);
    
    for iCluster=1:nClusters
        
        iiCluster = iiCluster + 1;
        
        idx_voxels = find(cluster_img==iiCluster);
        
        nVoxels = length(idx_voxels);
        
        for iVoxel=1:nVoxels
            
            new_cluster_img(idx_voxels(iVoxel)) = iCluster;
            
        end
        
    end
    
end

scl_slope = 1;
scl_inter = 0;
dim = [size(cluster_img,1),size(cluster_img,2),size(cluster_img,3)];
dtype = 'FLOAT32';
offset = 0;
nifti_file.mat = load_cluster.mat;
nifti_file.mat_intent = load_cluster.mat_intent;
nifti_file.mat0 = load_cluster.mat0;
nifti_file.mat0_intent = load_cluster.mat0_intent;

fname = '758-clusters-same-index.nii';
descrip = '758-clusters-same-index';
input_data = new_cluster_img;
real_save_image;



