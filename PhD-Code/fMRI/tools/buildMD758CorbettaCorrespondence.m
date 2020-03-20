
function buildMD758CorbettaCorrespondence

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

MNI_size = [91 109 91];

iiSeed = 0;

for iNet=1:length(all_networks_labels)
    
    network_label = all_networks_labels{iNet};
    
    seeds = getFunctionalSeeds_v6(network_label);
    
    for iSeed=1:length(seeds.ROI)
        
        iiSeed = iiSeed + 1;
        
        Nets_And_Seeds{iiSeed,1} = seeds.ROI(iSeed).label;
       
        x_mm = seeds.ROI(iSeed).x;
        y_mm = seeds.ROI(iSeed).y;
        z_mm = seeds.ROI(iSeed).z;
        
        MNI_x_center = 45;
        MNI_y_center = 63;
        MNI_z_center = 36;

        size_voxels_mm = 2;

        x_center = MNI_x_center + round(x_mm/size_voxels_mm) + 1;
        y_center = MNI_y_center + round(y_mm/size_voxels_mm) + 1;
        z_center = MNI_z_center + round(z_mm/size_voxels_mm) + 1;
        
        Nets_And_Seeds{iiSeed,2} = sub2ind(MNI_size,x_center,y_center,z_center);
        
        all_seeds_center(iiSeed) = sub2ind(MNI_size,x_center,y_center,z_center);
        
    end
    
end

load('Z:\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info.mat')

nROI = 90;

for iROI=1:nROI
    
    nClusters = length(ROI(iROI).clusters);
    
    for iCluster=1:nClusters
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
        
        fall_in_these_nets = ismember(all_seeds_center,idx_voxels);
        
        idx_fall = find(fall_in_these_nets);
        
        if ~isempty(idx_fall)
           
            for iDX=1:length(idx_fall)
               
                ROI(iROI).clusters(iCluster).Nets{iDX} = Nets_And_Seeds{idx_fall(iDX),1};
                
            end
            
        end
        
    end
    
end

ROI_MD758Corbetta = ROI;

save('ROI_MD758Corbetta.mat','ROI_MD758Corbetta');

end

