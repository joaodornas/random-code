function seeds_overlap = seeIfAVoxelIsInsideASeed(Voxels)

seeds_overlap = cell.empty;
ioverlap = 0;

MNI_x_center = 45;
MNI_y_center = 63;
MNI_z_center = 36;

size_voxels_mm = 2;
ROI_radius_mm = 6 - size_voxels_mm;
ROI_radius_voxels = ROI_radius_mm/size_voxels_mm;

all_networks = {'DAN' 'VAN' 'SMN' 'VIS' 'FPC' 'LAN' 'DMN' 'AUD'};

for iNet=1:length(all_networks)
   
    seeds = getFunctionalSeeds_v5(all_networks{iNet});
        
    for iROI=1:length(seeds.ROI)

        x_center = MNI_x_center + round(seeds.ROI(iROI).x/size_voxels_mm) + 1;
        y_center = MNI_y_center + round(seeds.ROI(iROI).y/size_voxels_mm) + 1;
        z_center = MNI_z_center + round(seeds.ROI(iROI).z/size_voxels_mm) + 1;

        %xgv = (x_center-ROI_radius_voxels):(x_center+ROI_radius_voxels);
        %ygv = (y_center-ROI_radius_voxels):(y_center+ROI_radius_voxels);
        %zgv = (z_center-ROI_radius_voxels):(z_center+ROI_radius_voxels);

        %[X,Y,Z] = meshgrid(xgv,ygv,zgv);

        %seed_mesh = [X(:),Y(:),Z(:)];

        % x = Voxels(:,1);
        % y = Voxels(:,2);
        % z = Voxels(:,3);

        if ismember([x_center y_center z_center],Voxels,'rows');

            if isempty(seeds_overlap)
                ioverlap = ioverlap + 1;
                seeds_overlap{ioverlap} = seeds.ROI(iROI).label;
            else
                ioverlap = ioverlap + 1;
                seeds_overlap{ioverlap} = seeds.ROI(iROI).label;
            end
            
        end

    end
    
end

end
