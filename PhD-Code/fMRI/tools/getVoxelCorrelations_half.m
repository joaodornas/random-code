function  [area_rho, area_pval] = getVoxelCorrelations_half(RUN,AAL_img,idx_AAL,which_half)

idx_voxels = find(AAL_img==idx_AAL);

nVoxels = length(idx_voxels);
nTR = size(RUN,4);

disp(strcat('nVoxels:',int2str(nVoxels)));

area_rho = zeros(nVoxels,nVoxels);
area_pval = zeros(nVoxels,nVoxels);

area = zeros(nVoxels,nTR/2);

for iVoxel=1:nVoxels
    
    [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
    
    if strcmp(which_half,'FH')
        
        area(iVoxel,:) = RUN(idxx,idxy,idxz,1:(nTR/2));
        
    elseif strcmp(which_half,'SH')
        
        area(iVoxel,:) = RUN(idxx,idxy,idxz,((nTR/2)+1):end);
        
    end

end

area = area';

[area_rho,area_pval] = corr(area);

end


