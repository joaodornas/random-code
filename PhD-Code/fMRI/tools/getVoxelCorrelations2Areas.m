function  [rho, pval] = getVoxelCorrelations2Areas(RUN,AAL_img,idx_AAL_area1,idx_AAL_area2)

idx_voxels_area1 = find(AAL_img==idx_AAL_area1);
idx_voxels_area2 = find(AAL_img==idx_AAL_area2);

nVoxels_area1 = length(idx_voxels_area1);
nVoxels_area2 = length(idx_voxels_area2);
nTR = size(RUN,4);

disp(strcat('nVoxels:',int2str(nVoxels_area1+nVoxels_area2)));

rho = zeros(nVoxels_area1,nVoxels_area2);
pval = zeros(nVoxels_area1,nVoxels_area2);

area1 = zeros(nVoxels_area1,nTR);
area2 = zeros(nVoxels_area2,nTR);

for iVoxel=1:nVoxels_area1
    
    [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels_area1(iVoxel));
    
    area1(iVoxel,:) = RUN(idxx,idxy,idxz,:);

end

for iVoxel=1:nVoxels_area2
    
    [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels_area2(iVoxel));
    
    area2(iVoxel,:) = RUN(idxx,idxy,idxz,:);

end

tic
for iVoxel=1:nVoxels_area1
    
    if mod(iVoxel,500) == 0; disp(strcat('Voxels so far:',int2str(iVoxel))); end
   
%    [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
    
%    voxel_X = RUN(idxx,idxy,idxz,:);
 
    voxel_X = area1(iVoxel,:);
    
    rho_col = [];
    pval_col = [];
    
    parfor iiVoxel=1:nVoxels_area2
        
%        [iidxx,iidxy,iidxz] = ind2sub(size(AAL_img),idx_voxels(columns(iiVoxel)));
        
%        voxel_Y = RUN(iidxx,iidxy,iidxz,:);
 
        voxel_Y = area2(iiVoxel,:);

        [rho_mat,pval_mat] = corrcoef(voxel_X,voxel_Y);
        
        rho(iVoxel,iiVoxel) = rho_mat(1,2);
        pval(iVoxel,iiVoxel) = pval_mat(1,2);
        
    end
    
end
toc

end


