function  [area_rho, area_pval] = getVoxelCorrelations_corr(RUN,AAL_img,idx_AAL)

idx_voxels = find(AAL_img==idx_AAL);

nVoxels = length(idx_voxels);
nTR = size(RUN,4);

disp(strcat('nVoxels:',int2str(nVoxels)));

area_rho = zeros(nVoxels,nVoxels);
area_pval = zeros(nVoxels,nVoxels);

area = zeros(nVoxels,nTR);

for iVoxel=1:nVoxels
    
    [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
    
    area(iVoxel,:) = squeeze(RUN(idxx,idxy,idxz,:));

end

tic
[area_rho, area_pval] = corr(area');
toc

% tic
% for iVoxel=1:nVoxels
%     
%     if mod(iVoxel,500) == 0; disp(strcat('Voxels so far:',int2str(iVoxel))); end
%    
% %    [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
%     
% %    voxel_X = RUN(idxx,idxy,idxz,:);
%  
%     voxel_X = area(iVoxel,:);
%     
%     columns = iVoxel:nVoxels;
%     
%     rho_col = [];
%     pval_col = [];
%     
%     for iiVoxel=1:length(columns)
%         
% %        [iidxx,iidxy,iidxz] = ind2sub(size(AAL_img),idx_voxels(columns(iiVoxel)));
%         
% %        voxel_Y = RUN(iidxx,iidxy,iidxz,:);
%  
%         voxel_Y = area(columns(iiVoxel),:);
% 
%         [rho,pval] = corrcoef(voxel_X,voxel_Y);
%         
%         rho_col(iiVoxel) = rho(1,2);
%         pval_col(iiVoxel) = pval(1,2);
%         
%     end
%     
%     area_rho(iVoxel,iVoxel:nVoxels) = rho_col;
%     area_pval(iVoxel,iVoxel:nVoxels) = pval_col;
%         
%     area_rho(iVoxel:nVoxels,iVoxel) = rho_col';
%     area_pval(iVoxel:nVoxels,iVoxel) = pval_col';
%     
% end
% toc

end


