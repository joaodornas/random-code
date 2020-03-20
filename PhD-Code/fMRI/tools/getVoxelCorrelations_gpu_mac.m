function  [area_rho, area_pval] = getVoxelCorrelations_gpu_mac(RUN,AAL_img,idx_AAL)

idx_voxels = find(AAL_img==idx_AAL);

nVoxels = length(idx_voxels);
nTR = size(RUN,4);

disp(strcat('nVoxels:',int2str(nVoxels)));

area_rho = zeros(nVoxels,nVoxels);
area_pval = zeros(nVoxels,nVoxels);

area = zeros(nVoxels,nTR);

for iVoxel=1:nVoxels
    
    [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
    
    area(iVoxel,:) = RUN(idxx,idxy,idxz,:);

end

gpu_area = gpuArray(area);

tic
for iVoxel=1:nVoxels
    
    if mod(iVoxel,500) == 0; disp(strcat('Voxels so far:',int2str(iVoxel))); end
   
%    [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
    
%    voxel_X = RUN(idxx,idxy,idxz,:);
 
    voxel_X = gpu_area(iVoxel,:);
    
    columns = iVoxel:nVoxels;
    
    rho_col = [];
    pval_col = [];
    
    parfor iiVoxel=1:length(columns)
        
%        [iidxx,iidxy,iidxz] = ind2sub(size(AAL_img),idx_voxels(columns(iiVoxel)));
        
%        voxel_Y = RUN(iidxx,iidxy,iidxz,:);
 
        voxel_Y = gpu_area(columns(iiVoxel),:);

        gpu_rho = corrcoef(voxel_X,voxel_Y)
        rho = gather(gpu_rho(1,2));
        
%         rho_col(iiVoxel) = gather(rho(1,2));
%         pval_col(iiVoxel) = gather(pval(1,2));

        %rho = corr2(voxel_X,voxel_Y);
            
        t = (rho*sqrt(nTR-2))/sqrt(1-rho^2);
        pval = 1-tcdf(t,nTR-1);
        %t = (gather(rho)*sqrt(nTR-2))/sqrt(1-gather(rho)^2);
        %pval = 1-tcdf(t,nTR-1);
        
        rho_col(iiVoxel) = rho;
	    pval_col(iiVoxel) = pval;
        
    end
    
    area_rho(iVoxel,iVoxel:nVoxels) = rho_col;
    area_pval(iVoxel,iVoxel:nVoxels) = pval_col;
        
    area_rho(iVoxel:nVoxels,iVoxel) = rho_col';
    area_pval(iVoxel:nVoxels,iVoxel) = pval_col';
    
end
toc

end


