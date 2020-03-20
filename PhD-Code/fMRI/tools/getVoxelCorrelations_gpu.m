function  [area_rho, area_pval] = getVoxelCorrelations_gpu(RUN,AAL_img,idx_AAL)

idx_voxels = find(AAL_img==idx_AAL);

nVoxels = length(idx_voxels);
nTR = size(RUN,4);

disp(strcat('nVoxels:',int2str(nVoxels)));

area_rho = zeros(nVoxels,nVoxels);
area_pval = zeros(nVoxels,nVoxels);

area = gpuArray(zeros(nVoxels,nTR));

for iVoxel=1:nVoxels
    
    [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
    
    area(iVoxel,:) = RUN(idxx,idxy,idxz,:);

end

area = area';
gpuArea = gpuArray(area);

gpuRHO = corr(gpuArea);
rho = gather(gpuRHO);            

t = rho./(sqrt((1-(rho.^2))./(nTR-2)));
pval = (1-tcdf(t,nTR-1));


area_rho = rho;
area_pval = pval;

end


