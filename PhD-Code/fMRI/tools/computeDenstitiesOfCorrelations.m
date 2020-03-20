function [dAll_mat, dAll_pos_mat, dAll_neg_mat, voxels_all_density,voxels_pos_density,voxels_neg_density] = computeDenstitiesOfCorrelations(rho_mat,pval_mat)

pcriterion = 0.01;

rho_mat( pval_mat > pcriterion ) = 0;

rho_mat( rho_mat == 1 ) = 0.99;

All_mat = rho_mat ~= 0;

All_pos_mat = rho_mat > 0;

All_neg_mat = rho_mat < 0;

dAll_mat = double(All_mat);
dAll_pos_mat = double(All_pos_mat);
dAll_neg_mat = double(All_neg_mat);

voxels_all_density = sum(dAll_mat,1)./size(rho_mat,1);
voxels_pos_density = sum(dAll_pos_mat,1)./size(rho_mat,1);
voxels_neg_density = sum(dAll_neg_mat,1)./size(rho_mat,1);

end

