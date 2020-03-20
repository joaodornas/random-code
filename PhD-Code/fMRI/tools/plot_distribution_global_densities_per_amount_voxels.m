
load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\GLOBAL\densities\outside-volumes\2nd\allVoxelsDensities.mat');

nVoxels = 160990;
nRuns = 32;

% track_pos = [];
% track_neg = [];
% passive_pos = [];
% passive_neg = [];
% rest_pos = [];
% rest_neg = [];
% for iVoxel=1:nVoxels
%     
%    track_pos = [track_pos, allVoxelsDensities.track_pos(iVoxel,1:end-1)]; 
%    track_neg = [track_neg, allVoxelsDensities.track_neg(iVoxel,1:end-1)]; 
%    passive_pos = [passive_pos, allVoxelsDensities.passive_pos(iVoxel,1:end-1)]; 
%    passive_neg = [passive_neg, allVoxelsDensities.passive_neg(iVoxel,1:end-1)]; 
%    rest_pos = [rest_pos, allVoxelsDensities.rest_pos(iVoxel,1:end-1)]; 
%    rest_neg = [rest_neg, allVoxelsDensities.rest_neg(iVoxel,1:end-1)]; 
%     
% end

track_pos = reshape(allVoxelsDensities.track_pos(:,1:end-1),1,nVoxels*nRuns);
track_neg = reshape(allVoxelsDensities.track_neg(:,1:end-1),1,nVoxels*nRuns);
passive_pos = reshape(allVoxelsDensities.passive_pos(:,1:end-1),1,nVoxels*nRuns);
passive_neg = reshape(allVoxelsDensities.passive_neg(:,1:end-1),1,nVoxels*nRuns);
rest_pos = reshape(allVoxelsDensities.rest_pos(:,1:end-1),1,nVoxels*nRuns);
rest_neg = reshape(allVoxelsDensities.rest_neg(:,1:end-1),1,nVoxels*nRuns);

track_pos(isnan(track_pos)) = 0;
track_neg(isnan(track_neg)) = 0;
passive_pos(isnan(passive_pos)) = 0;
passive_neg(isnan(passive_neg)) = 0;
rest_pos(isnan(rest_pos)) = 0;
rest_neg(isnan(rest_neg)) = 0;

%%% PLOT DISTRIBUTION DENSITY

% min_val = min([track_pos(:);passive_pos(:)]);
% max_val = max([track_pos(:);passive_pos(:)]);
% plotDistributionDensity(min_val,max_val,track_pos,passive_pos,'Track-POS','Passive-POS','Track-Passive-POS','density','# of voxels');
% 
% min_val = min([track_neg(:);passive_neg(:)]);
% max_val = max([track_neg(:);passive_neg(:)]);
% plotDistributionDensity(min_val,max_val,track_neg,passive_neg,'Track-NEG','Passive-NEG','Track-Passive-NEG','density','# of voxels');
% 
% min_val = min([passive_pos(:);rest_pos(:)]);
% max_val = max([passive_pos(:);rest_pos(:)]);
% plotDistributionDensity(min_val,max_val,passive_pos,rest_pos,'Passive-POS','Rest-POS','Passive-Rest-POS','density','# of voxels');
% 
% min_val = min([passive_neg(:);rest_neg(:)]);
% max_val = max([passive_neg(:);rest_neg(:)]);
% plotDistributionDensity(min_val,max_val,passive_neg,rest_neg,'Passive-NEG','Rest-NEG','Passive-Rest-NEG','density','# of voxels');

%%% PLOT CUMULATIVE DISTRIBUTION with KOLMOGOROV-SMIRNOV TEST

% plotCumulativeDistributionKSTest2(track_pos,passive_pos,'Track-POS','Passive-POS','Track-Passive-POS','density','Cumulative Probability');
% 
% plotCumulativeDistributionKSTest2(track_neg,passive_neg,'Track-NEG','Passive-NEG','Track-Passive-NEG','density','Cumulative Probability');
% 
% plotCumulativeDistributionKSTest2(passive_pos,rest_pos,'Passive-POS','Rest-POS','Passive-Rest-POS','density','Cumulative Probability');
% 
% plotCumulativeDistributionKSTest2(passive_neg,rest_neg,'Passive-NEG','Rest-NEG','Passive-Rest-NEG','density','Cumulative Probability');


%%% PLOT CUMULATIVE DISTRIBUTION with KOLMOGOROV-SMIRNOV TEST with FUNCTIONAL PARCELLATION



