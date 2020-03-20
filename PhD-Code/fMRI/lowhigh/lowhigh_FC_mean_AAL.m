function lowhigh_FC_mean_AAL

%%% SETTINGS

settings_jan_0805;
%settings_elena_2905;

doTheMath(settings);

%plotResults(settings);
%plotResultsHalf(settings);

%plotFPO(settings);
%plotFPOHalf(settings);

%rhoBtwRuns(settings);

end

function doTheMath(settings)

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

%%% LOAD DATA

lowhigh_load_all_data_FSL;

% [mean_area_MOT4Run1_voxel,mean_area_MOT2Run1_voxel,mean_area_RestingStateRun1_voxel] = getMeanParcellation(MOT4Run1,MOT2Run1,RestingStateRun1);
% [mean_area_MOT4Run2_voxel,mean_area_MOT2Run2_voxel,mean_area_RestingStateRun2_voxel] = getMeanParcellation(MOT4Run2,MOT2Run2,RestingStateRun2);

% [rho_MOT4Run1,pval_MOT4Run1,rho_MOT2Run1,pval_MOT2Run1,rho_RestingStateRun1,pval_RestingStateRun1] = getCorrelations(mean_area_MOT4Run1_voxel,mean_area_MOT2Run1_voxel,mean_area_RestingStateRun1_voxel);
% [rho_MOT4Run2,pval_MOT4Run2,rho_MOT2Run2,pval_MOT2Run2,rho_RestingStateRun2,pval_RestingStateRun2] = getCorrelations(mean_area_MOT4Run2_voxel,mean_area_MOT2Run2_voxel,mean_area_RestingStateRun2_voxel);
% 
% save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL','.mat'),'rho_MOT4Run1','pval_MOT4Run1','rho_MOT2Run1','pval_MOT2Run1','rho_RestingStateRun1','pval_RestingStateRun1','rho_MOT4Run2','pval_MOT4Run2','rho_MOT2Run2','pval_MOT2Run2','rho_RestingStateRun2','pval_RestingStateRun2','-v7.3');

% [rho_MOT4Run1,pval_MOT4Run1,rho_MOT2Run1,pval_MOT2Run1,rho_RestingStateRun1,pval_RestingStateRun1] = getCorrelations_corrcoef(mean_area_MOT4Run1_voxel,mean_area_MOT2Run1_voxel,mean_area_RestingStateRun1_voxel);
% [rho_MOT4Run2,pval_MOT4Run2,rho_MOT2Run2,pval_MOT2Run2,rho_RestingStateRun2,pval_RestingStateRun2] = getCorrelations_corrcoef(mean_area_MOT4Run2_voxel,mean_area_MOT2Run2_voxel,mean_area_RestingStateRun2_voxel);
% 
% save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL-corrcoef','.mat'),'rho_MOT4Run1','pval_MOT4Run1','rho_MOT2Run1','pval_MOT2Run1','rho_RestingStateRun1','pval_RestingStateRun1','rho_MOT4Run2','pval_MOT4Run2','rho_MOT2Run2','pval_MOT2Run2','rho_RestingStateRun2','pval_RestingStateRun2','-v7.3');

% [rho_MOT4Run1_FH,pval_MOT4Run1_FH,rho_MOT2Run1_FH,pval_MOT2Run1_FH,rho_RestingStateRun1_FH,pval_RestingStateRun1_FH, rho_MOT4Run1_SH,pval_MOT4Run1_SH,rho_MOT2Run1_SH,pval_MOT2Run1_SH,rho_RestingStateRun1_SH,pval_RestingStateRun1_SH] = getCorrelationsHalf(mean_area_MOT4Run1_voxel,mean_area_MOT2Run1_voxel,mean_area_RestingStateRun1_voxel);
% [rho_MOT4Run2_FH,pval_MOT4Run2_FH,rho_MOT2Run2_FH,pval_MOT2Run2_FH,rho_RestingStateRun2_FH,pval_RestingStateRun2_FH, rho_MOT4Run2_SH,pval_MOT4Run2_SH,rho_MOT2Run2_SH,pval_MOT2Run2_SH,rho_RestingStateRun2_SH,pval_RestingStateRun2_SH] = getCorrelationsHalf(mean_area_MOT4Run2_voxel,mean_area_MOT2Run2_voxel,mean_area_RestingStateRun2_voxel);
% 
% save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL-Half','.mat'),'rho_MOT4Run1_FH','pval_MOT4Run1_FH','rho_MOT2Run1_FH','pval_MOT2Run1_FH','rho_RestingStateRun1_FH','pval_RestingStateRun1_FH','rho_MOT4Run1_SH','pval_MOT4Run1_SH','rho_MOT2Run1_SH','pval_MOT2Run1_SH','rho_RestingStateRun1_SH','pval_RestingStateRun1_SH','rho_MOT4Run2_FH','pval_MOT4Run2_FH','rho_MOT2Run2_FH','pval_MOT2Run2_FH','rho_RestingStateRun2_FH','pval_RestingStateRun2_FH','rho_MOT4Run2_SH','pval_MOT4Run2_SH','rho_MOT2Run2_SH','pval_MOT2Run2_SH','rho_RestingStateRun2_SH','pval_RestingStateRun2_SH','-v7.3');

[mean_area_MOT4Run1_voxel_FPO,mean_area_MOT2Run1_voxel_FPO,mean_area_RestingStateRun1_voxel_FPO] = getMeanParcellationFPO(MOT4Run1,MOT2Run1,RestingStateRun1);
[mean_area_MOT4Run2_voxel_FPO,mean_area_MOT2Run2_voxel_FPO,mean_area_RestingStateRun2_voxel_FPO] = getMeanParcellationFPO(MOT4Run2,MOT2Run2,RestingStateRun2);

[rho_MOT4Run1_FH,pval_MOT4Run1_FH,rho_MOT2Run1_FH,pval_MOT2Run1_FH,rho_RestingStateRun1_FH,pval_RestingStateRun1_FH, rho_MOT4Run1_SH,pval_MOT4Run1_SH,rho_MOT2Run1_SH,pval_MOT2Run1_SH,rho_RestingStateRun1_SH,pval_RestingStateRun1_SH] = getCorrelationsHalfFPO(mean_area_MOT4Run1_voxel_FPO,mean_area_MOT2Run1_voxel_FPO,mean_area_RestingStateRun1_voxel_FPO);
[rho_MOT4Run2_FH,pval_MOT4Run2_FH,rho_MOT2Run2_FH,pval_MOT2Run2_FH,rho_RestingStateRun2_FH,pval_RestingStateRun2_FH, rho_MOT4Run2_SH,pval_MOT4Run2_SH,rho_MOT2Run2_SH,pval_MOT2Run2_SH,rho_RestingStateRun2_SH,pval_RestingStateRun2_SH] = getCorrelationsHalfFPO(mean_area_MOT4Run2_voxel_FPO,mean_area_MOT2Run2_voxel_FPO,mean_area_RestingStateRun2_voxel_FPO);

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL-FPO-Half','.mat'),'rho_MOT4Run1_FH','pval_MOT4Run1_FH','rho_MOT2Run1_FH','pval_MOT2Run1_FH','rho_RestingStateRun1_FH','pval_RestingStateRun1_FH','rho_MOT4Run1_SH','pval_MOT4Run1_SH','rho_MOT2Run1_SH','pval_MOT2Run1_SH','rho_RestingStateRun1_SH','pval_RestingStateRun1_SH','rho_MOT4Run2_FH','pval_MOT4Run2_FH','rho_MOT2Run2_FH','pval_MOT2Run2_FH','rho_RestingStateRun2_FH','pval_RestingStateRun2_FH','rho_MOT4Run2_SH','pval_MOT4Run2_SH','rho_MOT2Run2_SH','pval_MOT2Run2_SH','rho_RestingStateRun2_SH','pval_RestingStateRun2_SH','-v7.3');

end

function [mean_area_MOT4_voxel,mean_area_MOT2_voxel,mean_area_RestingState_voxel] = getMeanParcellation(MOT4,MOT2,RestingState)

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
%load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

%%% GET ROI MEAN TIME SERIES

nTR = size(MOT4,4);

mean_area_MOT4_voxel = zeros(nNodes,nTR);
mean_area_MOT2_voxel = zeros(nNodes,nTR);
mean_area_RestingState_voxel = zeros(nNodes,nTR);

%ROI_Zeros = zeros(nNodes,4);

for iNode=1:nNodes

    idx_ROI = AAL_ROI(iNode).ID;
    idx_ROI_Label = AAL_ROI(iNode).Nom_L;
    idx_ROI_Label = strrep(idx_ROI_Label,'_','-');
    
    idx_voxels = find(AAL_img == idx_ROI);
    
    nVoxels = length(idx_voxels);
    
    area_MOT4 = zeros(nVoxels,nTR);
    area_MOT2 = zeros(nVoxels,nTR);
    area_RestingState = zeros(nVoxels,nTR);
    
    iz4 = 0;
    iz2 = 0;
    izr = 0;
    
    zeros_MOT4 = [];
    zeros_MOT2 = [];
    zeros_RestingState = [];
    
    for iVoxel=1:nVoxels
       
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
            
        area_MOT4(iVoxel,:) = MOT4(idxx,idxy,idxz,:);
        
        area_MOT2(iVoxel,:) = MOT2(idxx,idxy,idxz,:);
            
        area_RestingState(iVoxel,:) = RestingState(idxx,idxy,idxz,:);
        
        if sum(area_MOT4(iVoxel,:)) == 0
            
            iz4 = iz4 + 1;
            zeros_MOT4(iz4) = iVoxel;
            
        end
        
        if sum(area_MOT2(iVoxel,:)) == 0
            
            iz2 = iz2 + 1;
            zeros_MOT2(iz2) = iVoxel;
            
        end
        
        if sum(area_RestingState(iVoxel,:)) == 0
            
            izr = izr + 1;
            zeros_RestingState(izr) = iVoxel;
            
        end
        

    end
    
    area_MOT4(zeros_MOT4,:) = [];
    area_MOT2(zeros_MOT2,:) = [];
    area_RestingState(zeros_RestingState,:) = [];
    
    mn_MOT4 = mean(area_MOT4,1);
    mn_MOT2 = mean(area_MOT2,1);
    mn_RestingState = mean(area_RestingState,1);
      
    mean_area_MOT4_voxel(iNode,:) = mn_MOT4(:);
    mean_area_MOT2_voxel(iNode,:) = mn_MOT2(:);
    mean_area_RestingState_voxel(iNode,:) = mn_RestingState(:);
    
end

mean_area_MOT4_voxel = mean_area_MOT4_voxel';
mean_area_MOT2_voxel = mean_area_MOT2_voxel';
mean_area_RestingState_voxel = mean_area_RestingState_voxel';

end

function [mean_area_MOT4_voxel,mean_area_MOT2_voxel,mean_area_RestingState_voxel] = getMeanParcellationFPO(MOT4,MOT2,RestingState)

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26 27 28];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [17 18 19 20 57 58 59 60 61 62 63 64 65 66 67 68 69 70];

aal_idx = [idx_frontal,idx_occipital,idx_parietal];

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
%load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

%nNodes = 90;
nNodes = length(aal_idx);

%%% GET ROI MEAN TIME SERIES

nTR = size(MOT4,4);

mean_area_MOT4_voxel = zeros(nNodes,nTR);
mean_area_MOT2_voxel = zeros(nNodes,nTR);
mean_area_RestingState_voxel = zeros(nNodes,nTR);

%ROI_Zeros = zeros(nNodes,4);

for iNode=1:nNodes

    idx_ROI = AAL_ROI(aal_idx(iNode)).ID;
    idx_ROI_Label = AAL_ROI(aal_idx(iNode)).Nom_L;
    idx_ROI_Label = strrep(idx_ROI_Label,'_','-');
    
    idx_voxels = find(AAL_img == idx_ROI);
    
    nVoxels = length(idx_voxels);
    
    area_MOT4 = zeros(nVoxels,nTR);
    area_MOT2 = zeros(nVoxels,nTR);
    area_RestingState = zeros(nVoxels,nTR);
    
    iz4 = 0;
    iz2 = 0;
    izr = 0;
    
    zeros_MOT4 = [];
    zeros_MOT2 = [];
    zeros_RestingState = [];
    
    for iVoxel=1:nVoxels
       
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
            
        area_MOT4(iVoxel,:) = MOT4(idxx,idxy,idxz,:);
        
        area_MOT2(iVoxel,:) = MOT2(idxx,idxy,idxz,:);
            
        area_RestingState(iVoxel,:) = RestingState(idxx,idxy,idxz,:);
        
        if sum(area_MOT4(iVoxel,:)) == 0
            
            iz4 = iz4 + 1;
            zeros_MOT4(iz4) = iVoxel;
            
        end
        
        if sum(area_MOT2(iVoxel,:)) == 0
            
            iz2 = iz2 + 1;
            zeros_MOT2(iz2) = iVoxel;
            
        end
        
        if sum(area_RestingState(iVoxel,:)) == 0
            
            izr = izr + 1;
            zeros_RestingState(izr) = iVoxel;
            
        end
        

    end
    
    area_MOT4(zeros_MOT4,:) = [];
    area_MOT2(zeros_MOT2,:) = [];
    area_RestingState(zeros_RestingState,:) = [];
    
    mn_MOT4 = mean(area_MOT4,1);
    mn_MOT2 = mean(area_MOT2,1);
    mn_RestingState = mean(area_RestingState,1);
      
    mean_area_MOT4_voxel(iNode,:) = mn_MOT4(:);
    mean_area_MOT2_voxel(iNode,:) = mn_MOT2(:);
    mean_area_RestingState_voxel(iNode,:) = mn_RestingState(:);
    
end

mean_area_MOT4_voxel = mean_area_MOT4_voxel';
mean_area_MOT2_voxel = mean_area_MOT2_voxel';
mean_area_RestingState_voxel = mean_area_RestingState_voxel';

end

function [rho_MOT4,pval_MOT4,rho_MOT2,pval_MOT2,rho_RestingState,pval_RestingState] = getCorrelations(mean_area_MOT4_voxel,mean_area_MOT2_voxel,mean_area_RestingState_voxel)

[rho_MOT4,pval_MOT4] = corr(mean_area_MOT4_voxel);
[rho_MOT2,pval_MOT2] = corr(mean_area_MOT2_voxel);
[rho_RestingState,pval_RestingState] = corr(mean_area_RestingState_voxel);

end

function [rho_MOT4,pval_MOT4,rho_MOT2,pval_MOT2,rho_RestingState,pval_RestingState] = getCorrelations_corrcoef(mean_area_MOT4_voxel,mean_area_MOT2_voxel,mean_area_RestingState_voxel)

[rho_MOT4,pval_MOT4] = corrcoef(mean_area_MOT4_voxel);
[rho_MOT2,pval_MOT2] = corrcoef(mean_area_MOT2_voxel);
[rho_RestingState,pval_RestingState] = corrcoef(mean_area_RestingState_voxel);

end

function [rho_MOT4_FH,pval_MOT4_FH,rho_MOT2_FH,pval_MOT2_FH,rho_RestingState_FH,pval_RestingState_FH, rho_MOT4_SH,pval_MOT4_SH,rho_MOT2_SH,pval_MOT2_SH,rho_RestingState_SH,pval_RestingState_SH] = getCorrelationsHalf(mean_area_MOT4_voxel,mean_area_MOT2_voxel,mean_area_RestingState_voxel)

nTR = 300;

mean_area_MOT4_voxel_FH = mean_area_MOT4_voxel(1:(nTR/2),:);
mean_area_MOT4_voxel_SH = mean_area_MOT4_voxel(((nTR/2)+1):end,:);

mean_area_MOT2_voxel_FH = mean_area_MOT2_voxel(1:(nTR/2),:);
mean_area_MOT2_voxel_SH = mean_area_MOT2_voxel(((nTR/2)+1):end,:);

mean_area_RestingState_voxel_FH = mean_area_RestingState_voxel(1:(nTR/2),:);
mean_area_RestingState_voxel_SH = mean_area_RestingState_voxel(((nTR/2)+1):end,:);

[rho_MOT4_FH,pval_MOT4_FH] = corr(mean_area_MOT4_voxel_FH);
[rho_MOT2_FH,pval_MOT2_FH] = corr(mean_area_MOT2_voxel_FH);
[rho_RestingState_FH,pval_RestingState_FH] = corr(mean_area_RestingState_voxel_FH);

[rho_MOT4_SH,pval_MOT4_SH] = corr(mean_area_MOT4_voxel_SH);
[rho_MOT2_SH,pval_MOT2_SH] = corr(mean_area_MOT2_voxel_SH);
[rho_RestingState_SH,pval_RestingState_SH] = corr(mean_area_RestingState_voxel_SH);

end

function [rho_MOT4_FH,pval_MOT4_FH,rho_MOT2_FH,pval_MOT2_FH,rho_RestingState_FH,pval_RestingState_FH, rho_MOT4_SH,pval_MOT4_SH,rho_MOT2_SH,pval_MOT2_SH,rho_RestingState_SH,pval_RestingState_SH] = getCorrelationsHalfFPO(mean_area_MOT4_voxel,mean_area_MOT2_voxel,mean_area_RestingState_voxel)

% idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26 27 28];
% idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
% idx_parietal = [17 18 19 20 57 58 59 60 61 62 63 64 65 66 67 68 69 70];
% 
% aal_idx = [idx_frontal,idx_occipital,idx_parietal];

nTR = 300;

mean_area_MOT4_voxel_FH = mean_area_MOT4_voxel(1:(nTR/2),:);
mean_area_MOT4_voxel_SH = mean_area_MOT4_voxel(((nTR/2)+1):end,:);

mean_area_MOT2_voxel_FH = mean_area_MOT2_voxel(1:(nTR/2),:);
mean_area_MOT2_voxel_SH = mean_area_MOT2_voxel(((nTR/2)+1):end,:);

mean_area_RestingState_voxel_FH = mean_area_RestingState_voxel(1:(nTR/2),:);
mean_area_RestingState_voxel_SH = mean_area_RestingState_voxel(((nTR/2)+1):end,:);
 
% new_mean_area_MOT4_voxel_FH = zeros(nTR/2,length(aal_idx));
% for iidx=1:length(aal_idx)
%     new_mean_area_MOT4_voxel_FH(:,iidx) = mean_area_MOT4_voxel_FH(:,aal_idx(iidx)); 
% end
% 
% new_mean_area_MOT2_voxel_FH = zeros(nTR/2,length(aal_idx));
% for iidx=1:length(aal_idx)
%     new_mean_area_MOT2_voxel_FH(:,iidx) = mean_area_MOT2_voxel_FH(:,aal_idx(iidx)); 
% end
% 
% new_mean_area_RestingState_voxel_FH = zeros(nTR/2,length(aal_idx));
% for iidx=1:length(aal_idx)
%     new_mean_area_RestingState_voxel_FH(:,iidx) = mean_area_RestingState_voxel_FH(:,aal_idx(iidx)); 
% end
% 
% 
% new_mean_area_MOT4_voxel_SH = zeros(nTR/2,length(aal_idx));
% for iidx=1:length(aal_idx)
%     new_mean_area_MOT4_voxel_SH(:,iidx) = mean_area_MOT4_voxel_SH(:,aal_idx(iidx)); 
% end
% 
% new_mean_area_MOT2_voxel_SH = zeros(nTR/2,length(aal_idx));
% for iidx=1:length(aal_idx)
%     new_mean_area_MOT2_voxel_SH(:,iidx) = mean_area_MOT2_voxel_SH(:,aal_idx(iidx)); 
% end
% 
% new_mean_area_RestingState_voxel_SH = zeros(nTR/2,length(aal_idx));
% for iidx=1:length(aal_idx)
%     new_mean_area_RestingState_voxel_SH(:,iidx) = mean_area_RestingState_voxel_SH(:,aal_idx(iidx)); 
% end

% [rho_MOT4_FH,pval_MOT4_FH] = corr(new_mean_area_MOT4_voxel_FH);
% [rho_MOT2_FH,pval_MOT2_FH] = corr(new_mean_area_MOT2_voxel_FH);
% [rho_RestingState_FH,pval_RestingState_FH] = corr(new_mean_area_RestingState_voxel_FH);
% 
% [rho_MOT4_SH,pval_MOT4_SH] = corr(new_mean_area_MOT4_voxel_SH);
% [rho_MOT2_SH,pval_MOT2_SH] = corr(new_mean_area_MOT2_voxel_SH);
% [rho_RestingState_SH,pval_RestingState_SH] = corr(new_mean_area_RestingState_voxel_SH);

[rho_MOT4_FH,pval_MOT4_FH] = corr(mean_area_MOT4_voxel_FH);
[rho_MOT2_FH,pval_MOT2_FH] = corr(mean_area_MOT2_voxel_FH);
[rho_RestingState_FH,pval_RestingState_FH] = corr(mean_area_RestingState_voxel_FH);

[rho_MOT4_SH,pval_MOT4_SH] = corr(mean_area_MOT4_voxel_SH);
[rho_MOT2_SH,pval_MOT2_SH] = corr(mean_area_MOT2_voxel_SH);
[rho_RestingState_SH,pval_RestingState_SH] = corr(mean_area_RestingState_voxel_SH);


end

function plotResults(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL','.mat'));

plotCorrelation(settings,rho_MOT4Run1,pval_MOT4Run1,'High-Run1');
plotCorrelation(settings,rho_MOT4Run2,pval_MOT4Run2,'High-Run2');

plotCorrelation(settings,rho_MOT2Run1,pval_MOT2Run1,'Low-Run1');
plotCorrelation(settings,rho_MOT2Run2,pval_MOT2Run2,'Low-Run2');

plotCorrelation(settings,rho_RestingStateRun1,pval_RestingStateRun1,'RestingState-Run1');
plotCorrelation(settings,rho_RestingStateRun2,pval_RestingStateRun2,'RestingState-Run2');

end

function plotResultsHalf(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL-Half','.mat'));

plotCorrelation(settings,rho_MOT4Run1_FH,pval_MOT4Run1_FH,'High-Run1-FH');
plotCorrelation(settings,rho_MOT4Run2_FH,pval_MOT4Run2_FH,'High-Run2-FH');

plotCorrelation(settings,rho_MOT2Run1_FH,pval_MOT2Run1_FH,'Low-Run1-FH');
plotCorrelation(settings,rho_MOT2Run2_FH,pval_MOT2Run2_FH,'Low-Run2-FH');

plotCorrelation(settings,rho_RestingStateRun1_FH,pval_RestingStateRun1_FH,'RestingState-Run1-FH');
plotCorrelation(settings,rho_RestingStateRun2_FH,pval_RestingStateRun2_FH,'RestingState-Run2-FH');

plotCorrelation(settings,rho_MOT4Run1_SH,pval_MOT4Run1_SH,'High-Run1-FH');
plotCorrelation(settings,rho_MOT4Run2_SH,pval_MOT4Run2_SH,'High-Run2-FH');

plotCorrelation(settings,rho_MOT2Run1_SH,pval_MOT2Run1_SH,'Low-Run1-FH');
plotCorrelation(settings,rho_MOT2Run2_SH,pval_MOT2Run2_SH,'Low-Run2-FH');

plotCorrelation(settings,rho_RestingStateRun1_SH,pval_RestingStateRun1_SH,'RestingState-Run1-FH');
plotCorrelation(settings,rho_RestingStateRun2_SH,pval_RestingStateRun2_SH,'RestingState-Run2-FH');

end

function plotFPO(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL','.mat'));

plotCorrelationFPO(settings,rho_MOT4Run1,pval_MOT4Run1,'High-Run1');
plotCorrelationFPO(settings,rho_MOT4Run2,pval_MOT4Run2,'High-Run2');

plotCorrelationFPO(settings,rho_MOT2Run1,pval_MOT2Run1,'Low-Run1');
plotCorrelationFPO(settings,rho_MOT2Run2,pval_MOT2Run2,'Low-Run2');

plotCorrelationFPO(settings,rho_RestingStateRun1,pval_RestingStateRun1,'RestingState-Run1');
plotCorrelationFPO(settings,rho_RestingStateRun2,pval_RestingStateRun2,'RestingState-Run2');

end

function plotFPOHalf(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL-Half','.mat'));

plotCorrelationFPO(settings,rho_MOT4Run1_FH,pval_MOT4Run1_FH,'High-Run1-FH');
plotCorrelationFPO(settings,rho_MOT4Run2_FH,pval_MOT4Run2_FH,'High-Run2-FH');

plotCorrelationFPO(settings,rho_MOT2Run1_FH,pval_MOT2Run1_FH,'Low-Run1-FH');
plotCorrelationFPO(settings,rho_MOT2Run2_FH,pval_MOT2Run2_FH,'Low-Run2-FH');

plotCorrelationFPO(settings,rho_RestingStateRun1_FH,pval_RestingStateRun1_FH,'RestingState-Run1-FH');
plotCorrelationFPO(settings,rho_RestingStateRun2_FH,pval_RestingStateRun2_FH,'RestingState-Run2-FH');

plotCorrelationFPO(settings,rho_MOT4Run1_SH,pval_MOT4Run1_SH,'High-Run1-SH');
plotCorrelationFPO(settings,rho_MOT4Run2_SH,pval_MOT4Run2_SH,'High-Run2-SH');

plotCorrelationFPO(settings,rho_MOT2Run1_SH,pval_MOT2Run1_SH,'Low-Run1-SH');
plotCorrelationFPO(settings,rho_MOT2Run2_SH,pval_MOT2Run2_SH,'Low-Run2-SH');

plotCorrelationFPO(settings,rho_RestingStateRun1_SH,pval_RestingStateRun1_SH,'RestingState-Run1-SH');
plotCorrelationFPO(settings,rho_RestingStateRun2_SH,pval_RestingStateRun2_SH,'RestingState-Run2-SH');

end

function plotCorrelation(settings,rho,pval,label)

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26 27 28];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [17 18 19 20 57 58 59 60 61 62 63 64 65 66 67 68 69 70];
idx_temporal = [21 22 29 30 39 40 55 56 79 80 81 82 83 84 85 86 87 88 89 90];

nNode = 90;

CLIM = [-1 1];

pcriterion = 0.001;

kk = find( pval > pcriterion );
rho(kk) = zeros(size(kk));

f = figure;

clrmp = colormap('jet');
clrmp(32,:) = [1 1 1];

imagesc(rho,CLIM);
caxis(CLIM);
title(label);
colormap(clrmp);

h = colorbar;
set(h,'YLim',CLIM, 'YTick',[-1 -0.5 0 0.5 1], 'PlotBoxAspectRatio', [1 20 1]);

hold on;

nNode = nNode + 1;

line_position = idx_frontal(1);
plot( 1 + line_position*[1 1], [0 nNode], 'k', 'LineWidth', 2 );
plot( [0 nNode], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = idx_frontal(16);
plot( 1 + line_position*[1 1], [0 nNode], 'k', 'LineWidth', 2 );
plot( [0 nNode], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = idx_frontal(17);
plot( 1 + line_position*[1 1], [0 nNode], 'k', 'LineWidth', 2 );
plot( [0 nNode], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = idx_frontal(end);
plot( 1 + line_position*[1 1], [0 nNode], 'k', 'LineWidth', 2 );
plot( [0 nNode], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = idx_parietal(1);
plot( 1 + line_position*[1 1], [0 nNode], 'k', 'LineWidth', 2 );
plot( [0 nNode], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = idx_parietal(4);
plot( 1 + line_position*[1 1], [0 nNode], 'k', 'LineWidth', 2 );
plot( [0 nNode], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = idx_parietal(5);
plot( 1 + line_position*[1 1], [0 nNode], 'k', 'LineWidth', 2 );
plot( [0 nNode], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = idx_parietal(end);
plot( 1 + line_position*[1 1], [0 nNode], 'k', 'LineWidth', 2 );
plot( [0 nNode], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = idx_occipital(1);
plot( 1 + line_position*[1 1], [0 nNode], 'k', 'LineWidth', 2 );
plot( [0 nNode], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = idx_occipital(end);
plot( 1 + line_position*[1 1], [0 nNode], 'k', 'LineWidth', 2 );
plot( [0 nNode], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

Tick = [round(idx_frontal(1)+(idx_frontal(16)-idx_frontal(1))/2) round(idx_parietal(1)+(idx_parietal(4)-idx_parietal(1))/2) round(idx_frontal(17)+(idx_frontal(end)-idx_frontal(17))/2) round(idx_occipital(1)+(idx_occipital(end)-idx_occipital(1))/2) round(idx_parietal(5)+(idx_parietal(end)-idx_parietal(5))/2)];
TickLabel = { 'F' ; 'P'; 'F'; 'O' ; 'P' };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL','-',label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL','-',label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL','-',label,'.pdf'));


end

function plotCorrelationFPO(settings,rho,pval,label)

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26 27 28];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [17 18 19 20 57 58 59 60 61 62 63 64 65 66 67 68 69 70];

F_n = length(idx_frontal);
O_n = length(idx_occipital);
P_n = length(idx_parietal);

aal_idx = [idx_frontal,idx_occipital,idx_parietal];

nNode = 90;

CLIM = [-1 1];

pcriterion = 0.001;

kk = find( pval > pcriterion );
rho(kk) = zeros(size(kk));

new_rho = zeros(size(aal_idx,2),size(aal_idx,2));

for iidx=1:length(aal_idx)
    
    new_rho(iidx,:) = rho(aal_idx(iidx),aal_idx);
    
end

f = figure;

clrmp = colormap('jet');
clrmp(32,:) = [1 1 1];

imagesc(new_rho,CLIM);
caxis(CLIM);
title(label);
colormap(clrmp);

h = colorbar;
set(h,'YLim',CLIM, 'YTick',[-1 -0.5 0 0.5 1], 'PlotBoxAspectRatio', [1 20 1]);

hold on;

nNode = nNode + 1;

line_position = F_n;
plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

line_position = F_n + P_n;
plot( 1 + line_position*[1 1], [0 F_n+P_n+O_n], 'k', 'LineWidth', 2 );
plot( [0 F_n+P_n+O_n], 1 + line_position*[1 1], 'k', 'LineWidth', 2 );

Tick = [F_n/2 (F_n+P_n/2) (F_n+P_n+O_n/2)];
TickLabel = { 'F' ; 'P' ; 'O' };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL-FPO','-',label,'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL-FPO','-',label,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL-FPO','-',label,'.pdf'));


end

function rhoBtwRuns(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-mean-AAL','.mat'));

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26 27 28];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [17 18 19 20 57 58 59 60 61 62 63 64 65 66 67 68 69 70];

F_n = length(idx_frontal);
O_n = length(idx_occipital);
P_n = length(idx_parietal);

aal_idx = [idx_frontal,idx_occipital,idx_parietal];

new_rho_MOT4Run1 = zeros(size(aal_idx,2),size(aal_idx,2));
new_rho_MOT4Run2 = zeros(size(aal_idx,2),size(aal_idx,2));

new_rho_MOT2Run1 = zeros(size(aal_idx,2),size(aal_idx,2));
new_rho_MOT2Run2 = zeros(size(aal_idx,2),size(aal_idx,2));

new_rho_RestingStateRun1 = zeros(size(aal_idx,2),size(aal_idx,2));
new_rho_RestingStateRun2 = zeros(size(aal_idx,2),size(aal_idx,2));

for iidx=1:length(aal_idx)
    
    new_rho_MOT4Run1(iidx,:) = rho_MOT4Run1(aal_idx(iidx),aal_idx);
    new_rho_MOT4Run2(iidx,:) = rho_MOT4Run2(aal_idx(iidx),aal_idx);
    
    new_rho_MOT2Run1(iidx,:) = rho_MOT2Run1(aal_idx(iidx),aal_idx);
    new_rho_MOT2Run2(iidx,:) = rho_MOT2Run2(aal_idx(iidx),aal_idx);
    
    new_rho_RestingStateRun1(iidx,:) = rho_RestingStateRun1(aal_idx(iidx),aal_idx);
    new_rho_RestingStateRun2(iidx,:) = rho_RestingStateRun2(aal_idx(iidx),aal_idx);
    
end


rho_MOT4 = corr([new_rho_MOT4Run1(:),new_rho_MOT4Run2(:)]);
rho_MOT2 = corr([new_rho_MOT2Run1(:),new_rho_MOT2Run2(:)]);
rho_RestingState = corr([new_rho_RestingStateRun1(:),new_rho_RestingStateRun2(:)]);

disp(strcat('MOT4 correlation btw. runs:',num2str(rho_MOT4)));
disp(strcat('MOT2 correlation btw. runs:',num2str(rho_MOT2)));
disp(strcat('RestingState correlation btw. runs:',num2str(rho_RestingState)));


end
