
disp('Loading Data');

nRuns = 32;

all_Track_pos = [];
all_Track_neg = [];
for iRun=1:nRuns
     load(strcat('Global-Track-Run-',int2str(iRun),'.mat'));
     all_Track_pos = [all_Track_pos; Track.run_pos];
     all_Track_neg = [all_Track_neg; Track.run_neg];
     clear Track
end

all_Passive_pos = [];
all_Passive_neg = [];
for iRun=(nRuns+1):nRuns*2
     load(strcat('Global-Passive-Run-',int2str(iRun),'.mat'));
     all_Passive_pos = [all_Passive_pos; Passive.run_pos];
     all_Passive_neg = [all_Passive_neg; Passive.run_neg];
     clear Passive
end

all_Rest_pos = [];
all_Rest_neg = [];
for iRun=(nRuns*2+1):nRuns*3
     load(strcat('Global-Rest-Run-',int2str(iRun),'.mat'));
     all_Rest_pos = [all_Rest_pos; Rest.run_pos];
     all_Rest_neg = [all_Rest_neg; Rest.run_neg];
     clear Rest
end

nTotalVoxels = 160990;
nROI = 90;

voxels_densities_Track_pos = zeros(nTotalVoxels,nRuns);
voxels_densities_Track_neg = zeros(nTotalVoxels,nRuns);

voxels_densities_Passive_pos = zeros(nTotalVoxels,nRuns);
voxels_densities_Passive_neg = zeros(nTotalVoxels,nRuns);

voxels_densities_Rest_pos = zeros(nTotalVoxels,nRuns);
voxels_densities_Rest_neg = zeros(nTotalVoxels,nRuns);

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

iiVoxel = 0;
for iROI=1:nROI
   
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    idx_voxels = find(AAL_img==AAL_ROI(iROI).ID);
    
    nVoxels = length(idx_voxels);
    
    disp(strcat('nVoxes=',int2str(nVoxels)));
    
    for iVoxel=1:nVoxels
       
        iiVoxel = iiVoxel + 1;
        
        %%% TRACK - NEG
        
        rois = squeeze(all_Track_neg(:,2));
        local_voxels = squeeze(all_Track_neg(:,3));
        
        idx_rois = find(rois==iROI);
        idx_local_voxels = find(local_voxels==iVoxel);
        
        idx_this_voxel = idx_rois(ismember(idx_rois,idx_local_voxels));
        
        tmp = squeeze(all_Track_neg(idx_this_voxel,4));
        if length(tmp) < nRuns
            
            tmp(end+1:end+nRuns-length(tmp)) = NaN;
            
        end
        
        voxels_densities_Track_neg(iiVoxel,:) = tmp(:);
        
        %%% TRACK - POS
        
        rois = squeeze(all_Track_pos(:,2));
        local_voxels = squeeze(all_Track_pos(:,3));
        
        idx_rois = find(rois==iROI);
        idx_local_voxels = find(local_voxels==iVoxel);
        
        idx_this_voxel = idx_rois(ismember(idx_rois,idx_local_voxels));
        
        tmp = squeeze(all_Track_pos(idx_this_voxel,4));
        if length(tmp) < nRuns
            
            tmp(end+1:end+nRuns-length(tmp)) = NaN;
            
        end
        
        voxels_densities_Track_pos(iiVoxel,:) = tmp(:);
        
        %%% PASSIVE - NEG
        
        rois = squeeze(all_Passive_neg(:,2));
        local_voxels = squeeze(all_Passive_neg(:,3));
        
        idx_rois = find(rois==iROI);
        idx_local_voxels = find(local_voxels==iVoxel);
        
        idx_this_voxel = idx_rois(ismember(idx_rois,idx_local_voxels));
        
        tmp = squeeze(all_Passive_neg(idx_this_voxel,4));
        if length(tmp) < nRuns
            
            tmp(end+1:end+nRuns-length(tmp)) = NaN;
            
        end
        
        voxels_densities_Passive_neg(iiVoxel,:) = tmp(:);
        
        %%% PASSIVE - POS
        
        rois = squeeze(all_Passive_pos(:,2));
        local_voxels = squeeze(all_Passive_pos(:,3));
        
        idx_rois = find(rois==iROI);
        idx_local_voxels = find(local_voxels==iVoxel);
        
        idx_this_voxel = idx_rois(ismember(idx_rois,idx_local_voxels));
        
        tmp = squeeze(all_Passive_pos(idx_this_voxel,4));
        if length(tmp) < nRuns
            
            tmp(end+1:end+nRuns-length(tmp)) = NaN;
            
        end
        
        voxels_densities_Passive_pos(iiVoxel,:) = tmp(:);
        
        %%% REST - NEG
        
        rois = squeeze(all_Rest_neg(:,2));
        local_voxels = squeeze(all_Rest_neg(:,3));
        
        idx_rois = find(rois==iROI);
        idx_local_voxels = find(local_voxels==iVoxel);
        
        idx_this_voxel = idx_rois(ismember(idx_rois,idx_local_voxels));
        
        tmp = squeeze(all_Rest_neg(idx_this_voxel,4));
        if length(tmp) < nRuns
            
            tmp(end+1:end+nRuns-length(tmp)) = NaN;
            
        end
        
        voxels_densities_Rest_neg(iiVoxel,:) = tmp(:);
        
        %%% REST - POS
        
        rois = squeeze(all_Rest_pos(:,2));
        local_voxels = squeeze(all_Rest_pos(:,3));
        
        idx_rois = find(rois==iROI);
        idx_local_voxels = find(local_voxels==iVoxel);
        
        idx_this_voxel = idx_rois(ismember(idx_rois,idx_local_voxels));
        
        tmp = squeeze(all_Rest_pos(idx_this_voxel,4));
        if length(tmp) < nRuns
            
            tmp(end+1:end+nRuns-length(tmp)) = NaN;
            
        end
        
        voxels_densities_Rest_pos(iiVoxel,:) = tmp(:);
        
    end
    
    
end