
function check_overlaps_densities
     
%%% NO Threshold
nothreshfiles{1} = 'LHR-Correlation-Contrast-PR-Neg-z';
nothreshfiles{2} = 'LHR-Correlation-Contrast-PR-Pos-z';

nothreshfiles{3} = 'LHR-Correlation-Contrast-RP-Neg-z';
nothreshfiles{4} = 'LHR-Correlation-Contrast-RP-Pos-z';

nothreshfiles{5} = 'LHR-Correlation-Contrast-PT-Neg-z';
nothreshfiles{6} = 'LHR-Correlation-Contrast-PT-Pos-z';

nothreshfiles{7} = 'LHR-Correlation-Contrast-TP-Neg-z';
nothreshfiles{8} = 'LHR-Correlation-Contrast-TP-Pos-z';

nothreshfolder = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\densities\contrasts-z';

%%% YES Threshold
yesthreshfiles{1} = 'LHR-Correlation-Contrast-PR-Neg-z-clu';
yesthreshfiles{2} = 'LHR-Correlation-Contrast-PR-Pos-z-clu';

yesthreshfiles{3} = 'LHR-Correlation-Contrast-RP-Neg-z-clu';
yesthreshfiles{4} = 'LHR-Correlation-Contrast-RP-Pos-z-clu';

yesthreshfiles{5} = 'LHR-Correlation-Contrast-PT-Neg-z-clu';
yesthreshfiles{6} = 'LHR-Correlation-Contrast-PT-Pos-z-clu';

yesthreshfiles{7} = 'LHR-Correlation-Contrast-TP-Neg-z-clu';
yesthreshfiles{8} = 'LHR-Correlation-Contrast-TP-Pos-z-clu';

yesthreshfolder = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\densities\clusters-dlh=1.5';

% passive_increase_nothresh = getOverlaps(nothreshfiles{1},nothreshfiles{2},nothreshfolder);
% passive_decrease_nothresh = getOverlaps(nothreshfiles{3},nothreshfiles{4},nothreshfolder);
% track_increase_nothresh = getOverlaps(nothreshfiles{7},nothreshfiles{8},nothreshfolder);
% track_decrease_nothresh = getOverlaps(nothreshfiles{5},nothreshfiles{6},nothreshfolder);
% 
% passive_increase_yesthresh = getOverlaps(yesthreshfiles{1},yesthreshfiles{2},yesthreshfolder);
% passive_decrease_yesthresh = getOverlaps(yesthreshfiles{3},yesthreshfiles{4},yesthreshfolder);
% track_increase_yesthresh = getOverlaps(yesthreshfiles{7},yesthreshfiles{8},yesthreshfolder);
% track_decrease_yesthresh = getOverlaps(yesthreshfiles{5},yesthreshfiles{6},yesthreshfolder);
% 
% save(strcat('LHR-All-Subjects-Densities-Overlaps.mat'),'passive_increase_nothresh','passive_decrease_nothresh','track_increase_nothresh','track_decrease_nothresh','passive_increase_yesthresh','passive_decrease_yesthresh','track_increase_yesthresh','track_decrease_yesthresh');

%%% TRACK
overlaps_nupd_track = getOverlapsCombinations(nothreshfiles{7},nothreshfiles{6},nothreshfolder,'T-nupd');
overlaps_ndpu_track = getOverlapsCombinations(nothreshfiles{5},nothreshfiles{8},nothreshfolder,'T-ndpu');
overlaps_nupu_track = getOverlapsCombinations(nothreshfiles{7},nothreshfiles{8},nothreshfolder,'T-nupu');
overlaps_ndpd_track = getOverlapsCombinations(nothreshfiles{5},nothreshfiles{6},nothreshfolder,'T-ndpd');

%%% PASSIVE
overlaps_nupd_passive = getOverlapsCombinations(nothreshfiles{1},nothreshfiles{4},nothreshfolder,'P-nupd');
overlaps_ndpu_passive = getOverlapsCombinations(nothreshfiles{3},nothreshfiles{2},nothreshfolder,'P-ndpu');
overlaps_nupu_passive = getOverlapsCombinations(nothreshfiles{1},nothreshfiles{2},nothreshfolder,'P-nupu');
overlaps_ndpd_passive = getOverlapsCombinations(nothreshfiles{3},nothreshfiles{4},nothreshfolder,'P-ndpd');

save(strcat('LHR-All-Subjects-Densities-Overlaps-Combinations.mat'),'overlaps_nupd_track','overlaps_ndpu_track','overlaps_nupu_track','overlaps_ndpd_track','overlaps_nupd_passive','overlaps_ndpu_passive','overlaps_nupu_passive','overlaps_ndpd_passive');


end

function overlaps = getOverlaps(img_name1,img_name2,folder)

img_neg_file = nifti(strcat(img_name1,'.nii'));
img_neg_file.dat.fname = strcat(folder,'\',img_neg_file.dat.fname);

img_pos_file = nifti(strcat(img_name2,'.nii'));
img_pos_file.dat.fname = strcat(folder,'\',img_pos_file.dat.fname);

img_neg = img_neg_file.dat(:,:,:);
img_pos = img_pos_file.dat(:,:,:);

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_ROI = 1:90;

overlaps = cell(90+1,1+1);

for iROI=idx_ROI
    
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    nROIvoxels = 0;
    
    for iVoxel=1:length(idx_voxels)
        
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
        
        neg_value = img_neg(idxx,idxy,idxz);
        pos_value = img_pos(idxx,idxy,idxz);
        
        if (neg_value ~= 0) && (pos_value ~= 0)
            
            nROIvoxels = nROIvoxels + 1;
            
        end
        
    end
    
    overlaps{1,2} = '#voxels';
    overlaps{1+iROI,1} = area_label;
    overlaps{1+iROI,2} = nROIvoxels;
    
end

end

function overlaps = getOverlapsCombinations(img_name1,img_name2,folder,overlap_label)

img_neg_file = nifti(strcat(img_name1,'.nii'));
img_neg_file.dat.fname = strcat(folder,'\',img_neg_file.dat.fname);

img_pos_file = nifti(strcat(img_name2,'.nii'));
img_pos_file.dat.fname = strcat(folder,'\',img_pos_file.dat.fname);

img_neg = img_neg_file.dat(:,:,:);
img_pos = img_pos_file.dat(:,:,:);

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_ROI = 1:90;

overlaps = cell(90+1,1+1);

overlap_img = zeros(size(AAL_img));

for iROI=idx_ROI
    
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    nROIvoxels = 0;
    
    for iVoxel=1:length(idx_voxels)
        
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
        
        neg_value = img_neg(idxx,idxy,idxz);
        pos_value = img_pos(idxx,idxy,idxz);
        
%         if ( (neg_value ~= 0) && (pos_value ~= 0) ) || ( (neg_value ~= 0) && (pos_value == 0) ) || ( (neg_value == 0) && (pos_value ~= 0) )
            
%         if ( (neg_value ~= 0) && (pos_value == 0) )
% 
%             nROIvoxels = nROIvoxels + 1;
% 
%             overlap_img(idxx,idxy,idxz) = neg_value;
%             
%         end
% 
%         if ( (neg_value == 0) && (pos_value ~= 0) )
% 
%             nROIvoxels = nROIvoxels + 1;
% 
%             overlap_img(idxx,idxy,idxz) = pos_value;
%             
%         end

        if ( (neg_value ~= 0) && (pos_value ~= 0) )

            nROIvoxels = nROIvoxels + 1;

            overlap_img(idxx,idxy,idxz) = min([neg_value pos_value]);
            
        end
        
    end
    
    overlaps{1,2} = '#voxels';
    overlaps{1+iROI,1} = area_label;
    overlaps{1+iROI,2} = nROIvoxels;
    
end

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;
dtype = 'FLOAT32';
offset = 0;
dim = load_aal.dat.dim;

descrip = 'Overlaps';
experiment_label = 'LHR';

fname = strcat(experiment_label,'-',overlap_label,'.nii');
input_data = overlap_img;
real_save_image;

end
