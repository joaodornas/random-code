function lowhigh_kmeans_voxels_AAL

%% SETTINGS

settings_jan_0805;
%settings_elena_2905;

doTheMath(settings);

return


function doTheMath(settings)

%% LOAD DATA

lowhigh_load_all_data_FSL;

MOT4 = (MOT4Run1 + MOT4Run2)./2;
MOT2 = (MOT2Run1 + MOT2Run2)./2;
RestingState = (RestingStateRun1 + RestingStateRun2)./2;

%% AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

disp('Frontal');
idx_nodes = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];
doKMeans(MOT4,MOT2,RestingState,AAL_img,AAL_ROI,settings,idx_nodes,'Frontal');

disp('Occipital');
idx_nodes = [49 50 51 52 53 54];
doKMeans(MOT4,MOT2,RestingState,AAL_img,AAL_ROI,settings,idx_nodes,'Occipital');

disp('Parietal');
idx_nodes = [59 60 61 62];
doKMeans(MOT4,MOT2,RestingState,AAL_img,AAL_ROI,settings,idx_nodes,'Parietal');

disp('Temporal');
idx_nodes = [81 82 83 84 85 86 87 88 89 90];
doKMeans(MOT4,MOT2,RestingState,AAL_img,AAL_ROI,settings,idx_nodes,'Temporal');

return

function doKMeans(MOT4,MOT2,RestingState,AAL_img,AAL_ROI,settings,idx_nodes,kind)

%% KMEANS

k = 4;

output_MOT4 = zeros(size(AAL_img));
output_MOT2 = zeros(size(AAL_img));
output_RestingState = zeros(size(AAL_img));

iAAL = 0;

disp('process voxels');
for iNode=1:length(idx_nodes)
   
    idx = idx_nodes(iNode);
    
    idx_AAL = AAL_ROI(idx).ID;
    
    idx_voxels = find(AAL_img==idx_AAL);
    
    for iVoxel=1:length(idx_voxels)
        
        iAAL = iAAL + 1;
    
        idx = idx_voxels(iVoxel);
        
        [idxx(iAAL), idxy(iAAL), idxz(iAAL)] = ind2sub(size(AAL_img),idx);
        
        area_MOT4(iAAL,:) = MOT4(idxx(iAAL),idxy(iAAL),idxz(iAAL),:);
        area_MOT2(iAAL,:) = MOT2(idxx(iAAL),idxy(iAAL),idxz(iAAL),:);
        area_RestingState(iAAL,:) = RestingState(idxx(iAAL),idxy(iAAL),idxz(iAAL),:);
    
    end
    
end

disp('compute kmeans');
clusters_area_MOT4 = kmeans(area_MOT4,k);
clusters_area_MOT2 = kmeans(area_MOT2,k);
clusters_area_RestingState = kmeans(area_RestingState,k);

% ma = k;
% mi = 1;
% col = jet;
% C_MOT4 = squeeze(ind2rgb(floor(((clusters_area_MOT4(:)-mi)/(ma-mi))*size(col,1)),col));
% C_MOT2 = squeeze(ind2rgb(floor(((clusters_area_MOT2(:)-mi)/(ma-mi))*size(col,1)),col));
% C_RestingState = squeeze(ind2rgb(floor(((clusters_area_RestingState(:)-mi)/(ma-mi))*size(col,1)),col));

C_MOT4 = clusters_area_MOT4;
C_MOT2 = clusters_area_MOT2;
C_RestingState = clusters_area_RestingState;

for iAAL=1:length(idxx)
    
    output_MOT4(idxx(iAAL),idxy(iAAL),idxz(iAAL)) = C_MOT4(iAAL);
    output_MOT2(idxx(iAAL),idxy(iAAL),idxz(iAAL)) = C_MOT2(iAAL);
    output_RestingState(idxx(iAAL),idxy(iAAL),idxz(iAAL)) = C_RestingState(iAAL);

end

run = 1;

dim = [settings.mot4.FSL.run(run).dim(1), settings.mot4.FSL.run(run).dim(2),settings.mot4.FSL.run(run).dim(3)];

nifti_file = settings.mot4.FSL.run(run).nifti_file;
dtype = settings.mot4.FSL.run(run).dtype;
offset = settings.mot4.FSL.run(run).offset;
scl_slope = settings.mot4.FSL.run(run).scl_slope;
scl_inter = settings.mot4.FSL.run(run).scl_inter;
    
dtype = 'FLOAT32';
offset = 0;

descrip = 'k-means';
 
fname = strcat('Low-High-',settings.folders.subject,'-MOT4-',kind,'-k-4-means.nii');

input_data = output_MOT4; 

lowhigh_save_image;     

fname = strcat('Low-High-',settings.folders.subject,'-MOT2-',kind,'-k-4-means.nii');

input_data = output_MOT2; 

lowhigh_save_image;     

fname = strcat('Low-High-',settings.folders.subject,'-RestingState-',kind,'-k-4-means.nii');

input_data = output_RestingState; 

lowhigh_save_image;     

return





