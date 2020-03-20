

%%% LOAD SETTINGS

pvalue = 0.05;

settings_jan_0805;


%%% DO THE MATH

% Load Data
disp('get Data');

% FSL
lowhigh_load_all_data_FSL;

% SPM
clear MOT4Run1 MOT4Run1 MOT2Run1 MOT2Run2 RestingStateRun1 RestingStateRun2
get_at_this_preprocessed_step = settings.folders.smooth.name;
prefix_for_the_preprocessed_step = settings.folders.smooth.prefix;

lowhigh_load_all_data;


% Common Mask
disp('get Common Mask');

masks(1).mask = mask_MOT4Run1;
masks(2).mask = mask_MOT4Run2;
masks(3).mask = mask_MOT2Run1;
masks(4).mask = mask_MOT2Run2;
masks(5).mask = mask_RestingStateRun1;
masks(6).mask =  mask_RestingStateRun2;

common_mask = get_common_mask(masks);

% Mean Between Runs
disp('Mean Between Runs');

MOT4 = ( MOT4Run1 + MOT4Run2 ) ./ 2;
MOT2 = ( MOT2Run1 + MOT2Run2 ) ./ 2;
RestingState = ( RestingStateRun1 + RestingStateRun2 ) ./ 2;

MOT4(common_mask == 0) = 0;
MOT2(common_mask == 0) = 0;
RestingState(common_mask == 0) = 0;

% Load AAL
disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_labels = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_labels.ROI;

nNodes = 90;
nStructures = length(AAL_ROI);
nStructures = nNodes;

% Reshape
disp('Reshape');

nTR = size(MOT4,4);

MOT4_vec = reshape(MOT4,[size(MOT4,1)*size(MOT4,2)*size(MOT4,3),nTR]);
MOT2_vec = reshape(MOT2,[size(MOT2,1)*size(MOT2,2)*size(MOT2,3),nTR]);
RestingState_vec = reshape(RestingState,[size(RestingState,1)*size(RestingState,2)*size(RestingState,3),nTR]);

MOT4_vec = MOT4_vec';
MOT2_vec = MOT2_vec';
RestingState_vec = RestingState_vec';

% MA
disp('MA');

% MOT4_MA = char.empty;
% MOT2_MA = char.empty;
% RestingState_MA = char.empty;
% 
% MOT4_MA{1,1} = 'STRUCTURES';
% MOT2_MA{1,1} = 'STRUCTURES';
% RestingState_MA{1,1} = 'STRUCTURES';

MOT4_MA_mean = char.empty;
MOT2_MA_mean = char.empty;
RestingState_MA_mean = char.empty;

MOT4_MA_mean{1,1} = 'STRUCTURES';
MOT2_MA_mean{1,1} = 'STRUCTURES';
RestingState_MA_mean{1,1} = 'STRUCTURES';

MOT4_MA_mean_mat = zeros(nStructures,nStructures);
MOT2_MA_mean_mat = zeros(nStructures,nStructures);
RestingState_MA_mean_mat = zeros(nStructures,nStructures);

MOT4_MA_pval_mat = zeros(nStructures,nStructures);
MOT2_MA_pval_mat = zeros(nStructures,nStructures);
RestingState_MA_pval_mat = zeros(nStructures,nStructures);

for iStructure=1:nStructures
    
    disp(strcat('iStructure:',int2str(iStructure)));
   
    idx_structure = AAL_ROI(1,iStructure).ID;
    label_structure = AAL_ROI(1,iStructure).Nom_L;
    
    idx_voxels_structure = find(AAL_img == idx_structure);
    n_voxels = length(idx_voxels_structure);
    
%     MOT4_MA{iStructure+1,1} = label_structure;
%     MOT2_MA{iStructure+1,1} = label_structure;
%     RestingState_MA{iStructure+1,1} = label_structure;
    
    MOT4_MA_mean{iStructure+1,1} = label_structure;
    MOT2_MA_mean{iStructure+1,1} = label_structure;
    RestingState_MA_mean{iStructure+1,1} = label_structure;
     
    for iiStructure=iStructure:nStructures
        
        disp(strcat('iiStructure:',int2str(iiStructure)));
        
        iidx_structure = AAL_ROI(1,iiStructure).ID;
        ilabel_structure = AAL_ROI(1,iiStructure).Nom_L;
        
        iidx_voxels_structure = find(AAL_img == iidx_structure);
        in_voxels = length(iidx_voxels_structure);
        
%         MOT4_MA{1,iiStructure+1} = ilabel_structure;
%         MOT2_MA{1,iiStructure+1} = ilabel_structure;
%         RestingState_MA{1,iiStructure+1} = ilabel_structure;
   
        MOT4_MA_mean{1,iiStructure+1} = ilabel_structure;
        MOT2_MA_mean{1,iiStructure+1} = ilabel_structure;
        RestingState_MA_mean{1,iiStructure+1} = ilabel_structure;
   
        MOT4_MA_voxel = 0;
        MOT2_MA_voxel = 0;
        RestingState_MA_voxel = 0;
        
        disp('MOT4 corr');
        
        tic;
        
        from_voxel = zeros(nTR,n_voxels);
        to_voxel = zeros(nTR,in_voxels);
        
        from_voxel = MOT4_vec(:,idx_voxels_structure);
        to_voxel = MOT4_vec(:,iidx_voxels_structure);
        
        from_voxel_mean = mean(from_voxel,2);
        to_voxel_mean = mean(to_voxel,2);
        
        [r_mean,p_mean] = corr(from_voxel_mean,to_voxel_mean);
        
        MOT4_MA = p_mean<pvalue;
        MOT4_pval = p_mean;
        
%         [rho_mat, pval_mat] = corr([from_voxel,to_voxel]);
%         
%         pval_between = pval_mat(1:n_voxels,(n_voxels+1):end);
%         
%         MOT4_MA_voxel = length(find(pval_between<pvalue)) / (n_voxels*in_voxels);
        
        toc
        
        disp('MOT2 corr');
        
        tic;
        
        from_voxel = zeros(nTR,n_voxels);
        to_voxel = zeros(nTR,in_voxels);
        
        from_voxel = MOT2_vec(:,idx_voxels_structure);
        to_voxel = MOT2_vec(:,iidx_voxels_structure);
        
        from_voxel_mean = mean(from_voxel,2);
        to_voxel_mean = mean(to_voxel,2);
        
        [r_mean,p_mean] = corr(from_voxel_mean,to_voxel_mean);
        
        MOT2_MA = p_mean<pvalue;
        MOT2_pval = p_mean;
        
%         [rho_mat, pval_mat] = corr([from_voxel,to_voxel]);
%         
%         pval_between = pval_mat(1:n_voxels,(n_voxels+1):end);
%         
%         MOT2_MA_voxel = length(find(pval_between<pvalue)) / (n_voxels*in_voxels);
        
        toc
        
        disp('RestingState corr');
        
        tic;
        
        from_voxel = zeros(nTR,n_voxels);
        to_voxel = zeros(nTR,in_voxels);

        from_voxel = RestingState_vec(:,idx_voxels_structure);
        to_voxel = RestingState_vec(:,iidx_voxels_structure);
        
        from_voxel_mean = mean(from_voxel,2);
        to_voxel_mean = mean(to_voxel,2);
        
        [r_mean,p_mean] = corr(from_voxel_mean,to_voxel_mean);
        
        RestingState_MA = p_mean<pvalue;
        RestingState_pval = p_mean;
        
%         [rho_mat, pval_mat] = corr([from_voxel,to_voxel]);
%         
%         pval_between = pval_mat(1:n_voxels,(n_voxels+1):end);
%         
%         RestingState_MA_voxel = length(find(pval_between<pvalue)) / (n_voxels*in_voxels);

        toc
%         
%         MOT4_MA{iStructure+1,iiStructure+1} = MOT4_MA_voxel;
%         MOT2_MA{iStructure+1,iiStructure+1} = MOT2_MA_voxel;
%         RestingState_MA{iStructure+1,iiStructure+1} = RestingState_MA_voxel;
        
        MOT4_MA_mean{iStructure+1,iiStructure+1} = MOT4_MA;
        MOT2_MA_mean{iStructure+1,iiStructure+1} = MOT2_MA;
        RestingState_MA_mean{iStructure+1,iiStructure+1} = RestingState_MA;
        
        MOT4_MA_mean_mat(iStructure,iiStructure) = MOT4_MA;
        MOT2_MA_mean_mat(iStructure,iiStructure) = MOT2_MA;
        RestingState_MA_mean_mat(iStructure,iiStructure) = RestingState_MA;
        
        MOT4_MA_pval_mat(iStructure,iiStructure) = MOT4_pval;
        MOT2_MA_pval_mat(iStructure,iiStructure) = MOT2_pval;
        RestingState_MA_pval_mat(iStructure,iiStructure) = RestingState_pval;
        
    end
    
end
    
    