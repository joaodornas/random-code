function lowhigh_functional_parcellations

analysis_label = 'Functional-Parcels';

settings_jan_0805;
% settings_elena_2905;

all_networks_labels = {'DAN','VAN','VIS','LAN','DMN','FPC','SMN','AUD'};

% network_seeds = getFunctionalSeeds('DAN');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
% network_seeds = getFunctionalSeeds('VAN');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
% network_seeds = getFunctionalSeeds('VIS');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
% network_seeds = getFunctionalSeeds('LAN');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
% network_seeds = getFunctionalSeeds('DMN');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
% network_seeds = getFunctionalSeeds('AUD');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
% network_seeds = getFunctionalSeeds('FPC');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
% network_seeds = getFunctionalSeeds('SMN');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);

% network_label = 'DAN';
% generate3DVolumes(settings,analysis_label,network_label);
% network_label = 'VAN';
% generate3DVolumes(settings,analysis_label,network_label);

% network_label = 'DAN';
% network_seeds = getFunctionalSeeds('DAN');
% generateAllVolumes(settings,analysis_label,network_label,network_seeds);
% network_label = 'VAN';
% network_seeds = getFunctionalSeeds('VAN');
% generateAllVolumes(settings,analysis_label,network_label,network_seeds);

% network_label = 'DAN';
% computePopulationLevelMaps(settings,analysis_label,network_label);
% network_label = 'VAN';
% computePopulationLevelMaps(settings,analysis_label,network_label);
% network_label = 'VIS';
% computePopulationLevelMaps(settings,analysis_label,network_label);
% network_label = 'LAN';
% computePopulationLevelMaps(settings,analysis_label,network_label);
% network_label = 'DMN';
% computePopulationLevelMaps(settings,analysis_label,network_label);
% network_label = 'AUD';
% computePopulationLevelMaps(settings,analysis_label,network_label);
% network_label = 'SMN';
% computePopulationLevelMaps(settings,analysis_label,network_label);
% network_label = 'FPC';
% computePopulationLevelMaps(settings,analysis_label,network_label);

% network_label = 'DAN';
% network_seeds = getFunctionalSeeds('DAN');
% getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);
% network_label = 'VAN';
% network_seeds = getFunctionalSeeds('VAN');
% getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);
% network_label = 'VIS';
% network_seeds = getFunctionalSeeds('VIS');
% getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);
% network_label = 'LAN';
% network_seeds = getFunctionalSeeds('LAN');
% getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);
network_label = 'DMN';
network_seeds = getFunctionalSeeds('DMN');
getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);
network_label = 'FPC';
network_seeds = getFunctionalSeeds('FPC');
getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);
network_label = 'SMN';
network_seeds = getFunctionalSeeds('SMN');
getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);
network_label = 'AUD';
network_seeds = getFunctionalSeeds('AUD');
getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);


end

function compute_correlation_maps_per_condition(settings,analysis_label,network_seeds)

experiment_label = settings.codes.experiment;
subject_label = settings.codes.subject;
analysis_step_label = 'Corr-Map';

%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom_residual_voxel;

lowhigh_load_all_data_FSL;

disp('High-Run-1-Correlation-Map');
high_run1_correlation_map = getCorrelationMap(network_seeds,MOT4Run1,mask_MOT4Run1);

disp('High-Run-2-Correlation-Map');
high_run2_correlation_map = getCorrelationMap(network_seeds,MOT4Run2,mask_MOT4Run2);

disp('RestingState-Run-1-Correlation-Map');
rest_run1_correlation_map = getCorrelationMap(network_seeds,RestingStateRun1,mask_RestingStateRun1);

disp('RestingState-Run-2-Correlation-Map');
rest_run2_correlation_map = getCorrelationMap(network_seeds,RestingStateRun2,mask_RestingStateRun2);

save(strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_seeds.network_label,'.mat'),'high_run1_correlation_map','high_run2_correlation_map','rest_run1_correlation_map','rest_run2_correlation_map');

end

function mean_ROI = getROISeedMean(ROI,Run)
        
nTR = size(Run,4);

MNI_x_center = 45;
MNI_y_center = 63;
MNI_z_center = 36;

size_voxels_mm = 2;
ROI_radius_mm = 6 - size_voxels_mm;
ROI_radius_voxels = ROI_radius_mm/size_voxels_mm;

x_center = MNI_x_center + round(ROI.x/size_voxels_mm);
y_center = MNI_y_center + round(ROI.y/size_voxels_mm);
z_center = MNI_z_center + round(ROI.z/size_voxels_mm);

xgv = (x_center-ROI_radius_voxels):(x_center+ROI_radius_voxels);
ygv = (y_center-ROI_radius_voxels):(y_center+ROI_radius_voxels);
zgv = (z_center-ROI_radius_voxels):(z_center+ROI_radius_voxels);

[X,Y,Z] = meshgrid(xgv,ygv,zgv);

idx_voxels = find(X);
nVoxels = length(idx_voxels);

ROI_voxels = zeros(nVoxels,nTR);

for iVoxel=1:nVoxels
   
    ROI_voxels(iVoxel,:) = Run(X(iVoxel),Y(iVoxel),Z(iVoxel),:);
    
end

mean_ROI = mean(ROI_voxels,1)';

end

function correlation_map = getCorrelationMap(network_seeds,Run,mask)

nTR = size(Run,4);
nROIs = length(network_seeds.ROI);

idx_voxels = find(mask);
nVoxels = length(idx_voxels);
brain = zeros(nTR,nVoxels);

for iVoxel=1:nVoxels
    [idxx,idxy,idxz] = ind2sub(size(mask),idx_voxels(iVoxel));
    brain(1:end,iVoxel) = Run(idxx,idxy,idxz,:);
end

mean_ROI = zeros(nTR,nROIs);

for iROI=1:nROIs
    
    mean_ROI(1:nTR,iROI) = getROISeedMean(network_seeds.ROI(iROI),Run);
    
end

for iROI=1:nROIs
    
    disp(network_seeds.ROI(iROI).label);
   
    seed = mean_ROI(:,iROI);
    
    rho = zeros(nVoxels,1);
    rho_fisher = zeros(nVoxels,1);
    rho_zscore = zeros(nVoxels,1);
    rho_zscore_nTR = zeros(nVoxels,1);
    pval = zeros(nVoxels,1);
    
    for iVoxel=1:nVoxels
        
        voxel = brain(:,iVoxel);
        
        [rho(iVoxel),pval(iVoxel)] = corr(seed,voxel);
    
    end
    
    rho_fisher = 0.5 * log( (1 + rho)./(1 - rho) );
    
    acf = autocorr(seed,nTR-1);
    df = 1 / ( (1/nTR) + (2/nTR)*sum(acf.*acf) );
    
    rho_zscore_nTR = rho_fisher ./ sqrt( 1/(nTR-3) );
    
    rho_zscore = rho_fisher ./ sqrt(1/(df-3)) ;
   
    correlation_map.network_label = network_seeds.network_label;
    correlation_map.ROI(iROI).seed_label = network_seeds.ROI(iROI).label;
    correlation_map.ROI(iROI).seed_mean = seed;
    correlation_map.ROI(iROI).rho = rho;
    correlation_map.ROI(iROI).rho_fisher = rho_fisher;
    correlation_map.ROI(iROI).rho_zscore = rho_zscore;
    correlation_map.ROI(iROI).rho_zscore_nTR = rho_zscore_nTR;
    correlation_map.ROI(iROI).pval = pval;
    correlation_map.ROI(iROI).df = df;
    correlation_map.idx_voxels = idx_voxels;
    
end

end

function generate3DVolumes(settings,analysis_label,network_label)

experiment_label = settings.codes.experiment;
subject_label = settings.codes.subject;
analysis_step_label = 'Corr-Map';

load(strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'.mat'));

load_aal = nifti('ROI_MNI_V4.nii');

high_run1_3DVolume = get3DVolume(high_run1_correlation_map);

high_run2_3DVolume = get3DVolume(high_run2_correlation_map);

rest_run1_3DVolume = get3DVolume(rest_run1_correlation_map);

rest_run2_3DVolume = get3DVolume(rest_run2_correlation_map);

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = load_aal.dat.dim;

descrip = network_label;

fname = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-High-Run-1','.nii');
input_data = high_run1_3DVolume; 
lowhigh_save_image;

fname = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-High-Run-2','.nii');
input_data = high_run2_3DVolume; 
lowhigh_save_image;

fname = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-Rest-Run-1','.nii');
input_data = rest_run1_3DVolume; 
lowhigh_save_image;

fname = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-Rest-Run-2','.nii');
input_data = rest_run2_3DVolume; 
lowhigh_save_image;
    

end

function ThreeDVolume = get3DVolume(correlation_map)

idx_voxels = correlation_map.idx_voxels;

MNI_dim = [91 109 91];

nROI = length(correlation_map.ROI);

ThreeDVolume = zeros(MNI_dim);

pvalue = 0.01;
threshold = abs(norminv(pvalue/2,0,1));

for iROI=1:nROI
   
    idx_valid_voxels = find(correlation_map.ROI(iROI).rho_zscore > threshold | correlation_map.ROI(iROI).rho_zscore < -threshold);
    
    for iVoxel=1:length(idx_valid_voxels)
        
        [idxx,idxy,idxz] = ind2sub(MNI_dim,idx_voxels(idx_valid_voxels(iVoxel)));
        
        ThreeDVolume(idxx,idxy,idxz) = ThreeDVolume(idxx,idxy,idxz) + iROI;
        
    end
    
end

end

function ROI = getAllVolumes(correlation_map)

idx_voxels = correlation_map.idx_voxels;

MNI_dim = [91 109 91];

nROI = length(correlation_map.ROI);

for iROI=1:nROI
    
    Rho = zeros(MNI_dim);
    Rho_fisher = zeros(MNI_dim);
    Rho_zscore = zeros(MNI_dim);
    Rho_zscore_nTR = zeros(MNI_dim);
   
    for iVoxel=1:length(idx_voxels)
        
        [idxx,idxy,idxz] = ind2sub(MNI_dim,idx_voxels(iVoxel));
        
        Rho(idxx,idxy,idxz) = correlation_map.ROI(iROI).rho(iVoxel);
        Rho_fisher(idxx,idxy,idxz) = correlation_map.ROI(iROI).rho_fisher(iVoxel);
        Rho_zscore(idxx,idxy,idxz) = correlation_map.ROI(iROI).rho_zscore(iVoxel);
        Rho_zscore_nTR(idxx,idxy,idxz) = correlation_map.ROI(iROI).rho_zscore_nTR(iVoxel);
        
    end
    
    ROI(iROI).Rho = Rho;
    ROI(iROI).Rho_fisher = Rho_fisher;
    ROI(iROI).Rho_zscore = Rho_zscore;
    ROI(iROI).Rho_zscore_nTR = Rho_zscore_nTR;
    
end

end

function generateAllVolumes(settings,analysis_label,network_label,network_seeds)

experiment_label = settings.codes.experiment;
subject_label = settings.codes.subject;
analysis_step_label = 'Corr-Map';

load(strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'.mat'));

ConditionRun = 'High-Run-1';
ROI = getAllVolumes(high_run1_correlation_map);
saveAllVolumes(settings,ROI,ConditionRun,analysis_label,network_label,network_seeds);

ConditionRun = 'High-Run-2';
ROI = getAllVolumes(high_run2_correlation_map);
saveAllVolumes(settings,ROI,ConditionRun,analysis_label,network_label,network_seeds);

ConditionRun = 'Rest-Run-1';
ROI = getAllVolumes(rest_run1_correlation_map);
saveAllVolumes(settings,ROI,ConditionRun,analysis_label,network_label,network_seeds);

ConditionRun = 'Rest-Run-2';
ROI = getAllVolumes(rest_run2_correlation_map);
saveAllVolumes(settings,ROI,ConditionRun,analysis_label,network_label,network_seeds);

end

function saveAllVolumes(settings,ROI,ConditionRun,analysis_label,network_label,network_seeds)

experiment_label = settings.codes.experiment;
subject_label = settings.codes.subject;
analysis_step_label = 'Corr-Map';

nROI = length(ROI);
load_aal = nifti('ROI_MNI_V4.nii');

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = load_aal.dat.dim;

descrip = network_label;

for iROI=1:nROI

    fname = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-',ConditionRun,'-','Rho','-',network_seeds.ROI(iROI).label,'.nii');
    input_data = ROI(iROI).Rho; 
    lowhigh_save_image;

    fname = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-',ConditionRun,'-','Rho-Fisher','-',network_seeds.ROI(iROI).label,'.nii');
    input_data = ROI(iROI).Rho_fisher; 
    lowhigh_save_image;

    fname = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-',ConditionRun,'-','Rho-Zscore','-',network_seeds.ROI(iROI).label,'.nii');
    input_data = ROI(iROI).Rho_zscore; 
    lowhigh_save_image;

    fname = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-',ConditionRun,'-','Rho-Zscore-nTR','-',network_seeds.ROI(iROI).label,'.nii');
    input_data = ROI(iROI).Rho_zscore_nTR; 
    lowhigh_save_image;
    
end

end

function populationMap = computePopulationLevelMaps(settings,analysis_label,network_label)

experiment_label = settings.codes.experiment;
subject_label = settings.codes.subject;
analysis_step_label = 'Corr-Map';

load(strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'.mat'));

%%% RANDOM EFFECTS ANALYSIS

nROI = length(high_run1_correlation_map.ROI);
nRuns = 4;

idx_voxels_high_run1 = high_run1_correlation_map.idx_voxels;
idx_voxels_high_run2 = high_run2_correlation_map.idx_voxels;
idx_voxels_rest_run1 = rest_run1_correlation_map.idx_voxels;
idx_voxels_rest_run2 = rest_run2_correlation_map.idx_voxels;

all_voxels_and_runs = [idx_voxels_high_run1(:);idx_voxels_high_run2(:);idx_voxels_rest_run1(:);idx_voxels_rest_run2(:)];

[a,b] = hist(all_voxels_and_runs(:),unique(all_voxels_and_runs(:)));
common_voxels = b(find(a==nRuns));
nCommonVoxels = length(common_voxels);

idx_common_high_run1 = zeros(1,nCommonVoxels);
idx_common_high_run2 = zeros(1,nCommonVoxels);
idx_common_rest_run1 = zeros(1,nCommonVoxels);
idx_common_rest_run2 = zeros(1,nCommonVoxels);

for iVoxel=1:nCommonVoxels
    
    idx_common_high_run1(iVoxel) = find(high_run1_correlation_map.idx_voxels == common_voxels(iVoxel));
    idx_common_high_run2(iVoxel) = find(high_run2_correlation_map.idx_voxels == common_voxels(iVoxel));
    idx_common_rest_run1(iVoxel) = find(rest_run1_correlation_map.idx_voxels == common_voxels(iVoxel));
    idx_common_rest_run2(iVoxel) = find(rest_run2_correlation_map.idx_voxels == common_voxels(iVoxel));
    
end

for iROI=1:nROI
    
    variance(iROI).high_run1 = nanvar(high_run1_correlation_map.ROI(iROI).rho_fisher(idx_common_high_run1));

    variance(iROI).high_run2 = nanvar(high_run2_correlation_map.ROI(iROI).rho_fisher(idx_common_high_run2));

    variance(iROI).rest_run1 = nanvar(rest_run1_correlation_map.ROI(iROI).rho_fisher(idx_common_rest_run1));

    variance(iROI).rest_run2 = nanvar(rest_run2_correlation_map.ROI(iROI).rho_fisher(idx_common_rest_run2));

end

sample_variance = zeros(nROI,nCommonVoxels);
sigma_teta = zeros(nROI,nCommonVoxels);

for iVoxel=1:nCommonVoxels
    
   for iROI=1:nROI
    
        high_run1_rho_fisher = high_run1_correlation_map.ROI(iROI).rho_fisher(idx_common_high_run1(iVoxel));
        high_run2_rho_fisher = high_run2_correlation_map.ROI(iROI).rho_fisher(idx_common_high_run2(iVoxel));
        rest_run1_rho_fisher = rest_run1_correlation_map.ROI(iROI).rho_fisher(idx_common_rest_run1(iVoxel));
        rest_run2_rho_fisher = rest_run2_correlation_map.ROI(iROI).rho_fisher(idx_common_rest_run2(iVoxel));
        
        sample_variance(iROI,iVoxel) = var([high_run1_rho_fisher,high_run2_rho_fisher,rest_run1_rho_fisher,rest_run2_rho_fisher]); 
        
        my_sigma_teta = sample_variance(iROI,iVoxel) - (variance(iROI).high_run1 + variance(iROI).high_run2 + variance(iROI).rest_run1 + variance(iROI).rest_run2)/nRuns;
        
        if my_sigma_teta < 0; my_sigma_teta = 0; end
    
        sigma_teta(iROI,iVoxel) = my_sigma_teta;
   end
    
end

all_teta_star = zeros(nROI,nCommonVoxels);

roi_t = struct('voxels_t',zeros(1,nCommonVoxels),'voxels_z',zeros(1,nCommonVoxels));

for iROI=1:nROI
    
    for iVoxel=1:nCommonVoxels
        
        high_run1_rho_fisher = high_run1_correlation_map.ROI(iROI).rho_fisher(idx_common_high_run1(iVoxel));
        high_run2_rho_fisher = high_run2_correlation_map.ROI(iROI).rho_fisher(idx_common_high_run2(iVoxel));
        rest_run1_rho_fisher = rest_run1_correlation_map.ROI(iROI).rho_fisher(idx_common_rest_run1(iVoxel));
        rest_run2_rho_fisher = rest_run2_correlation_map.ROI(iROI).rho_fisher(idx_common_rest_run2(iVoxel));
        
        w_high_run1 = 1 / ( variance(iROI).high_run1 + sigma_teta(iROI,iVoxel) );
        w_high_run2 = 1 / ( variance(iROI).high_run2 + sigma_teta(iROI,iVoxel) );
        w_rest_run1 = 1 / ( variance(iROI).rest_run1 + sigma_teta(iROI,iVoxel) );
        w_rest_run2 = 1 / ( variance(iROI).rest_run2 + sigma_teta(iROI,iVoxel) );
        
        teta_star = ( (high_run1_rho_fisher * w_high_run1) + (high_run2_rho_fisher * w_high_run2) + (rest_run1_rho_fisher * w_rest_run1) + (rest_run2_rho_fisher * w_rest_run2) ) / (w_high_run1 + w_high_run2 + w_rest_run1 + w_rest_run2); 
        
        all_teta_star(iROI,iVoxel) = teta_star;
        
        rho_t(iROI).voxels_t(iVoxel) = teta_star / sqrt( 1 / ( w_high_run1 + w_high_run2 + w_rest_run1 + w_rest_run2 ) );

    end

    my_std = nanstd(rho_t(iROI).voxels_t);
    my_mean = nanmean(rho_t(iROI).voxels_t);
    
    rho_t(iROI).voxels_z = ( rho_t(iROI).voxels_t - my_mean ) / my_std;
    
end

populationMap.network_label = network_label;
populationMap.nROI = nROI;
populationMap.nRuns = nRuns;
populationMap.idx_voxels_high_run1 = idx_voxels_high_run1;
populationMap.idx_voxels_high_run2 = idx_voxels_high_run2;
populationMap.idx_voxels_rest_run1 = idx_voxels_rest_run1;
populationMap.idx_voxels_rest_run2 = idx_voxels_rest_run2;
populationMap.common_voxels = common_voxels;
populationMap.idx_common_high_run1 = idx_common_high_run1;
populationMap.idx_common_high_run2 = idx_common_high_run2;
populationMap.idx_common_rest_run1 = idx_common_rest_run1;
populationMap.idx_common_rest_run2 = idx_common_rest_run2;
populationMap.variance = variance;
populationMap.sample_variance = sample_variance;
populationMap.sigma_teta = sigma_teta;
populationMap.all_teta_star = all_teta_star;
populationMap.rho_t = rho_t;

analysis_step_label = 'Pop-Map';

save(strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'.mat'),'populationMap');

end

function getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds)

experiment_label = settings.codes.experiment;
subject_label = settings.codes.subject;
analysis_step_label = 'Pop-Map';

load(strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'.mat'));

MNI_dim = [91 109 91];
nROI = populationMap.nROI;

load_aal = nifti('ROI_MNI_V4.nii');

common_voxels = populationMap.common_voxels;

for iROI=1:nROI
   
    disp(network_seeds.ROI(iROI).label);

    seed_brain = zeros(MNI_dim);
    
    for iVoxel=1:length(common_voxels)
        
        [idxx,idxy,idxz] = ind2sub(MNI_dim,common_voxels(iVoxel));
       
        seed_brain(idxx,idxy,idxz) = populationMap.rho_t(iROI).voxels_z(iVoxel);
        
    end

    nifti_file = load_aal;
    offset = load_aal.dat.offset;
    scl_slope = load_aal.dat.scl_slope;
    scl_inter = load_aal.dat.scl_inter;

    dtype = 'FLOAT32';
    offset = 0;

    dim = load_aal.dat.dim;

    descrip = network_label;

    fname = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-',network_seeds.ROI(iROI).label,'.nii');
    input_data = seed_brain; 
    lowhigh_save_image;
    
end

end
