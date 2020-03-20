
function real_functional_parcellations

% settings_subj1_2210;
% doTheMath(settings);
% % checkOverlaps(settings);
% % computeParcellationVolume(settings);
% % paintParcellationVolumeSlices(settings);
% clear settings;
% 
% settings_subj2_2610;
% doTheMath(settings);
% % checkOverlaps(settings);
% % computeParcellationVolume(settings);
% % paintParcellationVolumeSlices(settings);
% clear settings;
% 
% settings_subj3_0311;
% doTheMath(settings);
% % checkOverlaps(settings);
% % computeParcellationVolume(settings);
% % paintParcellationVolumeSlices(settings);
% clear settings;
% 
% settings_subj4_0211;
% doTheMath(settings);
% % checkOverlaps(settings);
% % computeParcellationVolume(settings);
% % paintParcellationVolumeSlices(settings);
% clear settings;
% 
% settings_subj5_0211;
% doTheMath(settings);
% % checkOverlaps(settings);
% % computeParcellationVolume(settings);
% % paintParcellationVolumeSlices(settings);
% clear settings;
% 
% settings_subj6_2411;
% doTheMath(settings);
% % checkOverlaps(settings);
% % computeParcellationVolume(settings);
% % paintParcellationVolumeSlices(settings);
% clear settings;
% 
% settings_subj7_1401;
% doTheMath(settings);
% % checkOverlaps(settings);
% % computeParcellationVolume(settings);
% % paintParcellationVolumeSlices(settings);
% clear settings;
% 
% settings_subj8_1401;
% doTheMath(settings);
% % checkOverlaps(settings);
% % computeParcellationVolume(settings);
% % paintParcellationVolumeSlices(settings);
% clear settings;

%all_settings = getAllSettings;
% all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

% doTheMathGroup;
% computeParcellationVolumeGroupLevelPerNet(all_networks_labels);
%computeParcellationVolumeGroupLevel;

% a_few_nets = {'LAN','DMN'};
% %doTheMathGroup;
% a_few_nets = {'LAN'};
% computeParcellationVolumeGroupLevelPerNet(a_few_nets);
% a_few_nets = {'DMN'};
% computeParcellationVolumeGroupLevelPerNet(a_few_nets);

% getMD758densityPerNET;

paintParcellationVolumeSlicesFromMySimpleModel;

end


function doTheMathSubject(settings)

analysis_label = 'Functional-Parcels';

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

%%% CORRELATION MAPS

% network_seeds = getFunctionalSeeds_v5('DAN');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
% network_seeds = getFunctionalSeeds_v5('VAN');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
% network_seeds = getFunctionalSeeds_v5('VIS');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
network_seeds = getFunctionalSeeds_v5('LAN');
compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
% network_seeds = getFunctionalSeeds_v5('DMN');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
% network_seeds = getFunctionalSeeds_v5('FPC');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
% network_seeds = getFunctionalSeeds_v5('SMN');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
% network_seeds = getFunctionalSeeds_v5('AUD');
% compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);

%%% POPULATION MAPS - SUBJECT LEVEL

% network_label = 'DAN';
% computePopulationLevelMaps(settings,analysis_label,network_label);
% network_label = 'VAN';
% computePopulationLevelMaps(settings,analysis_label,network_label);
% network_label = 'VIS';
% computePopulationLevelMaps(settings,analysis_label,network_label);
network_label = 'LAN';
computePopulationLevelMaps(settings,analysis_label,network_label);
% network_label = 'DMN';
% computePopulationLevelMaps(settings,analysis_label,network_label);
% network_label = 'FPC';
% computePopulationLevelMaps(settings,analysis_label,network_label);
% network_label = 'SMN';
% computePopulationLevelMaps(settings,analysis_label,network_label);
% network_label = 'AUD';
% computePopulationLevelMaps(settings,analysis_label,network_label);

%%% VOLUMES FROM POPULATION MAPS - SUBJECT LEVEL

% network_label = 'DAN';
% network_seeds = getFunctionalSeeds_v5('DAN');
% getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);
% network_label = 'VAN';
% network_seeds = getFunctionalSeeds_v5('VAN');
% getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);
% network_label = 'VIS';
% network_seeds = getFunctionalSeeds_v5('VIS');
% getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);
network_label = 'LAN';
network_seeds = getFunctionalSeeds_v5('LAN');
getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);
% network_label = 'DMN';
% network_seeds = getFunctionalSeeds_v5('DMN');
% getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);
% network_label = 'FPC';
% network_seeds = getFunctionalSeeds_v5('FPC');
% getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);
% network_label = 'SMN';
% network_seeds = getFunctionalSeeds_v5('SMN');
% getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);
% network_label = 'AUD';
% network_seeds = getFunctionalSeeds_v5('AUD');
% getVolumesFromPopulationLevelMaps(settings,analysis_label,network_label,network_seeds);

end

function doTheMathGroup

analysis_label = 'Functional-Parcels';

%all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

all_settings = getAllSettings;

%%% CORRELATION MAPS

for iSettings=1:length(all_settings)
    
    settings = all_settings(iSettings).settings;

%     network_seeds = getFunctionalSeeds_v5('DAN');
%     compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
%     network_seeds = getFunctionalSeeds_v5('VAN');
%     compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
%     network_seeds = getFunctionalSeeds_v5('VIS');
%     compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
    network_seeds = getFunctionalSeeds_v5('LAN');
    compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
    network_seeds = getFunctionalSeeds_v5('DMN');
    compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
%     network_seeds = getFunctionalSeeds_v5('FPC');
%     compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
%     network_seeds = getFunctionalSeeds_v5('SMN');
%     compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
    % network_seeds = getFunctionalSeeds_v5('AUD');
    % compute_correlation_maps_per_condition(settings,analysis_label,network_seeds);
    
end

%%% POPULATION MAPS - GROUP LEVEL

% network_label = 'DAN';
% computePopulationLevelMapsGroupLevel(all_settings,network_label);
% network_label = 'VAN';
% computePopulationLevelMapsGroupLevel(all_settings,network_label);
% network_label = 'VIS';
% computePopulationLevelMapsGroupLevel(all_settings,network_label);
network_label = 'LAN';
computePopulationLevelMapsGroupLevel(all_settings,network_label);
network_label = 'DMN';
computePopulationLevelMapsGroupLevel(all_settings,network_label);
% network_label = 'FPC';
% computePopulationLevelMapsGroupLevel(all_settings,network_label);
% network_label = 'SMN';
% computePopulationLevelMapsGroupLevel(all_settings,network_label);
% network_label = 'AUD';
% computePopulationLevelMapsGroupLevel(all_settings,network_label);

%%% VOLUMES FROM POPULATION MAPS - GROUP LEVEL

% network_label = 'DAN';
% network_seeds = getFunctionalSeeds_v5('DAN');
% getVolumesFromPopulationLevelMapsGroupLevel(analysis_label,network_label,network_seeds);
% network_label = 'VAN';
% network_seeds = getFunctionalSeeds_v5('VAN');
% getVolumesFromPopulationLevelMapsGroupLevel(analysis_label,network_label,network_seeds);
% network_label = 'VIS';
% network_seeds = getFunctionalSeeds_v5('VIS');
% getVolumesFromPopulationLevelMapsGroupLevel(analysis_label,network_label,network_seeds);
network_label = 'LAN';
network_seeds = getFunctionalSeeds_v5('LAN');
getVolumesFromPopulationLevelMapsGroupLevel(analysis_label,network_label,network_seeds);
network_label = 'DMN';
network_seeds = getFunctionalSeeds_v5('DMN');
getVolumesFromPopulationLevelMapsGroupLevel(analysis_label,network_label,network_seeds);
% network_label = 'FPC';
% network_seeds = getFunctionalSeeds_v5('FPC');
% getVolumesFromPopulationLevelMapsGroupLevel(analysis_label,network_label,network_seeds);
% network_label = 'SMN';
% network_seeds = getFunctionalSeeds_v5('SMN');
% getVolumesFromPopulationLevelMapsGroupLevel(analysis_label,network_label,network_seeds);
% network_label = 'AUD';
% network_seeds = getFunctionalSeeds_v5('AUD');
% getVolumesFromPopulationLevelMapsGroupLevel(analysis_label,network_label,network_seeds);

end

%%% MAIN FUNCTIONS

function compute_correlation_maps_per_condition(settings,analysis_label,network_seeds)

experiment_label = settings.codes.experiment;
subject_label = settings.codes.subject;
analysis_step_label = 'Corr-Map';

%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

kind = 'RestingState';
for irun=1:4
    [RestingState(irun).run, RestingState(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
end

for irun=1:4
    
    disp(strcat('RestingState-Run-',int2str(irun),'-Correlation-Map'));
    Correlation_Map(irun).maps = getCorrelationMap(network_seeds,RestingState(irun).run,RestingState(irun).mask);

end

save(strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_seeds.network_label,'.mat'),'Correlation_Map');

end

function mean_ROI = getROISeedMean(ROI,Run)
        
nTR = size(Run,4);

MNI_x_center = 45;
MNI_y_center = 63;
MNI_z_center = 36;

size_voxels_mm = 2;
ROI_radius_mm = 6 - size_voxels_mm;
ROI_radius_voxels = ROI_radius_mm/size_voxels_mm;

x_center = MNI_x_center + round(ROI.x/size_voxels_mm) + 1;
y_center = MNI_y_center + round(ROI.y/size_voxels_mm) + 1;
z_center = MNI_z_center + round(ROI.z/size_voxels_mm) + 1;

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

%%% SUBJECT LEVEL

function populationMap = computePopulationLevelMaps(settings,analysis_label,network_label)

experiment_label = settings.codes.experiment;
subject_label = settings.codes.subject;
analysis_step_label = 'Corr-Map';

load(strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'.mat'));

%%% RANDOM EFFECTS ANALYSIS

nROI = length(Correlation_Map(1).maps.ROI);
nRuns = length(Correlation_Map);

all_voxels_and_runs = [];
for irun=1:nRuns
    idx_voxels(irun).idx = Correlation_Map(irun).maps.idx_voxels;
    all_voxels_and_runs = [all_voxels_and_runs;idx_voxels(irun).idx(:)];
end

[a,b] = hist(all_voxels_and_runs(:),unique(all_voxels_and_runs(:)));
common_voxels = b(find(a==nRuns));
nCommonVoxels = length(common_voxels);

for irun=1:nRuns
    idx_common(irun).common_voxels = zeros(1,nCommonVoxels);
end

for iVoxel=1:nCommonVoxels
    
    for irun=1:nRuns
        idx_common(irun).common_voxels(iVoxel) = find(Correlation_Map(irun).maps.idx_voxels == common_voxels(iVoxel));
    end
    
end

for iROI=1:nROI
    
    for irun=1:nRuns
        
        variance(iROI).run(irun).var = nanvar(Correlation_Map(irun).maps.ROI(iROI).rho_fisher(idx_common(irun).common_voxels));

    end
    
end

sample_variance = zeros(nROI,nCommonVoxels);
sigma_teta = zeros(nROI,nCommonVoxels);

for iVoxel=1:nCommonVoxels
    
   for iROI=1:nROI
    
       for irun=1:nRuns
           
            run(irun).rho_fisher = Correlation_Map(irun).maps.ROI(iROI).rho_fisher(idx_common(irun).common_voxels(iVoxel));
        
       end
       
       sample_variance(iROI,iVoxel) = var([run(1).rho_fisher,run(2).rho_fisher,run(3).rho_fisher,run(4).rho_fisher]); 
        
       my_sigma_teta = sample_variance(iROI,iVoxel) - (variance(iROI).run(1).var + variance(iROI).run(2).var + variance(iROI).run(3).var + variance(iROI).run(4).var)/nRuns;
        
       if my_sigma_teta < 0; my_sigma_teta = 0; end
    
       sigma_teta(iROI,iVoxel) = my_sigma_teta;
   
   end
    
end

all_teta_star = zeros(nROI,nCommonVoxels);

rho_t = struct('voxels_t',zeros(1,nCommonVoxels),'voxels_z',zeros(1,nCommonVoxels));

for iROI=1:nROI
    
    for iVoxel=1:nCommonVoxels
        
        run1_rho_fisher = Correlation_Map(1).maps.ROI(iROI).rho_fisher(idx_common(1).common_voxels(iVoxel));
        run2_rho_fisher = Correlation_Map(2).maps.ROI(iROI).rho_fisher(idx_common(2).common_voxels(iVoxel));
        run3_rho_fisher = Correlation_Map(3).maps.ROI(iROI).rho_fisher(idx_common(3).common_voxels(iVoxel));
        run4_rho_fisher = Correlation_Map(4).maps.ROI(iROI).rho_fisher(idx_common(4).common_voxels(iVoxel));
        
        w_run1 = 1 / ( variance(iROI).run(1).var + sigma_teta(iROI,iVoxel) );
        w_run2 = 1 / ( variance(iROI).run(2).var + sigma_teta(iROI,iVoxel) );
        w_run3 = 1 / ( variance(iROI).run(3).var + sigma_teta(iROI,iVoxel) );
        w_run4 = 1 / ( variance(iROI).run(4).var + sigma_teta(iROI,iVoxel) );
        
        teta_star = ( (run1_rho_fisher * w_run1) + (run2_rho_fisher * w_run2) + (run3_rho_fisher * w_run3) + (run4_rho_fisher * w_run4) ) / (w_run1 + w_run2 + w_run3 + w_run4); 
        
        all_teta_star(iROI,iVoxel) = teta_star;
        
        rho_t(iROI).voxels_t(iVoxel) = teta_star / sqrt( 1 / (w_run1 + w_run2 + w_run3 + w_run4) );

    end

    my_std = nanstd(rho_t(iROI).voxels_t);
    my_mean = nanmean(rho_t(iROI).voxels_t);
    
    rho_t(iROI).voxels_z = ( rho_t(iROI).voxels_t - my_mean ) / my_std;
    
end

populationMap.network_label = network_label;
populationMap.nROI = nROI;
populationMap.nRuns = nRuns;
populationMap.idx_voxels = idx_voxels(1:irun);
populationMap.common_voxels = common_voxels;
populationMap.idx_common = idx_common;
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
    real_save_image;
    
end

end

function checkOverlaps(settings)

experiment_label = settings.codes.experiment;
subject_label = settings.codes.subject;
analysis_label = 'Functional-Parcels';
analysis_step_label = 'Pop-Map';

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','DMN','AUD'};

nNet = length(all_networks_labels);

for iNet=1:nNet

    network_label = all_networks_labels{iNet};
    
    network_seeds = getFunctionalSeeds_v5(all_networks_labels{iNet});
    
    nROI = length(network_seeds.ROI);
    
    for iROI=1:nROI
        
        seed_volume_filename = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-',network_seeds.ROI(iROI).label,'-clu');
        seed_volume = nifti(strcat(seed_volume_filename,'.nii'));
        seed_volume.dat.fname = strcat(settings.folders.main,'/',settings.folders.experiment,'/',settings.folders.subject,'/','output','/','Functional_Parcellation','/','Seeds_Volumes_Clusters','/',seed_volume_filename,'.nii');
        
        net(iNet).ROI(iROI).idx_voxels = find(seed_volume.dat(:,:,:));
        if size(net(iNet).ROI(iROI).idx_voxels,1) ~= 1, net(iNet).ROI(iROI).idx_voxels = net(iNet).ROI(iROI).idx_voxels'; end
    
    end
    
end

all_voxels = [];
nAllROI = 0;
for iNet=1:nNet
    
    nROI = length(net(iNet).ROI);
    
    for iROI=1:nROI
        
        nAllROI = nAllROI + 1;
        
        all_voxels = [all_voxels,net(iNet).ROI(iROI).idx_voxels];
        
    end
    
end

[a,b] = hist(all_voxels(:),unique(all_voxels(:)));
common_voxels = b(find(a>1));
nCommonVoxels = length(common_voxels);

save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-','overlaps','.mat'),'common_voxels','nCommonVoxels');

end

function computeParcellationVolume(settings)

experiment_label = settings.codes.experiment;
subject_label = settings.codes.subject;
analysis_label = 'Functional-Parcels';
analysis_step_label = 'Pop-Map';

%all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','DMN','AUD'};
all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','DMN'};

nNet = length(all_networks_labels);

color = zeros(nNet,3);

color(1,:) = [0,0,128];
color(2,:) = [0.5,0,0.9];
color(3,:) = [0,191,255];
color(4,:) = [0,100,0];
color(5,:) = [255,255,0];
color(6,:) = [0.91,0.41,0.17];
color(7,:) = [255,0,0];
color(8,:) = [255,222,173];

for iNet=1:nNet

    network_label = all_networks_labels{iNet};
    
    network_seeds = getFunctionalSeeds_v5(all_networks_labels{iNet});
    
    nROI = length(network_seeds.ROI);
    
    for iROI=1:nROI
        
        seed_volume_filename = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-',network_seeds.ROI(iROI).label,'-clu');
        seed_volume = nifti(strcat(seed_volume_filename,'.nii'));
        seed_volume.dat.fname = strcat(settings.folders.main,'/',settings.folders.experiment,'/',settings.folders.subject,'/','output','/','Functional_Parcellation','/','after_correction','/','Seeds_Volumes_Clusters','/',seed_volume_filename,'.nii');
        
        net(iNet).ROI(iROI).img = seed_volume.dat(:,:,:);
        net(iNet).ROI(iROI).idx_voxels = find(seed_volume.dat(:,:,:));
        if size(net(iNet).ROI(iROI).idx_voxels,1) ~= 1, net(iNet).ROI(iROI).idx_voxels = net(iNet).ROI(iROI).idx_voxels'; end
    
    end
    
end

standard_volume_filename = 'avg152T1.nii';
standard_volume = nifti(standard_volume_filename);
standard_volume.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\spm8\canonical\',standard_volume_filename);
standard = standard_volume.dat(:,:,:);

all_voxels = [];
nAllROI = 0;
for iNet=1:nNet
    
    nROI = length(net(iNet).ROI);
    
    for iROI=1:nROI
        
        nAllROI = nAllROI + 1;
        
        all_voxels = [all_voxels,net(iNet).ROI(iROI).idx_voxels];
        
    end
    
end

unique_voxels = unique(all_voxels(:));

ParcellationVolume = zeros(size(net(1).ROI(1).img));

img_2_color = zeros(size(net(1).ROI(1).img,2),size(net(1).ROI(1).img,1),3);
img_47_color = zeros(size(net(1).ROI(1).img,2),size(net(1).ROI(1).img,1),3);

standard_2 = zeros(size(net(1).ROI(1).img,1),size(net(1).ROI(1).img,2));
standard_47 = zeros(size(net(1).ROI(1).img,1),size(net(1).ROI(1).img,2));

standard_2(:,:) = standard(:,:,34);
standard_47(:,:) = standard(:,:,63);

standard_2 = permute(standard_2,[2 1]);
standard_47 = permute(standard_47,[2 1]);

standard_2 = flip(standard_2);
standard_47 = flip(standard_47);

for iVoxel=1:length(unique_voxels)
    
   [idxx,idxy,idxz] = ind2sub(size(net(1).ROI(1).img),unique_voxels(iVoxel));
    
   zstat = 0;
   
   for iNet=1:nNet
       
       nROI = length(net(iNet).ROI);
       
       for iROI=1:nROI
           
           new_zstat = net(iNet).ROI(iROI).img(idxx,idxy,idxz);
           
           if new_zstat > zstat
               
               zstat = new_zstat;
               
               ParcellationVolume(idxx,idxy,idxz) = iNet*100 + iROI;
               
               if idxz == 63, img_47_color(idxy,idxx,:) = color(iNet,:); end
               if idxz == 34, img_2_color(idxy,idxx,:) = color(iNet,:); end    
               
           end
           
       end
       
   end
   
end

img_47_color = imresize(flip(img_47_color),[918,956],'nearest');
img_2_color = imresize(flip(img_2_color),[918,956],'nearest');

standard_47 = imresize(standard_47,[918,956],'nearest');
standard_2 = imresize(standard_2,[918,956],'nearest');

imwrite(img_47_color,strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-+47-color','.png'));
imwrite(img_2_color,strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-+2-color','.png'));
   
imwrite(standard_47,'standard-+47.png');
imwrite(standard_2,'standard-+2.png');

h = imfuse(standard_47,img_47_color,'blend');
k = imfuse(standard_2,img_2_color,'blend');

imwrite(h,strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-+47-blend','.png'));
imwrite(k,strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-+2-blend','.png'));

close all;

nifti_file = seed_volume;
offset = seed_volume.dat.offset;
scl_slope = seed_volume.dat.scl_slope;
scl_inter = seed_volume.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = seed_volume.dat.dim;

descrip = analysis_label;

fname = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'.nii');
input_data = ParcellationVolume; 
real_save_image;


end

function paintParcellationVolumeSlices(settings)

experiment_label = settings.codes.experiment;
subject_label = settings.codes.subject;
analysis_label = 'Functional-Parcels';
analysis_step_label = 'Pop-Map';

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','DMN','AUD'};

nNet = length(all_networks_labels);

color = zeros(nNet,3);

color(1,:) = [0,0,128];
color(2,:) = [160,32,240];
color(3,:) = [0,191,255];
color(4,:) = [0,100,0];
color(5,:) = [255,255,0];
color(6,:) = [255,165,0];
color(7,:) = [255,0,0];
color(8,:) = [255,222,173];


img_2 = imread(strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-+2','.png'));
img_47 = imread(strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-+47','.png'));

img_2_2D = img_2(:,:,1);
img_47_2D = img_47(:,:,1);

img_2_color = img_2;
img_47_color = img_47;

for iNet=1:nNet

   network_seeds = getFunctionalSeeds_v5(all_networks_labels{iNet});
   nROI = length(network_seeds.ROI);

   for iROI=1:nROI

      idx_ROI = iNet*100 + iROI;
      
      idx_pixels = find(img_2_2D == idx_ROI);
      for idx=1:length(idx_pixels)
        [idxx,idxy] = ind2sub(size(img_2_2D),idx_pixels(idx));
        img_2_color(idxx,idxy,:) = color(iNet,:);
      end
      
      idx_pixels = find(img_47_2D == idx_ROI);
      for idx=1:length(idx_pixels)
        [idxx,idxy] = ind2sub(size(img_47_2D),idx_pixels(idx));
        img_47_color(idxx,idxy,:) = color(iNet,:);
      end

   end

end
 
imwrite(img_47_color,strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-+47-color','.png'));
imwrite(img_2_color,strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-+2-color','.png'));

end

%%% GROUP LEVEL

function populationMapGroupLevel = computePopulationLevelMapsGroupLevel(all_settings,network_label)

disp('POPULATUION MAP GROUP LEVEL');
disp(network_label);

analysis_label = 'Functional-Parcels';
analysis_step_label = 'Corr-Map';

nSubjects = length(all_settings);
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    experiment_label = settings.codes.experiment;
    subject_label = settings.codes.subject;

    c = load(strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'.mat'));
    
    corMaps(iSubject).Correlation_Map = c.Correlation_Map;
    
end

%%% RANDOM EFFECTS ANALYSIS
disp('RANDOM EFFECTS ANALYSIS');

network_seeds = getFunctionalSeeds_v5(network_label);

nROI = length(network_seeds.ROI);
nRuns = length(corMaps(1).Correlation_Map);

%%% FIND UNIQUE VOXELS AMONG RUNS AND SUBJECTS
disp('FIND UNIQUE VOXELS AMONG RUNS AND SUBJECTS');

all_voxels_and_subjects_and_runs = [];

iirun = 0;

for iSubject=1:nSubjects
    
    for irun=1:nRuns
        
        iirun = iirun + 1;
   
        idx_voxels(iirun).idx = corMaps(iSubject).Correlation_Map(irun).maps.idx_voxels;
        all_voxels_and_subjects_and_runs = [all_voxels_and_subjects_and_runs;idx_voxels(iirun).idx(:)];
    
    end
    
end

unique_voxels = unique(all_voxels_and_subjects_and_runs);
nUniqueVoxels = length(unique_voxels);

%% ORGANIZE ALL VOXELS FOR ALL RUNS AND SUBJECTS
disp('ORGANIZE ALL VOXELS FOR ALL RUNS AND SUBJECTS');

for iROI=1:nROI
    
    iirun = 0;
    for iSubject=1:nSubjects

        for irun=1:nRuns

            iirun = iirun + 1;

            all_runs(iROI).run(iirun).rho_fisher = zeros(1,nUniqueVoxels);

        end

    end

end

for iROI=1:nROI

   iirun = 0;

   for iSubject=1:nSubjects

       for irun=1:nRuns

           iirun = iirun + 1;
           
           this_run_idx_voxels = corMaps(iSubject).Correlation_Map(irun).maps.idx_voxels;
           
           cross_idx = ismember(unique_voxels,this_run_idx_voxels);
           
           all_runs(iROI).run(iirun).rho_fisher(cross_idx) = corMaps(iSubject).Correlation_Map(irun).maps.ROI(iROI).rho_fisher(:);
        
       end

   end

end


%% ESTIMATE SAMPLE VARIANCE FOR EACH RUN/SUBJECT
disp('ESTIMATE SAMPLE VARIANCE FOR EACH RUN/SUBJECT');

for iROI=1:nROI
    
    iirun = 0;
    
    variance(iROI).sum_of_var = 0;
    
    for iSubject=1:nSubjects
    
        for irun=1:nRuns

            iirun = iirun + 1;
            
            variance(iROI).all_runs(iirun).var = nanvar(all_runs(iROI).run(iirun).rho_fisher);
            
            variance(iROI).sum_of_var = variance(iROI).sum_of_var + variance(iROI).all_runs(iirun).var;

        end
    
    end
    
end

%% ESTIMATE SAMPLE VARIANCE FOR EACH VOXEL
disp('ESTIMATE SAMPLE VARIANCE FOR EACH VOXEL');

sample_variance = zeros(nROI,nUniqueVoxels);

for iVoxel=1:nUniqueVoxels
    
   for iROI=1:nROI
       
       iirun = 0;
       
       for iSubject=1:nSubjects
           
           for irun=1:nRuns
               
               iirun = iirun + 1;
               
               this_voxel(iirun) = all_runs(iROI).run(iirun).rho_fisher(iVoxel);
               
           end
           
       end
 
       sample_variance(iROI,iVoxel) = nanvar(this_voxel);
       
   end
   
end
       
%% ESTIMATE OF SIGMA TETA      
disp('ESTIMATE OF SIGMA TETA');

sigma_teta = zeros(nROI,nUniqueVoxels);
   
for iROI=1:nROI
    
    for iVoxel=1:nUniqueVoxels
    
        sigma_teta(iROI,iVoxel) = sample_variance(iROI,iVoxel) - ( variance(iROI).sum_of_var / (nRuns * nSubjects) );

        if sigma_teta(iROI,iVoxel) < 0; sigma_teta(iROI,iVoxel) = 0; end
   
    end
   
end

%% ESTIMATE WEIGHTED VALUE
disp('ESTIMATE WEIGHTED VALUE');

for iROI=1:nROI
    
    for iVoxel=1:nUniqueVoxels
        
        iirun = 0;
        
        for iSubject=1:nSubjects
            
            for irun=1:nRuns
                
                iirun = iirun + 1;
       
                w_per_ROI(iROI).w(iirun) = 1 / ( variance(iROI).all_runs(iirun).var + sigma_teta(iROI,iVoxel) );
                
            end
            
        end
        
    end
    
end

%% ESTIMATE FINAL VALUE OF VOXEL 
disp('ESTIMATE FINAL VALUE OF VOXEL'); 

all_teta_star = zeros(nROI,nUniqueVoxels);

rho_t(1:nROI) = struct('voxels_t',zeros(1,nUniqueVoxels),'voxels_z',zeros(1,nUniqueVoxels));

for iROI=1:nROI
    
    for iVoxel=1:nUniqueVoxels
        
        iirun = 0;
        
        for iSubject=1:nSubjects
            
            for iRun=1:nRuns
                
                iirun = iirun + 1;
                
                this_voxel(iirun) = all_runs(iROI).run(iirun).rho_fisher(iVoxel);
                
            end
            
        end
        
        all_teta_star(iROI,iVoxel) = sum(w_per_ROI(iROI).w.*this_voxel) / sum(w_per_ROI(iROI).w) ;
        
        rho_t(iROI).voxels_t(iVoxel) = all_teta_star(iROI,iVoxel) / sqrt( 1 / sum(w_per_ROI(iROI).w) );
        
    end
    
end

for iROI=1:nROI
    
    my_std = nanstd(rho_t(iROI).voxels_t);
    my_mean = nanmean(rho_t(iROI).voxels_t);

    rho_t(iROI).voxels_z = ( rho_t(iROI).voxels_t - my_mean ) / my_std;

end

populationGroupMap.network_label = network_label;
populationGroupMap.nROI = nROI;
populationGroupMap.nRuns = nRuns;
populationGroupMap.unique_voxels = unique_voxels;
populationGroupMap.variance = variance;
populationGroupMap.sample_variance = sample_variance;
populationGroupMap.sigma_teta = sigma_teta;
populationGroupMap.all_teta_star = all_teta_star;
populationGroupMap.rho_t = rho_t;

analysis_step_label = 'Pop-Map';

disp('SAVE');

save(strcat(experiment_label,'-','All-Subjects','-',analysis_label,'-',analysis_step_label,'-',network_label,'.mat'),'populationGroupMap');

end

function getVolumesFromPopulationLevelMapsGroupLevel(analysis_label,network_label,network_seeds)

experiment_label = 'LHR';
subject_label = 'All-Subjects';
analysis_step_label = 'Pop-Map';

load(strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'.mat'));

MNI_dim = [91 109 91];
nROI = populationGroupMap.nROI;

load_aal = nifti('ROI_MNI_V4.nii');

common_voxels = populationGroupMap.unique_voxels;

for iROI=1:nROI
   
    disp(network_seeds.ROI(iROI).label);

    seed_brain = zeros(MNI_dim);
    
    for iVoxel=1:length(common_voxels)
        
        [idxx,idxy,idxz] = ind2sub(MNI_dim,common_voxels(iVoxel));
       
        seed_brain(idxx,idxy,idxz) = populationGroupMap.rho_t(iROI).voxels_z(iVoxel);
        
    end

    nifti_file = load_aal;
    offset = load_aal.dat.offset;
    scl_slope = load_aal.dat.scl_slope;
    scl_inter = load_aal.dat.scl_inter;

    dtype = 'FLOAT32';
    offset = 0;

    dim = load_aal.dat.dim;

    descrip = network_label;

    fname = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-',strrep(network_seeds.ROI(iROI).label,'/','-'),'.nii');
    input_data = seed_brain; 
    real_save_image;
    
end

end

function computeParcellationVolumeGroupLevel

experiment_label = 'LHR';
subject_label = 'All-Subjects';
analysis_label = 'Functional-Parcels';
analysis_step_label = 'Pop-Map';
seed_volume_folder = 'Z:\Dropbox (Uni Magdeburg)\_DATA\Parcellation\Functional_Parcellation\v4\Seeds_Clu';

%all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','DMN','AUD'};
all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','DMN'};

nNet = length(all_networks_labels);

color = zeros(nNet,3);

color(1,:) = [0,0,128];
color(2,:) = [0.5,0,0.9];
color(3,:) = [0,191,255];
color(4,:) = [0,100,0];
color(5,:) = [255,255,0];
color(6,:) = [0.91,0.41,0.17];
color(7,:) = [255,0,0];
color(8,:) = [255,222,173];

for iNet=1:nNet

    network_label = all_networks_labels{iNet};
    
    network_seeds = getFunctionalSeeds_v5(all_networks_labels{iNet});
    
    nROI = length(network_seeds.ROI);
    
    for iROI=1:nROI
        
        seed_volume_filename = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-',strrep(network_seeds.ROI(iROI).label,'/','-'),'-clu');
        seed_volume = nifti(strcat(seed_volume_filename,'.nii'));
        seed_volume.dat.fname = strcat(seed_volume_folder,'/',seed_volume_filename,'.nii');
        
        net(iNet).ROI(iROI).img = seed_volume.dat(:,:,:);
        net(iNet).ROI(iROI).idx_voxels = find(seed_volume.dat(:,:,:));
        if size(net(iNet).ROI(iROI).idx_voxels,1) ~= 1, net(iNet).ROI(iROI).idx_voxels = net(iNet).ROI(iROI).idx_voxels'; end
    
    end
    
end

standard_volume_filename = 'MNI152_T1_2mm.nii';
standard_volume = nifti(standard_volume_filename);
standard_volume.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\FSL\',standard_volume_filename);
standard = standard_volume.dat(:,:,:);

all_voxels = [];
nAllROI = 0;
for iNet=1:nNet
    
    nROI = length(net(iNet).ROI);
    
    for iROI=1:nROI
        
        nAllROI = nAllROI + 1;
        
        all_voxels = [all_voxels,net(iNet).ROI(iROI).idx_voxels];
        
    end
    
end

unique_voxels = unique(all_voxels(:));

ParcellationVolume = zeros(size(net(1).ROI(1).img));

img_2_color = zeros(size(net(1).ROI(1).img,2),size(net(1).ROI(1).img,1),3);
img_47_color = zeros(size(net(1).ROI(1).img,2),size(net(1).ROI(1).img,1),3);

standard_2 = zeros(size(net(1).ROI(1).img,1),size(net(1).ROI(1).img,2));
standard_47 = zeros(size(net(1).ROI(1).img,1),size(net(1).ROI(1).img,2));

standard_2(:,:) = standard(:,:,34);
standard_47(:,:) = standard(:,:,63);

standard_2 = permute(standard_2,[2 1]);
standard_47 = permute(standard_47,[2 1]);

standard_2 = flip(standard_2);
standard_47 = flip(standard_47);

for iVoxel=1:length(unique_voxels)
    
   [idxx,idxy,idxz] = ind2sub(size(net(1).ROI(1).img),unique_voxels(iVoxel));
    
   zstat = 0;
   
   for iNet=1:nNet
       
       nROI = length(net(iNet).ROI);
       
       network_seeds = getFunctionalSeeds_v5(all_networks_labels{iNet});
       
       for iROI=1:nROI
           
           new_zstat = net(iNet).ROI(iROI).img(idxx,idxy,idxz);
           
           if new_zstat > zstat
               
               zstat = new_zstat;
               
               ParcellationVolume(idxx,idxy,idxz) = network_seeds.ROI(iROI).idx;
               
               if idxz == 63, img_47_color(idxy,idxx,:) = color(iNet,:); end
               if idxz == 34, img_2_color(idxy,idxx,:) = color(iNet,:); end    
               
           end
           
       end
       
   end
   
end

img_47_color = imresize(flip(img_47_color),[918,956],'nearest');
img_2_color = imresize(flip(img_2_color),[918,956],'nearest');

standard_47 = imresize(standard_47,[918,956],'nearest');
standard_2 = imresize(standard_2,[918,956],'nearest');

imwrite(img_47_color,strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-+47-color','.png'));
imwrite(img_2_color,strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-+2-color','.png'));
   
imwrite(standard_47,'standard-+47.png');
imwrite(standard_2,'standard-+2.png');

h = imfuse(standard_47,img_47_color,'blend');
k = imfuse(standard_2,img_2_color,'blend');

imwrite(h,strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-+47-blend','.png'));
imwrite(k,strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-+2-blend','.png'));

close all;

nifti_file = seed_volume;
offset = seed_volume.dat.offset;
scl_slope = seed_volume.dat.scl_slope;
scl_inter = seed_volume.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = seed_volume.dat.dim;

descrip = analysis_label;

fname = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'.nii');
input_data = ParcellationVolume; 
real_save_image;

end

function computeParcellationVolumeGroupLevelPerNet(all_networks_labels)

experiment_label = 'LHR';
subject_label = 'All-Subjects';
analysis_label = 'Functional-Parcels';
analysis_step_label = 'Pop-Map';
%seed_volume_folder = 'Z:\Dropbox (Uni Magdeburg)\_DATA\Parcellation\Functional_Parcellation\v4\Seeds_Clu';
seed_volume_folder = 'Z:\Dropbox (Uni Magdeburg)\__tmp\functional_parcellation_v5';

nNet = length(all_networks_labels);

color = zeros(nNet,3);

color(1,:) = [0,0,128];
color(2,:) = [0.5,0,0.9];
color(3,:) = [0,191,255];
color(4,:) = [0,100,0];
color(5,:) = [255,255,0];
color(6,:) = [0.91,0.41,0.17];
color(7,:) = [255,0,0];
color(8,:) = [255,222,173];

for iNet=1:nNet

    network_label = all_networks_labels{iNet};
    
    network_seeds = getFunctionalSeeds_v5(all_networks_labels{iNet});
    
    nROI = length(network_seeds.ROI);
    
    for iROI=1:nROI
        
        seed_volume_filename = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-',strrep(network_seeds.ROI(iROI).label,'/','-'),'-clu');
        seed_volume = nifti(strcat(seed_volume_filename,'.nii'));
        seed_volume.dat.fname = strcat(seed_volume_folder,'/',seed_volume_filename,'.nii');
        
        net(iNet).ROI(iROI).img = seed_volume.dat(:,:,:);
        net(iNet).ROI(iROI).idx_voxels = find(seed_volume.dat(:,:,:));
        if size(net(iNet).ROI(iROI).idx_voxels,1) ~= 1, net(iNet).ROI(iROI).idx_voxels = net(iNet).ROI(iROI).idx_voxels'; end
    
    end
    
end

standard_volume_filename = 'MNI152_T1_2mm.nii';
standard_volume = nifti(standard_volume_filename);
standard_volume.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\FSL\',standard_volume_filename);
standard = standard_volume.dat(:,:,:);

all_voxels = [];
nAllROI = 0;
for iNet=1:nNet
    
    nROI = length(net(iNet).ROI);
    
    for iROI=1:nROI
        
        nAllROI = nAllROI + 1;
        
        all_voxels = [all_voxels,net(iNet).ROI(iROI).idx_voxels];
        
    end
    
end

unique_voxels = unique(all_voxels(:));

for iNet=1:nNet
    
    img_2_color = zeros(size(net(1).ROI(1).img,2),size(net(1).ROI(1).img,1),3);
    img_47_color = zeros(size(net(1).ROI(1).img,2),size(net(1).ROI(1).img,1),3);

    standard_2 = zeros(size(net(1).ROI(1).img,1),size(net(1).ROI(1).img,2));
    standard_47 = zeros(size(net(1).ROI(1).img,1),size(net(1).ROI(1).img,2));

    standard_2(:,:) = standard(:,:,34);
    standard_47(:,:) = standard(:,:,63);

    standard_2 = permute(standard_2,[2 1]);
    standard_47 = permute(standard_47,[2 1]);

    standard_2 = flip(standard_2);
    standard_47 = flip(standard_47);

   ParcellationVolume_zstat = zeros(size(net(1).ROI(1).img));
   ParcellationVolume_seed = zeros(size(net(1).ROI(1).img));
   
   network_seeds = getFunctionalSeeds_v5(all_networks_labels{iNet});
   
   network_label = network_seeds.network_label;
   
   for iVoxel=1:length(unique_voxels)
       
       [idxx,idxy,idxz] = ind2sub(size(net(1).ROI(1).img),unique_voxels(iVoxel));
       
       nROI = length(net(iNet).ROI);
       
       zstat = 0;
       
       for iROI=1:nROI
           
           new_zstat = net(iNet).ROI(iROI).img(idxx,idxy,idxz);
           
           if new_zstat > zstat
               
               zstat = new_zstat;
               
               ParcellationVolume_zstat(idxx,idxy,idxz) = zstat;
               ParcellationVolume_seed(idxx,idxy,idxz) = network_seeds.ROI(iROI).idx;
               
               if idxz == 63, img_47_color(idxy,idxx,:) = color(iNet,:); end
               if idxz == 34, img_2_color(idxy,idxx,:) = color(iNet,:); end    
               
           end
           
       end
       
   end
   
    img_47_color = imresize(flip(img_47_color),[918,956],'nearest');
    img_2_color = imresize(flip(img_2_color),[918,956],'nearest');

    standard_47 = imresize(standard_47,[918,956],'nearest');
    standard_2 = imresize(standard_2,[918,956],'nearest');

    imwrite(img_47_color,strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-+47-color','.png'));
    imwrite(img_2_color,strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-+2-color','.png'));

    imwrite(standard_47,'standard-+47.png');
    imwrite(standard_2,'standard-+2.png');

    h = imfuse(standard_47,img_47_color,'blend');
    k = imfuse(standard_2,img_2_color,'blend');

    imwrite(h,strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-+47-blend','.png'));
    imwrite(k,strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-+2-blend','.png'));

    nifti_file = seed_volume;
    offset = seed_volume.dat.offset;
    scl_slope = seed_volume.dat.scl_slope;
    scl_inter = seed_volume.dat.scl_inter;

    dtype = 'FLOAT32';
    offset = 0;

    dim = seed_volume.dat.dim;

    descrip = analysis_label;

    fname = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-','zstat','.nii');
    input_data = ParcellationVolume_zstat; 
    real_save_image;
    
    fname = strcat(experiment_label,'-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_label,'-','seed','.nii');
    input_data = ParcellationVolume_seed; 
    real_save_image;
   
end

end

%%% MD758 and CORBETTA correspondence

function getMD758densityPerNET

nROIs = 90;
MNI_size = [91 109 91];
nTotalClusters = 758;

disp('...functional networks');

networks = {'DAN' 'VAN' 'VIS' 'AUD' 'LAN' 'FPC' 'SMN' 'DMN'};
folder_nets = 'Z:\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually\';

DAN = nifti(strcat(folder_nets,'DAN-bin.nii'));
DAN.dat.fname = strcat(folder_nets,'DAN-bin.nii');
DAN_img = DAN.dat(:,:,:);

VAN = nifti(strcat(folder_nets,'VAN-bin.nii'));
VAN.dat.fname = strcat(folder_nets,'VAN-bin.nii');
VAN_img = VAN.dat(:,:,:);

VIS = nifti(strcat(folder_nets,'VIS-bin.nii'));
VIS.dat.fname = strcat(folder_nets,'VIS-bin.nii');
VIS_img = VIS.dat(:,:,:);

AUD = nifti(strcat(folder_nets,'AUD-bin.nii'));
AUD.dat.fname = strcat(folder_nets,'AUD-bin.nii');
AUD_img = AUD.dat(:,:,:);

LAN = nifti(strcat(folder_nets,'LAN-bin.nii'));
LAN.dat.fname = strcat(folder_nets,'LAN-bin.nii');
LAN_img = LAN.dat(:,:,:);

FPC = nifti(strcat(folder_nets,'FPC-bin.nii'));
FPC.dat.fname = strcat(folder_nets,'FPC-bin.nii');
FPC_img = FPC.dat(:,:,:);

SMN = nifti(strcat(folder_nets,'SMN-bin.nii'));
SMN.dat.fname = strcat(folder_nets,'SMN-bin.nii');
SMN_img = SMN.dat(:,:,:);

DMN = nifti(strcat(folder_nets,'DMN-bin.nii'));
DMN.dat.fname = strcat(folder_nets,'DMN-bin.nii');
DMN_img = DMN.dat(:,:,:);

folder_fun = 'Z:\_DATA\Parcellation\758-Cluster\';

MD = nifti(strcat(folder_fun,'LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii'));
MD.dat.fname = strcat(folder_fun,'LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
MD_img = MD.dat(:,:,:);

load('Z:\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat');

mem_DAN = zeros(1,nTotalClusters);
mem_VAN = zeros(1,nTotalClusters);
mem_VIS = zeros(1,nTotalClusters);
mem_LAN = zeros(1,nTotalClusters);
mem_FPC = zeros(1,nTotalClusters);
mem_DMN = zeros(1,nTotalClusters);
mem_SMN = zeros(1,nTotalClusters);
mem_AUD = zeros(1,nTotalClusters);
 
iiVoxel = 0;
iiCluster = 0;
for iROI=1:nROIs

    nClusters = ROI(iROI).nClusters;

    for iCluster=1:nClusters

        iiCluster = iiCluster + 1;

        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

        nVoxels = length(idx_voxels);
        
        func_voxels(iiCluster) = nVoxels;

        for iVoxel=1:nVoxels

            iiVoxel = iiVoxel + 1;

            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

            if DAN_img(idxx,idxy,idxz)

                mem_DAN(iiCluster) = mem_DAN(iiCluster) + 1;
                
            end

            if VAN_img(idxx,idxy,idxz)

                mem_VAN(iiCluster) = mem_VAN(iiCluster) + 1;

            end
            
           if VIS_img(idxx,idxy,idxz)

                mem_VIS(iiCluster) = mem_VIS(iiCluster) + 1;
 
            end

            if LAN_img(idxx,idxy,idxz)

                mem_LAN(iiCluster) = mem_LAN(iiCluster) + 1;

            end

            if FPC_img(idxx,idxy,idxz)

                mem_FPC(iiCluster) = mem_FPC(iiCluster) + 1;

            end

            if DMN_img(idxx,idxy,idxz)

                mem_DMN(iiCluster) = mem_DMN(iiCluster) + 1;

            end

            if SMN_img(idxx,idxy,idxz)

                mem_SMN(iiCluster) = mem_SMN(iiCluster) + 1;

            end

            if AUD_img(idxx,idxy,idxz)

                mem_AUD(iiCluster) = mem_AUD(iiCluster) + 1;

            end
                            

        end

    end

end

save('MD758-Corbetta-density.mat','func_voxels','mem_DAN','mem_VAN','mem_VIS','mem_LAN','mem_FPC','mem_DMN','mem_AUD','mem_SMN');

end

function plotMD758densityPerNET

load('MD758-Corbetta-density-simple-model.mat');
nTotalClusters = 758;

% color(1,:) = [0,0,128]./255;
% color(2,:) = [0.5,0,0.9];
% color(6,:) = [0,191,255]./255;
% color(3,:) = [0,100,0]./255;
% color(5,:) = [255,255,0]./255;
% color(4,:) = [0.91,0.41,0.17];
% color(8,:) = [255,0,0]./255;
% color(7,:) = [255,222,173]./255;

color{1} = 'b';
color{2} = 'b';
color{3} = 'b';
color{4} = 'g';
color{5} = 'y';
color{6} = 'c';
color{7} = 'r';
color{8} = 'k';

subplot(4,2,1);
idx_net = 1;
idx = find(mem_DAN);
b = bar(1:nTotalClusters,mem_DAN,color{idx_net});
xlim([1 nTotalClusters]);

subplot(4,2,2);
idx_net = 2;
idx = find(mem_VAN);
b = bar(1:nTotalClusters,mem_VAN,color{idx_net});
xlim([1 nTotalClusters]);

subplot(4,2,3);
idx_net = 3;
idx = find(mem_SMN);
b = bar(1:nTotalClusters,mem_SMN,color{idx_net});
xlim([1 nTotalClusters]);

subplot(4,2,4);
idx_net = 4;
idx = find(mem_VIS);
b = bar(1:nTotalClusters,mem_VIS,color{idx_net});
xlim([1 nTotalClusters]);

subplot(4,2,5);
idx_net = 5;
idx = find(mem_FPC);
b = bar(1:nTotalClusters,mem_FPC,color{idx_net});
xlim([1 nTotalClusters]);

subplot(4,2,6);
idx_net = 6;
idx = find(mem_LAN);
b = bar(1:nTotalClusters,mem_LAN,color{idx_net});
xlim([1 nTotalClusters]);

subplot(4,2,7);
idx_net = 7;
idx = find(mem_DMN);
b = bar(1:nTotalClusters,mem_DMN,color{idx_net});
xlim([1 nTotalClusters]);

subplot(4,2,8);
idx_net = 8;
idx = find(mem_AUD);
b = bar(1:nTotalClusters,mem_AUD,color{idx_net});
xlim([1 nTotalClusters]);

end

function paintParcellationVolumeSlicesFromMySimpleModel

MNI_size = [91 109 91];

standard_volume_filename = 'avg152T1.nii';
standard_volume = nifti(standard_volume_filename);
standard_volume.dat.fname = strcat('Z:\_TOOLBOX\spm8\canonical\',standard_volume_filename);
standard = standard_volume.dat(:,:,:);

networks_volume_filename = 'Functional-Network-no-overlap-expected.nii';
networks_volume = nifti(networks_volume_filename);
networks_volume.dat.fname = strcat('Z:\',networks_volume_filename);
networks = networks_volume.dat(:,:,:);

img_2_color = zeros(MNI_size(2),MNI_size(1),3);
img_47_color = zeros(MNI_size(2),MNI_size(1),3);

standard_2 = zeros(MNI_size(1),MNI_size(2));
standard_47 = zeros(MNI_size(1),MNI_size(2));

standard_2(:,:) = standard(:,:,34);
standard_47(:,:) = standard(:,:,63);

standard_2 = permute(standard_2,[2 1]);
standard_47 = permute(standard_47,[2 1]);

standard_2 = flip(standard_2);
standard_47 = flip(standard_47);

all_networks_labels = {'DAN','VAN','VIS','LAN','FPC','SMN','AUD','DMN'};

nNet = length(all_networks_labels);

color = zeros(nNet,3);

color(1,:) = [0,0,128];
color(2,:) = [160,32,240];
color(3,:) = [0,100,0];
color(4,:) = [255,165,0];
color(5,:) = [255,255,0];
color(6,:) = [0,191,255];
color(7,:) = [255,222,173];
color(8,:) = [255,0,0];

idx_voxels = find(networks);

for iVoxel=1:length(idx_voxels)
    
    [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
    
    iNet = networks(idxx,idxy,idxz);
    
    if iNet > 0
               
        if idxz == 63, img_47_color(idxy,idxx,:) = color(iNet,:); end
        if idxz == 34, img_2_color(idxy,idxx,:) = color(iNet,:); end    
               
    end

end
 
img_47_color = imresize(flip(img_47_color),[918,956],'nearest');
img_2_color = imresize(flip(img_2_color),[918,956],'nearest');

standard_47 = imresize(standard_47,[918,956],'nearest');
standard_2 = imresize(standard_2,[918,956],'nearest');

imwrite(img_47_color,strcat('Functional-Networks','-+47-color','.png'));
imwrite(img_2_color,strcat('Functional-Networks','-+2-color','.png'));

h = imfuse(standard_47,img_47_color,'blend');
k = imfuse(standard_2,img_2_color,'blend');

imwrite(h,strcat('Functional-Networks','-+47-blend','.png'));
imwrite(k,strcat('Functional-Networks','-+2-blend','.png'));

end

function getMD758densityPerNETSimpleMode

nROIs = 90;
MNI_size = [91 109 91];
nTotalClusters = 758;

disp('...functional networks');

all_networks_labels = {'DAN','VAN','VIS','LAN','FPC','SMN','AUD','DMN'};
networks_volume_filename = 'Functional-Network-no-overlap-expected.nii';
networks_volume = nifti(networks_volume_filename);
networks_volume.dat.fname = strcat('Z:\',networks_volume_filename);
networks = networks_volume.dat(:,:,:);

folder_fun = 'Z:\_DATA\Parcellation\758-Cluster\';

MD = nifti(strcat(folder_fun,'LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii'));
MD.dat.fname = strcat(folder_fun,'LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
MD_img = MD.dat(:,:,:);

load('Z:\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat');

mem_DAN = zeros(1,nTotalClusters);
mem_VAN = zeros(1,nTotalClusters);
mem_VIS = zeros(1,nTotalClusters);
mem_LAN = zeros(1,nTotalClusters);
mem_FPC = zeros(1,nTotalClusters);
mem_DMN = zeros(1,nTotalClusters);
mem_SMN = zeros(1,nTotalClusters);
mem_AUD = zeros(1,nTotalClusters);
 
iiVoxel = 0;
iiCluster = 0;
for iROI=1:nROIs

    nClusters = ROI(iROI).nClusters;

    for iCluster=1:nClusters

        iiCluster = iiCluster + 1;

        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

        nVoxels = length(idx_voxels);
        
        func_voxels(iiCluster) = nVoxels;

        for iVoxel=1:nVoxels

            iiVoxel = iiVoxel + 1;

            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

            iNet = networks(idxx,idxy,idxz);
            
            switch iNet
                
                case 1
                    
                    mem_DAN(iiCluster) = mem_DAN(iiCluster) + 1;
                    
                case 2
                    
                    mem_VAN(iiCluster) = mem_VAN(iiCluster) + 1;
                    
                case 3
                    
                    mem_VIS(iiCluster) = mem_VIS(iiCluster) + 1;
                    
                case 4
                    
                    mem_LAN(iiCluster) = mem_LAN(iiCluster) + 1;
                    
                case 5
                    
                    mem_FPC(iiCluster) = mem_FPC(iiCluster) + 1;
                    
                case 6
                    
                    mem_SMN(iiCluster) = mem_SMN(iiCluster) + 1;
                    
                case 7
                    
                    mem_AUD(iiCluster) = mem_AUD(iiCluster) + 1;
                    
                case 8
                    
                    mem_DMN(iiCluster) = mem_DMN(iiCluster) + 1;
                    
            end
            
        end

    end

end

save('MD758-Corbetta-density-simple-model.mat','func_voxels','mem_DAN','mem_VAN','mem_VIS','mem_LAN','mem_FPC','mem_DMN','mem_AUD','mem_SMN');

end
