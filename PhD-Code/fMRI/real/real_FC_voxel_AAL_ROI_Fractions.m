function real_FC_voxel_AAL_ROI_Fractions

%settings_subj1_2210;
%settings_subj2_2610;
%settings_subj3_0311;
%settings_subj4_0211;
%settings_subj5_0211;
%settings_subj6_2411;

% settings_subj1_2210;
% all_settings(1).settings = settings;
% clear settings
% settings_subj2_2610;
% all_settings(2).settings = settings;
% clear settings
% %settings_subj3_0311;
% settings_subj4_0211;
% all_settings(3).settings = settings;
% clear settings
% settings_subj5_0211;
% all_settings(4).settings = settings;
% clear settings
% %settings_subj6_2411;

Condition1_label = 'Track';
Condition2_label = 'Passive';
%computeAllFractionsFCDifferencePositiveNegative(settings,Condition1_label,Condition2_label);
%generate3DimgFractions(settings,Condition1_label,Condition2_label);
%generate3DimgFractionsForAllRunsPerSubject(settings,Condition1_label,Condition2_label);
%generate3DimgFractionsForAllRunsALLSubjects(all_settings,Condition1_label,Condition2_label);

Condition1_label = 'Passive';
Condition2_label = 'RestingState';
%computeAllFractionsFCDifferencePositiveNegative(settings,Condition1_label,Condition2_label);
%generate3DimgFractions(settings,Condition1_label,Condition2_label);
%generate3DimgFractionsForAllRunsPerSubject(settings,Condition1_label,Condition2_label);
%generate3DimgFractionsForAllRunsALLSubjects(all_settings,Condition1_label,Condition2_label);

Condition1_label = 'Track';
Condition2_label = 'RestingState';
%computeAllFractionsFCDifferencePositiveNegative(settings,Condition1_label,Condition2_label);
%generate3DimgFractions(settings,Condition1_label,Condition2_label);
%generate3DimgFractionsForAllRunsPerSubject(settings,Condition1_label,Condition2_label);
%generate3DimgFractionsForAllRunsALLSubjects(all_settings,Condition1_label,Condition2_label);

% settings_subj1_2210;
% doEverything(settings);
% clear settings
% settings_subj2_2610;
% doEverything(settings);
% clear settings
settings_subj3_0311;
doEverything(settings);
clear settings
% settings_subj4_0211;
% doEverything(settings);
% clear settings
% settings_subj5_0211;
% doEverything(settings);
% clear settings
settings_subj6_2411;
doEverything(settings);
clear settings

end

function doEverything(settings)

Condition1_label = 'Track';
Condition2_label = 'Passive';

computeAllFractionsFCDifferencePositiveNegative(settings,Condition1_label,Condition2_label);
generate3DimgFractions(settings,Condition1_label,Condition2_label);
%generate3DimgFractionsForAllRunsPerSubject(settings,Condition1_label,Condition2_label);

Condition1_label = 'Passive';
Condition2_label = 'RestingState';

computeAllFractionsFCDifferencePositiveNegative(settings,Condition1_label,Condition2_label);
generate3DimgFractions(settings,Condition1_label,Condition2_label);
%generate3DimgFractionsForAllRunsPerSubject(settings,Condition1_label,Condition2_label);

Condition1_label = 'Track';
Condition2_label = 'RestingState';

computeAllFractionsFCDifferencePositiveNegative(settings,Condition1_label,Condition2_label);
generate3DimgFractions(settings,Condition1_label,Condition2_label);
%generate3DimgFractionsForAllRunsPerSubject(settings,Condition1_label,Condition2_label);

%generateFractionsDistributionsAllContrasts(settings);

%plotFractionsDistributionsAllContrasts(settings)

end

function computeAllFractionsFCDifferencePositiveNegative(settings,Condition1_label,Condition2_label)

disp(strcat(Condition1_label,'-',Condition2_label));

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
%load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nROI = 90;
idx_areas = 1:nROI;

shouldIDoCluster = 0;

for iROI=1:nROI
    
    tic 
    
    label_ROI = AAL_ROI(idx_areas(iROI)).Nom_L;
    area_label = strrep(label_ROI,'_','-');
    disp(label_ROI);
    
    Condition1 = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-',Condition1_label,'-',area_label,'.mat'));
    Condition2 = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-',Condition2_label,'-',area_label,'.mat'));
    
    for irun=1:4

        rho_cond1 = eval(strcat('Condition1.FC_Voxels.run(irun).rho_',Condition1_label));
        pval_cond1 = eval(strcat('Condition1.FC_Voxels.run(irun).pval_',Condition1_label));

        rho_cond2 = eval(strcat('Condition2.FC_Voxels.run(irun).rho_',Condition2_label));
        pval_cond2 = eval(strcat('Condition2.FC_Voxels.run(irun).pval_',Condition2_label)); 

        [C1increaseFraction, C2increaseFraction, C1decreaseFraction, C2decreaseFraction, C1changes, C2changes, C1increase, C2increase, C1decrease, C2decrease, Dfisher, C1fisher, C2fisher, C1sorted, C2sorted, target_order, C1_source_order, C2_source_order, source_order, Ncomponent, cluster_assignment, Ncluster] = getFCDifferencePositiveNegativeFraction(rho_cond1,pval_cond1,rho_cond2,pval_cond2,shouldIDoCluster);

        FC_Fraction.run(irun).C1increaseFraction = C1increaseFraction;
        FC_Fraction.run(irun).C2increaseFraction = C2increaseFraction;
        
        FC_Fraction.run(irun).C1decreaseFraction = C1decreaseFraction;
        FC_Fraction.run(irun).C2decreaseFraction = C2decreaseFraction;
        
    end
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction-',Condition1_label(1),'-',Condition2_label(1),'-',area_label,'.mat'),'source_order','area_label','target_order','C1_source_order','C2_source_order','FC_Fraction');
    
    clear FC_Fraction
    
    toc
    
end

end

function generate3DimgFractions(settings,Condition1_label,Condition2_label)

disp(strcat(Condition1_label,'-',Condition2_label));

idx_areas = 1:90;

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

%nROI = 90;

for irun=1:4
    
    C1increaseFraction3D = zeros(size(AAL_img));
    C1decreaseFraction3D = zeros(size(AAL_img));

    C2increaseFraction3D = zeros(size(AAL_img));
    C2decreaseFraction3D = zeros(size(AAL_img));
    
    for iROI=1:length(idx_areas)
    
        label_ROI = AAL_ROI(idx_areas(iROI)).Nom_L;
        area_label = strrep(label_ROI,'_','-');
        disp(label_ROI);

        idx_AAL = AAL_ROI(idx_areas(iROI)).ID;
        idx_voxels = find(AAL_img == idx_AAL);

        load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction-',Condition1_label(1),'-',Condition2_label(1),'-',area_label,'.mat'));

        for iVoxel=1:length(idx_voxels)

            [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

            C1increaseFraction3D(idxx,idxy,idxz) = FC_Fraction.run(irun).C1increaseFraction(iVoxel);
            C2increaseFraction3D(idxx,idxy,idxz) = FC_Fraction.run(irun).C2increaseFraction(iVoxel);
            
            C1decreaseFraction3D(idxx,idxy,idxz) = FC_Fraction.run(irun).C1decreaseFraction(iVoxel);
            C2decreaseFraction3D(idxx,idxy,idxz) = FC_Fraction.run(irun).C2decreaseFraction(iVoxel);

        end
        
        clear FC_Fraction
    
    end
    
   
    scl_slope = 1;
    scl_inter = 0; 
    dim = [size(AAL_img,1),size(AAL_img,2),size(AAL_img,3)];
    dtype = 'FLOAT32';
    offset = 0;
    nifti_file.mat = load_aal.mat;
    nifti_file.mat_intent = load_aal.mat_intent;
    nifti_file.mat0 = load_aal.mat0;
    nifti_file.mat0_intent = load_aal.mat0_intent;

    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction-',Condition1_label(1),'-',Condition2_label(1),'-',Condition1_label(1),'i-R',int2str(irun),'.nii');
    descrip = strcat(Condition1_label,'increase','Fraction3D');
    input_data = C1increaseFraction3D; 
    real_save_image;
    
    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction-',Condition1_label(1),'-',Condition2_label(1),'-',Condition1_label(1),'d-R',int2str(irun),'.nii');
    descrip = strcat(Condition1_label,'decrease','Fraction3D');
    input_data = C1decreaseFraction3D; 
    real_save_image;
    
    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction-',Condition1_label(1),'-',Condition2_label(1),'-',Condition2_label(1),'i-R',int2str(irun),'.nii');
    descrip = strcat(Condition2_label,'increase','Fraction3D');
    input_data = C2increaseFraction3D; 
    real_save_image;
    
    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction-',Condition1_label(1),'-',Condition2_label(1),'-',Condition2_label(1),'d-R',int2str(irun),'.nii');
    descrip = strcat(Condition2_label,'decrease','Fraction3D');
    input_data = C2decreaseFraction3D; 
    real_save_image;
    
end

end

function generate3DimgFractionsForAllRunsPerSubject(settings,Condition1_label,Condition2_label)

disp(strcat(Condition1_label,'-',Condition2_label));

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

%%% LOAD RUNS

measure = {strcat(Condition1_label(1),'d'), strcat(Condition1_label(1),'i'), strcat(Condition2_label(1),'d'), strcat(Condition2_label(1),'i')};

for imeasure=1:length(measure)
    
    disp(measure{imeasure});
    
    for irun=1:4
        filename = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction-',Condition1_label(1),'-',Condition2_label(1),'-',measure{imeasure},'-R',int2str(irun),'.nii');
        load_run(irun).run = nifti(filename); 
        load_run(irun).run.dat.fname = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\','output','\','FC_Voxels_AAL_ROI','\','3DFractions','\',Condition1_label,'-',Condition2_label,'\',filename);
        runs(irun).run = load_run(irun).run.dat(:,:,:);
    end

    allruns_img = applyRandomEffectsToAllRuns(runs,AAL_img);
    
    scl_slope = 1;
    scl_inter = 0; 
    dim = [size(AAL_img,1),size(AAL_img,2),size(AAL_img,3)];
    dtype = 'FLOAT32';
    offset = 0;
    nifti_file.mat = load_aal.mat;
    nifti_file.mat_intent = load_aal.mat_intent;
    nifti_file.mat0 = load_aal.mat0;
    nifti_file.mat0_intent = load_aal.mat0_intent;

    fname = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction-',Condition1_label(1),'-',Condition2_label(1),'-',measure{imeasure},'-AllRuns','.nii');
    descrip = 'Fractions3DAllRuns';
    input_data = allruns_img; 
    real_save_image;
    
end

end

function generate3DimgFractionsForAllRunsALLSubjects(all_settings,Condition1_label,Condition2_label)

disp(strcat(Condition1_label,'-',Condition2_label));

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

%%% LOAD RUNS

measure = {strcat(Condition1_label(1),'d'), strcat(Condition1_label(1),'i'), strcat(Condition2_label(1),'d'), strcat(Condition2_label(1),'i')};

for imeasure=1:length(measure)
    
    disp(measure{imeasure});
    
    iirun = 0;
    for isubject=1:length(all_settings)
        
        settings = all_settings(isubject).settings;
    
        for irun=1:4
            iirun = iirun + 1;
            filename = strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction-',Condition1_label(1),'-',Condition2_label(1),'-',measure{imeasure},'-R',int2str(irun),'.nii');
            load_run(iirun).run = nifti(filename); 
            load_run(iirun).run.dat.fname = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\','output','\','FC_Voxels_AAL_ROI','\','3DFractions','\',Condition1_label,'-',Condition2_label,'\',filename);
            runs(iirun).run = load_run(iirun).run.dat(:,:,:);
        end
        
    end

    allruns_img = applyRandomEffectsToAllRuns(runs,AAL_img);
    
    scl_slope = 1;
    scl_inter = 0; 
    dim = [size(AAL_img,1),size(AAL_img,2),size(AAL_img,3)];
    dtype = 'FLOAT32';
    offset = 0;
    nifti_file.mat = load_aal.mat;
    nifti_file.mat_intent = load_aal.mat_intent;
    nifti_file.mat0 = load_aal.mat0;
    nifti_file.mat0_intent = load_aal.mat0_intent;

    fname = strcat(settings.codes.experiment,'-',int2str(length(all_settings)),'-Subjects','-','FC-Voxels-AAL-ROI-Fraction-',Condition1_label(1),'-',Condition2_label(1),'-',measure{imeasure},'-AllRuns','.nii');
    descrip = 'Fractions3DAllRuns';
    input_data = allruns_img; 
    real_save_image;
    
end

end

function allruns_img = applyRandomEffectsToAllRuns(runs,AAL_img)

idx_voxels = find(AAL_img);

for irun=1:length(runs)
         
    voxels = zeros(1,length(idx_voxels));
    
    for iVoxel=1:length(idx_voxels)
    
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

        voxels(iVoxel) = squeeze(runs(irun).run(idxx,idxy,idxz));
        
    end
   
    runszstat(irun).voxels = zscore(voxels);

end

for irun=1:length(runszstat)
    
   runszstat(irun).var = nanvar(runszstat(irun).voxels); 
    
end

for iVoxel=1:length(idx_voxels)

    all_runs_voxel = 0;
    for irun=1:length(runszstat)
        
        all_runs_voxel(irun) = runszstat(irun).voxels(iVoxel);
        
    end
    
    sample_variance(iVoxel) = var(all_runs_voxel);
    
    my_sigma_teta = 0;
    for irun=1:length(runszstat)
        
        my_sigma_teta = sample_variance(iVoxel) - runszstat(irun).var/length(runszstat);
        
    end
        
    if my_sigma_teta < 0; my_sigma_teta = 0; end
    
    sigma_teta(iVoxel) = my_sigma_teta;
    
end

all_teta_star = zeros(1,length(idx_voxels));

for iVoxel=1:length(idx_voxels)
   
    for irun=1:length(runszstat)
       
        runszstat(irun).singleVoxel = runszstat(irun).voxels(iVoxel);
        
        w_run(irun).w = 1 / ( runszstat(irun).var + sigma_teta(iVoxel) );
        
    end
    
    total_w = 0;
    for irun=1:length(runszstat)
    
        total_w = total_w + w_run(irun).w;
    
    end
    
    teta_star = 0;
    for irun=1:length(runszstat)
        
        teta_star = teta_star + ( (runszstat(irun).singleVoxel * w_run(irun).w) ) / (total_w); 
        
    end
    
    all_teta_star(iVoxel) = teta_star;
        
    allruns_voxels_t(iVoxel) = teta_star / sqrt( 1 / (total_w) );
    
end
    
    my_std = nanstd(allruns_voxels_t);
    my_mean = nanmean(allruns_voxels_t);
    
    allruns_voxels_z = ( allruns_voxels_t - my_mean ) / my_std;
    
    allruns_img = zeros(size(AAL_img));
    
    for iVoxel=1:length(idx_voxels)
         
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
        
        allruns_img(idxx,idxy,idxz) = allruns_voxels_z(iVoxel);
        
    end
    
end

function generateFractionsDistributionsAllContrasts(settings)

Condition1_label = 'Track';
Condition2_label = 'Passive';
generateFractionsDistributions(settings,Condition1_label,Condition2_label);

Condition1_label = 'Passive';
Condition2_label = 'RestingState';
generateFractionsDistributions(settings,Condition1_label,Condition2_label);

Condition1_label = 'Track';
Condition2_label = 'RestingState';
generateFractionsDistributions(settings,Condition1_label,Condition2_label);

end

function generateFractionsDistributions(settings,Condition1_label,Condition2_label)

disp(strcat(Condition1_label,'-',Condition2_label));

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
%load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nROI = 90;
idx_areas = 1:nROI;

nRun = 4;

for irun=1:nRun
    
    all_fractions(irun).C1increaseFractions = [];
    all_fractions(irun).C1decreaseFractions = [];
    all_fractions(irun).C2increaseFractions = [];
    all_fractions(irun).C2decreaseFractions = [];

end

for iROI=1:nROI
    
    tic 
    
    label_ROI = AAL_ROI(idx_areas(iROI)).Nom_L;
    area_label = strrep(label_ROI,'_','-');
    disp(label_ROI);
    
    load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction','-',Condition1_label(1),'-',Condition2_label(1),'-',area_label,'.mat'));

    for irun=1:nRun
       
         all_fractions(irun).C1increaseFractions = [all_fractions(irun).C1increaseFractions(:);FC_Fraction.run(irun).C1increaseFraction(:)];
         all_fractions(irun).C1decreaseFractions = [all_fractions(irun).C1decreaseFractions(:);FC_Fraction.run(irun).C1decreaseFraction(:)];
         all_fractions(irun).C2increaseFractions = [all_fractions(irun).C2increaseFractions(:);FC_Fraction.run(irun).C2increaseFraction(:)];
         all_fractions(irun).C2decreaseFractions = [all_fractions(irun).C2decreaseFractions(:);FC_Fraction.run(irun).C2decreaseFraction(:)];
        
    end
    
    clear FC_Fraction
    
end

for irun=1:nRun
    
    disp('getting ksdensity C1 Increase');
    
    [probability_distribution(irun).C1increase, interval(irun).C1increase] = getProbabilityDistribution(all_fractions(irun).C1increaseFractions(:));

    disp('getting ksdensity C1 Decrease');
    
    [probability_distribution(irun).C1decrease, interval(irun).C1decrease] = getProbabilityDistribution(all_fractions(irun).C1decreaseFractions(:));
    
    disp('getting ksdensity C2 Increase');
    
    [probability_distribution(irun).C2increase, interval(irun).C2increase] = getProbabilityDistribution(all_fractions(irun).C2increaseFractions(:));
    
    disp('getting ksdensity C2 Decrease');
    
    [probability_distribution(irun).C2decrease, interval(irun).C2decrease] = getProbabilityDistribution(all_fractions(irun).C2decreaseFractions(:));

end


save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction','-',Condition1_label(1),'-',Condition2_label(1),'-','Fractions-Density','.mat'),'probability_distribution','interval');

end

function [probability_distribution, interval] = getProbabilityDistribution(data)

min_data = 0;
%max_data = max(data(:));
max_data = 1;

step = 0.001;
interval = min_data:step:max_data;

k = ksdensity(data(:),interval);

probability_distribution = k .* step;

end

function plotFractionsDistributions(settings,Condition1_label,Condition2_label)


load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction','-',Condition1_label(1),'-',Condition2_label(1),'-','Fractions-Density','.mat'));

f = figure;

for irun=1:4
    
   subplot(2,2,irun); 
   
   interval(irun).C1increase(probability_distribution(irun).C1increase < 0.01) = [];
   probability_distribution(irun).C1increase(probability_distribution(irun).C1increase < 0.01) = []; 
   probability_distribution(irun).C1increase(interval(irun).C1increase < 0.01) = [];
   interval(irun).C1increase(interval(irun).C1increase < 0.01) = [];
   plot(interval(irun).C1increase,probability_distribution(irun).C1increase,'b');
   
   hold on
   
   interval(irun).C1decrease(probability_distribution(irun).C1decrease < 0.01) = [];   
   probability_distribution(irun).C1decrease(probability_distribution(irun).C1decrease < 0.01) = []; 
   probability_distribution(irun).C1decrease(interval(irun).C1decrease < 0.01) = [];
   interval(irun).C1decrease(interval(irun).C1decrease < 0.01) = [];
   plot(interval(irun).C1decrease,probability_distribution(irun).C1decrease,'r');
   
   interval(irun).C2increase(probability_distribution(irun).C2increase < 0.01) = [];
   probability_distribution(irun).C2increase(probability_distribution(irun).C2increase < 0.01) = []; 
   probability_distribution(irun).C2increase(interval(irun).C2increase < 0.01) = [];
   interval(irun).C2increase(interval(irun).C2increase < 0.01) = [];
   plot(interval(irun).C2increase,probability_distribution(irun).C2increase,'k');
   
   interval(irun).C2decrease(probability_distribution(irun).C2decrease < 0.01) = [];
   probability_distribution(irun).C2decrease(probability_distribution(irun).C2decrease < 0.01) = []; 
   probability_distribution(irun).C2decrease(interval(irun).C2decrease < 0.01) = [];
   interval(irun).C2decrease(interval(irun).C2decrease < 0.01) = [];
   plot(interval(irun).C2decrease,probability_distribution(irun).C2decrease,'y');
   
   max_interval = max([interval(irun).C1increase(:);interval(irun).C1decrease(:);interval(irun).C2increase(:);interval(irun).C2decrease(:)]);
   max_probability = max([probability_distribution(irun).C1increase(:);probability_distribution(irun).C1decrease(:);probability_distribution(irun).C2increase(:);probability_distribution(irun).C2decrease(:)]);
   
   xlim([0.01 max_interval]);
   ylim([0.01 max_probability]);
    
   title(strcat('Run:',int2str(irun)));
   
end

print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction','-',Condition1_label(1),'-',Condition2_label(1),'-','Fractions-Density','.jpg'));
print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction','-',Condition1_label(1),'-',Condition2_label(1),'-','Fractions-Density','.eps'));

g = figure;

plot(1:10,'b');
hold on
plot(1:10,'r');
plot(1:10,'k');
plot(1:10,'y');
legend({strcat(Condition1_label,'-','Increase'),strcat(Condition1_label,'-','Decrease'),strcat(Condition2_label,'-','Increase'),strcat(Condition2_label,'-','Decrease')});  

print(g,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction','-',Condition1_label(1),'-',Condition2_label(1),'-','Fractions-Density-legend','.jpg'));
print(g,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-Fraction','-',Condition1_label(1),'-',Condition2_label(1),'-','Fractions-Density-legend','.eps'));

end

function plotFractionsDistributionsAllContrasts(settings)

Condition1_label = 'Track';
Condition2_label = 'Passive';
plotFractionsDistributions(settings,Condition1_label,Condition2_label);

Condition1_label = 'Passive';
Condition2_label = 'RestingState';
plotFractionsDistributions(settings,Condition1_label,Condition2_label);

Condition1_label = 'Track';
Condition2_label = 'RestingState';
plotFractionsDistributions(settings,Condition1_label,Condition2_label);

end