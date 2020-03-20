function real_methods_paper_v1_check

%%% check datasets Magdeburg and Petra

function plotFCSubject(idx_ROI)

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

area_label = strrep(AAL_ROI(idx_ROI).Nom_L,'_','-');

nRuns = 4;
nSubjects = 8;
r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

for iSet=1:nSubjects

    load(strcat('LHR-SUBJ',int2str(iSet),'-FC-Voxels-AAL-ROI-corr-RestingState-',area_label,'.mat'));
    
    for iRun=1:nRuns
        
        get_all.corr(iRun,:,:) = FC_Voxels.run(iRun).rho_RestingState(:,:);
    
        single_rho = get_all.corr(iRun,:,:);
    
        get_all.fisher(iRun,:,:) = (1/2).*log((1 + single_rho)./(1 - single_rho));

    end
    
    get_all.m_z = squeeze(mean(get_all.fisher,1));
    get_all.s_z = squeeze(std(get_all.fisher));
    
    plotFCgeneral(get_all.m_z,get_all.s_z,strcat('LHR-SUBJ',int2str(iSet),'-','FC-Mean-',area_label));
    
    info(iSet).IM_total = TotalInformationContent( get_all.m_z(:), get_all.s_z, z_threshold, cv_threshold );
    
    clear FC_Voxels
    
end

save(strcat('LHR','-','FC-Mean-',area_label),'info');

end

function plotFCSubjectPETRA(idx_ROI)

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

area_label = strrep(AAL_ROI(idx_ROI).Nom_L,'_','-');

nRuns = 3;
nSubjects = 9;
r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

all_settings = getAllSettingsPetra;

for iSet=1:nSubjects
    
    settings = all_settings(iSet).settings;

    load(strcat('PETRA-',settings.codes.subject,'-FC-Voxels-AAL-ROI-corr-Petra-RestingState-',area_label,'.mat'));
    
    for iRun=1:nRuns
        
        get_all.corr(iRun,:,:) = FC_Voxels.run(iRun).rho_RestingState(:,:);
    
        single_rho = get_all.corr(iRun,:,:);
    
        get_all.fisher(iRun,:,:) = (1/2).*log((1 + single_rho)./(1 - single_rho));

    end
    
    get_all.m_z = squeeze(mean(get_all.fisher,1));
    get_all.s_z = squeeze(std(get_all.fisher));
    
    plotFCgeneral(get_all.m_z,get_all.s_z,strcat('PETRA-',settings.codes.subject,'-','FC-Mean-',area_label));
    
    info(iSet).IM_total = TotalInformationContent( get_all.m_z(:), get_all.s_z, z_threshold, cv_threshold );
    
    clear FC_Voxels
    
end

save(strcat('PETRA','-','FC-Mean-',area_label),'info');

end

function plotFCSubjectHCP(idx_ROI)

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

area_label = strrep(AAL_ROI(idx_ROI).Nom_L,'_','-');

nRuns = 4;
nSubjects = 8;
r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

all_settings = getAllSettingsHCP;

for iSet=1:nSubjects
    
    settings = all_settings(iSet).settings;

    for iRun=1:nRuns

        load(strcat('HCP','-','SUBJ',int2str(iSet),'-FC-Voxels-AAL-ROI-corr-HCP-RestingState-',area_label,'-',int2str(iRun),'.mat'));
    
        get_all.corr(iRun,:,:) = FC_Voxels.run(iRun).rho_RestingState(:,:);
    
        single_rho = get_all.corr(iRun,:,:);
    
        get_all.fisher(iRun,:,:) = (1/2).*log((1 + single_rho)./(1 - single_rho));

    end
    
    get_all.m_z = squeeze(mean(get_all.fisher,1));
    get_all.s_z = squeeze(std(get_all.fisher));
    
    plotFCgeneral(get_all.m_z,get_all.s_z,strcat('HCP','-','SUBJ',int2str(iSet),'-','FC-Mean-',area_label));
    
    info(iSet).IM_total = TotalInformationContent( get_all.m_z(:), get_all.s_z(:), z_threshold, cv_threshold );
    
    clear FC_Voxels
    
end

save(strcat('PETRA','-','FC-Mean-',area_label),'info');

end

function plotDistributionFCAALLHRPETRAHCP

%%% GET AAL - MAGDEBURG

load('AAL_ROI_mean_run_corr.mat');

nROI = 90;
nTotalRuns = 32;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                pairs_rest_rho(iRun,iPair) = single_rho;
            
            end
            
        end
        
    end
    
end

m_pairs_MD = mean(pairs_rest_rho,1);
m_pairs_z_MD = mean(pairs_rest_z,1);

%%% GET AAL - PETRA

load('AAL_ROI_mean_run_corr-Petra.mat');

nROI = 90;
nTotalRuns = 27;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                pairs_rest_rho(iRun,iPair) = single_rho;
            
            end
            
        end
        
    end
    
end

m_pairs_PETRA = mean(pairs_rest_rho,1);
m_pairs_z_PETRA = mean(pairs_rest_z,1);

%%% GET AAL - HCP

load('AAL_ROI_mean_run_corr-HCP.mat');

nROI = 90;
nTotalRuns = 32;
nTR = 1200;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                pairs_rest_rho(iRun,iPair) = single_rho;
            
            end
            
        end
        
    end
    
end

m_pairs_HCP = mean(pairs_rest_rho,1);
m_pairs_z_HCP = mean(pairs_rest_z,1);

max_val = max([m_pairs_z_MD(:);m_pairs_z_PETRA(:);m_pairs_z_HCP(:)]);
min_val = min([m_pairs_z_MD(:);m_pairs_z_PETRA(:);m_pairs_z_HCP(:)]);

plotDistributionDensity3(min_val,max_val,m_pairs_z_MD,m_pairs_z_PETRA,m_pairs_z_HCP,'Magdeburg','PETRA','HCP','Distribution of Correlations','Correlation Fisher Transformed','Probability');

end

function RestingState = getRestingStateFromMD

nRuns = 4;

all_settings = getAllSettings;

nSubjects = length(all_settings);

iiRun = 0;
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;

    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;

    kind = 'RestingState';
    for irun=1:nRuns
        iiRun = iiRun + 1;
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
end

end

function RestingState = getRestingStateFromPETRA

run = 1;
nSegments = 3;
nTR = 150;

all_settings = getAllSettingsPetra;

iRun = 0;
for iSet=1:length(all_settings)
    
    settings = all_settings(iSet).settings;

    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;

    [residual, settings] = real_get_data_FSL_Petra(settings,run,file,get_at_this_preprocessed_step);
    
    for iSeg=1:nSegments
        
        iRun = iRun + 1;
        
        RestingState(iRun).run(:,:,:,:) = residual(:,:,:,(1+(iSeg-1)*nTR):(nTR+(iSeg-1)*nTR));
        
    end

end

end

function RestingState = getRestingStateFromHCP

nRuns = 4;

all_settings = getAllSettingsHCP;
nSubjects = 8;

iiRun = 0;
for iSet=1:nSubjects
    
    settings = all_settings(iSet).settings;
    
    disp(strcat('Subj:',int2str(iSet)));

    for iRun=1:nRuns
        
        iiRun = iiRun + 1;
        
        [RestingState(iiRun).run, settings] = real_get_data_HCP(settings,iRun);
    
    end
    
end

end

function RestingState = getRestingStateFromMARTIN

nRuns = 4;

settings_martin;

for iRun=1:nRuns

    disp(strcat('Run:',int2str(iRun)));
    
    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;

    [RestingState(iRun).run, settings] = real_get_data_FSL_Martin(settings,iRun,file,get_at_this_preprocessed_step);
    
end

end

function getMeanSTDVoxels(RestingState,label)

nRuns = length(RestingState);

nTotalVoxels = 160990;

MNI_size = [91 109 91];

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

iiVoxel = 0;

mean_voxels = zeros(nRuns,nTotalVoxels);
std_voxels = zeros(nRuns,nTotalVoxels);

for iROI=1:nNodes
    
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    nVoxels = length(idx_voxels);
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(strcat(int2str(iROI),':',area_label,':',int2str(length(idx_voxels)),':voxels'));
    
    for iVoxel=1:nVoxels
        
        iiVoxel = iiVoxel + 1;
        
        [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
        
        for iRun=1:nRuns
            
            voxel_ts = RestingState(iRun).run(idxx,idxy,idxz,:);
            
            m_voxel_ts = nanmean(voxel_ts);
            s_voxel_ts = nanstd(voxel_ts);
            
            mean_voxels(iRun,iiVoxel) = m_voxel_ts;
            std_voxels(iRun,iiVoxel) = s_voxel_ts;
            
        end
        
    end
    
end

save(strcat('Mean-STD-Voxels-AAL-',label,'.mat'),'mean_voxels','std_voxels','-v7.3');

end

function plotSNRRestingStateMDPETRAHCP

%%% SNR

% MD

load('Mean-STD-Voxels-AAL-MD.mat');

vector_one = std_voxels(:);
vector_two = mean_voxels(:);
title_label = 'SNR-MD';

MyContourCorrelation_v5(vector_one,vector_two,title_label,1);

clear mean_voxels std_voxels

% PETRA

load('Mean-STD-Voxels-AAL-PETRA.mat');

vector_one = std_voxels(:);
vector_two = mean_voxels(:);
title_label = 'SNR-PETRA';

MyContourCorrelation_v5(vector_one,vector_two,title_label,1);

clear mean_voxels std_voxels

% HCP

load('Mean-STD-Voxels-AAL-HCP.mat');

vector_one = std_voxels(:);
vector_two = mean_voxels(:);
title_label = 'SNR-HCP';

MyContourCorrelation_v5(vector_one,vector_two,title_label,1);

clear mean_voxels std_voxels

end

function plotSNRRestingStateMDPETRAHCPv2

%%% SNR

% MD

load('Mean-STD-Voxels-AAL-MD.mat');

v(1).vector_one = std_voxels(:);
v(1).vector_two = mean_voxels(:);
v(1).title_label = 'MD';

clear mean_voxels std_voxels

% PETRA

load('Mean-STD-Voxels-AAL-PETRA.mat');

v(2).vector_one = std_voxels(:);
v(2).vector_two = mean_voxels(:);
v(2).title_label = 'PETRA';

clear mean_voxels std_voxels

% HCP

load('Mean-STD-Voxels-AAL-HCP.mat');

v(3).vector_one = std_voxels(:);
v(3).vector_two = mean_voxels(:);
v(3).title_label = 'HCP';

mysimpleplot(v);

clear mean_voxels std_voxels

end

function mysimpleplot(v)

color{1} = 'r*';
color{2} = 'b*';
color{3} = 'k*';

f = figure;

for iV=1:length(v)
    
    plot(v(iV).vector_two./v(iV).vector_one,v(iV).vector_one,color{iV},'MarkerSize',0.3);
    
    hold on;
    
end

legend({v(1).title_label,v(2).title_label,v(3).title_label});

end

function SNRinaVolume(label)

load(strcat('Mean-STD-Voxels-AAL-',label,'.mat'));

m_mean = nanmean(mean_voxels,1);
m_std = nanmean(std_voxels,1);

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

MNI_size = [91 109 91];

vol_mean = zeros(MNI_size);
vol_std = zeros(MNI_size);

iiVoxel = 0;
for iROI=1:nNodes
    
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    nVoxels = length(idx_voxels);
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(strcat(int2str(iROI),':',area_label,':',int2str(length(idx_voxels)),':voxels'));
    
    for iVoxel=1:nVoxels
        
        iiVoxel = iiVoxel + 1;
        
        [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
        
        vol_mean(idxx,idxy,idxz) = m_mean(iiVoxel);
        
        vol_std(idxx,idxy,idxz) = m_std(iiVoxel);
        
    end
    
end

nifti_file = load_aal;
offset = load_aal;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = load_aal.dat.dim;

descrip = 'paint';

% fname = strcat('Mean-Voxels-AAL-',label,'.nii');
% input_data = vol_mean; 
% real_save_image;
% 
% fname = strcat('STD-Voxels-AAL-',label,'.nii');
% input_data = vol_std; 
% real_save_image;

fname = strcat('SNR-Voxels-AAL-',label,'.nii');
input_data = vol_mean ./ vol_std; 
% real_save_image;

input_data(input_data==0) = [];
min(input_data(:))

end

%%% check dataset Martin

function getMeanTimeSeriesFrom90AALRunByRunMARTIN

nROI = 90;
nTotalRuns = 4;
nTR = 150;

settings_martin;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for iRun=1:nTotalRuns
    
    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;

    [run(iRun).residual, settings] = real_get_data_FSL_Martin(settings,iRun,file,get_at_this_preprocessed_step);

end
    
for iRun=1:nTotalRuns

    disp(strcat('Run:',int2str(iRun)));

    for iROI=1:nROI

        idx_ROI = AAL_ROI(iROI).ID;
        idx_ROI_Label = AAL_ROI(iROI).Nom_L;
        idx_ROI_Label = strrep(idx_ROI_Label,'_','-');

        idx_voxels = find(AAL_img == idx_ROI);

        nVoxels = length(idx_voxels);

        area_RestingState = zeros(nVoxels,nTR);

        izr = 0;

        zeros_RestingState = [];

        for iVoxel=1:nVoxels

            [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

            area_RestingState(iVoxel,:) = run(iRun).residual(idxx,idxy,idxz,:);

            if sum(area_RestingState(iVoxel,:)) == 0

                izr = izr + 1;
                zeros_RestingState(izr) = iVoxel;

            end

        end

        area_RestingState(zeros_RestingState,:) = [];

        mn_RestingState = mean(area_RestingState,1);

        mean_area_RestingState_voxel(iROI,:) = mn_RestingState(:);

    end

    mean_run(iRun).rest = mean_area_RestingState_voxel;

    clear mean_area_RestingState_voxel

end 

save('AAL_ROI_mean_run-Martin.mat','mean_run');

end

function getCorrFromMeanTimeSeriesFrom90AALRunByRunMARTIN

load('AAL_ROI_mean_run-Martin.mat');

nTotalRuns = 4;

for iRun=1:nTotalRuns
    
   run = mean_run(iRun).rest;
   run = run';
   rho(iRun).rest_corr = corr(run);
    
end

save('AAL_ROI_mean_run_corr-Martin.mat','rho');

end

function getAndPlotVarianceDistributionFromFC90AALMARTIN

load('AAL_ROI_mean_run_corr-Martin.mat');

nROI = 90;
nTotalRuns = 4;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            end
            
        end
        
    end
    
end

my_var_aal = var(pairs_rest);
my_cv_aal = std(pairs_rest) ./ abs(mean(pairs_rest,1));
my_mean_aal = mean(pairs_rest,1);
my_std_aal = std(pairs_rest);

vector_one = my_std_aal;
vector_two = my_mean_aal;
vector_three = my_cv_aal;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'AAL-Rest-Martin';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

fraction = FractionSignificantConsistent( vector_two, vector_one, z_threshold, cv_threshold )

end

function getMeanMD758MARTIN

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nROI = 90;
nTotalRuns = 4;
nTR = 150;

settings_martin;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for iRun=1:nTotalRuns
    
    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;

    [run(iRun).residual, settings] = real_get_data_FSL_Martin(settings,iRun,file,get_at_this_preprocessed_step);

end

for iRun=1:nTotalRuns

    disp(strcat('Run:',int2str(iRun)));
    
    iiCluster = 0;

    for iROI=1:nROI
        
        nClusters = length(ROI(iROI).clusters);
        
        for iCluster=1:nClusters
            
            iiCluster = iiCluster + 1;

            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);

            area_RestingState = zeros(nVoxels,nTR);

            izr = 0;

            zeros_RestingState = [];

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

                area_RestingState(iVoxel,:) = run(iRun).residual(idxx,idxy,idxz,:);

                if sum(area_RestingState(iVoxel,:)) == 0

                    izr = izr + 1;
                    zeros_RestingState(izr) = iVoxel;

                end

            end

            area_RestingState(zeros_RestingState,:) = [];

            mn_RestingState = mean(area_RestingState,1);

            mean_area_RestingState_voxel(iiCluster,:) = mn_RestingState(:);

        end
        
    end

    mean_run(iRun).rest = mean_area_RestingState_voxel;

    clear mean_area_RestingState_voxel

end 

save('MD758_ROI_mean_run-Martin.mat','mean_run');

end

function getMeanMD758SUBJECT5Again

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nROI = 90;
nTotalRuns = 4;
nTR = 150;
    
%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for iRun=1:nTotalRuns
    
    settings_subj5_again;

    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    file = settings.FSL.files.functional.residual;
    mask = settings.FSL.files.mask.warped;    
    kind = 'RestingState';

    [run(iRun).residual, settings] = real_get_data_FSL(settings,kind,iRun,file,mask,get_at_this_preprocessed_step);

    clear settings
    
end

for iRun=1:nTotalRuns

    disp(strcat('Run:',int2str(iRun)));
    
    iiCluster = 0;

    for iROI=1:nROI
        
        nClusters = length(ROI(iROI).clusters);
        
        for iCluster=1:nClusters
            
            iiCluster = iiCluster + 1;

            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);

            area_RestingState = zeros(nVoxels,nTR);

            izr = 0;

            zeros_RestingState = [];

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

                area_RestingState(iVoxel,:) = run(iRun).residual(idxx,idxy,idxz,:);

                if sum(area_RestingState(iVoxel,:)) == 0

                    izr = izr + 1;
                    zeros_RestingState(izr) = iVoxel;

                end

            end

            area_RestingState(zeros_RestingState,:) = [];

            mn_RestingState = mean(area_RestingState,1);

            mean_area_RestingState_voxel(iiCluster,:) = mn_RestingState(:);

        end
        
    end

    mean_run(iRun).rest = mean_area_RestingState_voxel;

    clear mean_area_RestingState_voxel

end 

save('MD758_ROI_mean_run-SUBJ5-Again.mat','mean_run');

end

function getCorrMD758MARTIN

load('MD758_ROI_mean_run-Martin.mat');

nTotalRuns = 4;

for iRun=1:nTotalRuns
    
   run = mean_run(iRun).rest;
   run = run';
   rho(iRun).rest_corr = corr(run);
    
end

save('MD758_ROI_mean_run_corr-Martin.mat','rho');

end

function getCorrMD758SUBJ5Again

load('MD758_ROI_mean_run-SUBJ5-Again.mat');

nTotalRuns = 4;

for iRun=1:nTotalRuns
    
   run = mean_run(iRun).rest;
   run = run';
   rho(iRun).rest_corr = corr(run);
    
end

save('MD758_ROI_mean_run_corr-SUBJ5-Again.mat','rho');

end

function getAndPlotVarianceMD758MARTIN

load('MD758_ROI_mean_run_corr-Martin.mat');

nROI = 758;
nTotalRuns = 4;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            end
            
        end
        
    end
    
end

my_var_aal = var(pairs_rest);
my_cv_aal = std(pairs_rest) ./ abs(mean(pairs_rest,1));
my_mean_aal = mean(pairs_rest,1);
my_std_aal = std(pairs_rest);

vector_one = my_std_aal;
vector_two = my_mean_aal;
vector_three = my_cv_aal;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'MD758-Rest-Martin';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

fraction = FractionSignificantConsistent( vector_two, vector_one, z_threshold, cv_threshold )

end

function plotMeanFCAALMARTIN

load('AAL_ROI_mean_run_corr-Martin.mat');

nROI = 90;
nTotalRuns = 4;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                pairs_rest_rho(iRun,iPair) = single_rho;
            
            end
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z,1);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ mean(pairs_rest_z,1);



mean_corr = zeros(nROI,nROI);
s_corr = zeros(nROI,nROI);

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            %if my_cv_aal(iPair) < cv_criterion & my_cv_aal(iPair) > 0 & p_pairs(iPair) < p_criterion
                
                mean_corr(iROI,iiROI) =  m_pairs_z(iPair);
                s_corr(iROI,iiROI) =  s_pairs_z(iPair);
                
            %end
            
        end
        
    end
    
end

label = 'FC-AAL-Martin';
plotFCgeneral(mean_corr,s_corr,label);

end

function plotMeanFCMD758MARTIN

load('MD758_ROI_mean_run_corr-Martin.mat');

nROI = 758;
nTotalRuns = 4;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                pairs_rest_rho(iRun,iPair) = single_rho;
            
            end
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z,1);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ mean(pairs_rest_z,1);



mean_corr = zeros(nROI,nROI);
s_corr = zeros(nROI,nROI);

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            %if my_cv_aal(iPair) < cv_criterion & my_cv_aal(iPair) > 0 & p_pairs(iPair) < p_criterion
                
                mean_corr(iROI,iiROI) =  m_pairs_z(iPair);
                s_corr(iROI,iiROI) =  s_pairs_z(iPair);
                
            %end
            
        end
        
    end
    
end

label = 'FC-MD758-Martin';
plotFCgeneral(mean_corr,s_corr,label);

end

function plotMeanFCMD758SUBJ5Again

load('MD758_ROI_mean_run_corr-SUBJ5-Again.mat');

nROI = 758;
nTotalRuns = 4;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                pairs_rest_rho(iRun,iPair) = single_rho;
            
            end
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z,1);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ mean(pairs_rest_z,1);



mean_corr = zeros(nROI,nROI);
s_corr = zeros(nROI,nROI);

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            %if my_cv_aal(iPair) < cv_criterion & my_cv_aal(iPair) > 0 & p_pairs(iPair) < p_criterion
                
                mean_corr(iROI,iiROI) =  m_pairs_z(iPair);
                s_corr(iROI,iiROI) =  s_pairs_z(iPair);
                
            %end
            
        end
        
    end
    
end

label = 'FC-MD758-SUBJ5-Again';
plotFCgeneral(mean_corr,s_corr,label);

end


%%% COMPARE PREPROCESS AND FMRI PROTOCOL - MARTIN DATASET

function [PETRA_custom, PETRA_melodic, MAGDEBURG_custom, MAGDEBURG_melodic] = loadAllDataMartinDataset

nTotalRuns = 4;
start_TR = 17;
end_TR = 166;

folder_PETRA_protocol = 'Z:\Dropbox (Uni Magdeburg)\_DATA\MARTIN\preprocessed\T2-Stimulus-RestingState\dn20_1590';
folder_MAGDEBURG_protocol = 'Z:\Dropbox (Uni Magdeburg)\_DATA\MARTIN\preprocessed\T2-Stimulus-RestingState\dn20_1622';
file_custom = 'FSL\custom\filtered_func_data_mcf_unwarp2standard-clean-voxel-res';
file_melodic = 'FSL\Melodic-Fieldmap.ica\filtered_func_data2standard-clean-voxel-res';

for iRun=1:nTotalRuns
    
    disp(strcat('Run:',int2str(iRun)));
    
    disp('PETRA - custom');

    total_file_path = strcat(folder_PETRA_protocol,'\Run-',int2str(iRun),'\',file_custom,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    PETRA_custom(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
    
    disp('PETRA - melodic');
    
    total_file_path = strcat(folder_PETRA_protocol,'\Run-',int2str(iRun),'\',file_melodic,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    PETRA_melodic(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
    
    disp('MAGDEBURG - custom');
    
    total_file_path = strcat(folder_MAGDEBURG_protocol,'\Run-',int2str(iRun),'\',file_custom,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    MAGDEBURG_custom(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
    
    disp('MAGDEBURG - melodic');
    
    total_file_path = strcat(folder_MAGDEBURG_protocol,'\Run-',int2str(iRun),'\',file_melodic,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    MAGDEBURG_melodic(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri

end

end

function [SUBJ5AgainCustom] = loadAllDataSUBJ5Dataset

nTotalRuns = 4;
start_TR = 17;
end_TR = 166;

folder_SUBJ5Again_protocol = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-5-preprocess-test\preprocessed\T2-Stimulus-RestingState';
file_custom = 'FSL\custom\filtered_func_data_mcf_unwarp2standard-clean-voxel-res';

    iRun = 1;
    total_file_path = strcat(folder_SUBJ5Again_protocol,'\Run-1-4-1','\',file_custom,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    SUBJ5AgainCustom(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
    
    iRun = 2;
    total_file_path = strcat(folder_SUBJ5Again_protocol,'\Run-2-8-1','\',file_custom,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    SUBJ5AgainCustom(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
    
    iRun = 3;
    total_file_path = strcat(folder_SUBJ5Again_protocol,'\Run-3-4-2','\',file_custom,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    SUBJ5AgainCustom(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
    
    iRun = 4;
    total_file_path = strcat(folder_SUBJ5Again_protocol,'\Run-4-8-2','\',file_custom,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    SUBJ5AgainCustom(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
   

end

function [SUBJ5AgainCustom] = loadAllDataSUBJ5DatasetPreMask

nTotalRuns = 4;
start_TR = 17;
end_TR = 166;

folder_SUBJ5Again_protocol = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-5-preprocess-test\preprocessed\T2-Stimulus-RestingState';
file_custom = 'FSL\custom\filtered_func_data_mcf_unwarp2standard-clean-voxel-pre-res';

    iRun = 1;
    total_file_path = strcat(folder_SUBJ5Again_protocol,'\Run-1-4-1','\',file_custom,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    SUBJ5AgainCustom(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
    
    iRun = 2;
    total_file_path = strcat(folder_SUBJ5Again_protocol,'\Run-2-8-1','\',file_custom,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    SUBJ5AgainCustom(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
    
    iRun = 3;
    total_file_path = strcat(folder_SUBJ5Again_protocol,'\Run-3-4-2','\',file_custom,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    SUBJ5AgainCustom(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
    
    iRun = 4;
    total_file_path = strcat(folder_SUBJ5Again_protocol,'\Run-4-8-2','\',file_custom,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    SUBJ5AgainCustom(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
   

end

function [SUBJ5AgainCustom] = loadAllDataSUBJ5DatasetOriginal

nTotalRuns = 4;
start_TR = 17;
end_TR = 166;

folder_SUBJ5Again_protocol = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-5-2-11-2015\preprocessed\T2-Stimulus-RestingState';
% file_custom = 'FSL\custom\filtered_func_data_mcf_unwarp2standard-clean-voxel-res';
file_custom = 'FSL\custom\filtered_func_data_mcf_unwarp2standard';

    iRun = 1;
    total_file_path = strcat(folder_SUBJ5Again_protocol,'\Run-1-4-1','\',file_custom,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    SUBJ5AgainCustom(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
    
    iRun = 2;
    total_file_path = strcat(folder_SUBJ5Again_protocol,'\Run-2-8-1','\',file_custom,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    SUBJ5AgainCustom(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
    
    iRun = 3;
    total_file_path = strcat(folder_SUBJ5Again_protocol,'\Run-3-4-2','\',file_custom,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    SUBJ5AgainCustom(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
    
    iRun = 4;
    total_file_path = strcat(folder_SUBJ5Again_protocol,'\Run-4-8-2','\',file_custom,'.nii');
    load_fmri = nifti(total_file_path);
    load_fmri.dat.fname = total_file_path;
    SUBJ5AgainCustom(iRun).run = load_fmri.dat(:,:,:,start_TR:end_TR);
    clear load_fmri
   

end

function getMeanMD758MARTINbothScans(run,label)

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nROI = 90;
nTotalRuns = length(run);
nTR = 150;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for iRun=1:nTotalRuns

    disp(strcat('Run:',int2str(iRun)));
    
    iiCluster = 0;

    for iROI=1:nROI
        
        nClusters = length(ROI(iROI).clusters);
        
        for iCluster=1:nClusters
            
            iiCluster = iiCluster + 1;

            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);

            area_RestingState = zeros(nVoxels,nTR);

            %izr = 0;

            %zeros_RestingState = [];

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

                area_RestingState(iVoxel,:) = run(iRun).run(idxx,idxy,idxz,:);

                %if sum(area_RestingState(iVoxel,:)) == 0

                    %izr = izr + 1;
                    %zeros_RestingState(izr) = iVoxel;

                %end

            end

            %area_RestingState(zeros_RestingState,:) = [];

            mn_RestingState = nanmean(area_RestingState,1);

            mean_area_RestingState_voxel(iiCluster,:) = mn_RestingState(:);

        end
        
    end

    mean_run(iRun).rest = mean_area_RestingState_voxel;

    clear mean_area_RestingState_voxel

end 

save(strcat('MD758_ROI_mean_run-',label,'.mat'),'mean_run');

end

function getCorrMD758MARTINbothScans(label)

load(strcat('MD758_ROI_mean_run-',label,'.mat'));

nTotalRuns = length(mean_run);

for iRun=1:nTotalRuns
    
   run = mean_run(iRun).rest;
   run = run';
   rho(iRun).rest_corr = corr(run);
    
end

save(strcat('MD758_ROI_mean_run_corr-',label,'.mat'),'rho');

end

function plotMeanFCMD758MARTINbothScans(label)

load(strcat('MD758_ROI_mean_run_corr-',label,'.mat'));

nROI = 758;
nTotalRuns = length(rho);
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                pairs_rest_rho(iRun,iPair) = single_rho;
            
            end
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z,1);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ mean(pairs_rest_z,1);



mean_corr = zeros(nROI,nROI);
s_corr = zeros(nROI,nROI);

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            %if my_cv_aal(iPair) < cv_criterion & my_cv_aal(iPair) > 0 & p_pairs(iPair) < p_criterion
                
                mean_corr(iROI,iiROI) =  m_pairs_z(iPair);
                s_corr(iROI,iiROI) =  s_pairs_z(iPair);
                
            %end
            
        end
        
    end
    
end

plot_label = strcat('FC-MD758-',label);
plotFCgeneral(mean_corr,s_corr,plot_label);

end


end

