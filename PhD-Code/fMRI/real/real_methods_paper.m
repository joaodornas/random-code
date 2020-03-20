
function real_methods_paper

% plotFCVoxelLevel;

% plotFCVoxelLevelSorted;

% getCorrFromMeanTimeSeriesFrom758ClustersRunByRun;

% getClustersOfVoxelsPer3DCoordinates;

% getCorrFromMeanTimeSeriesFromVoxelsFromClustersOf3DCoordinates;

% getAndPlotVarianceDistributionFromFC758Clusters;

% getAndPlotAllVarianceDistributions;

%%% C400 PARCELLATION

% getCorrFromC400Parcellation;
 
% getCorrFromC400ParcellationInsideParcelVoxelLevel;

% getAndPlotVarianceDistributionFromC400Parcellation;

% getAndPlotVarianceDistributionFromC400ParcellationInsideParcelVoxelLevel;

% getCorrFromC400ParcellationInsideParcelVoxelLevelsavepercluster;

%%% VAN ESSEN PARCELLATION

% getCorrFromVanEssenParcellation;

% getCorrFromVanEssenParcellationInsideParcelVoxelLevel;

% getAndPlotVarianceDistributionFromVanEssenParcellation;

% getAndPlotVarianceDistributionFromVanEssenParcellationInsideParcelVoxelLevel;

% plotMeanFCVanEssen;

%%% INDIVIDUAL VOXELS

% getCorrFromFromVoxelsFromASampleOf20VoxelsPerCluster;

% plotVarianceDistributionVoxelsFromASampleOf20VoxelsPerCluster('Track');
% plotVarianceDistributionVoxelsFromASampleOf20VoxelsPerCluster('Passive');
% plotVarianceDistributionVoxelsFromASampleOf20VoxelsPerCluster('Rest');

% idx_cluster = 329; %%% 45 voxels
% getCorrVoxelsInsideACluster(idx_cluster);
% 
% idx_cluster = 7; %%% 462 voxels
% getCorrVoxelsInsideACluster(idx_cluster);

% idx_cluster = 329; %%% 45 voxels
% plotFCVoxelsInsideACluster(idx_cluster);
%  
% idx_cluster = 7; %%% 462 voxels
% plotFCVoxelsInsideACluster(idx_cluster);
% 
% idx_cluster = 550; %%% 220 voxels
% plotFCVoxelsInsideACluster(idx_cluster);

% for iCluster=71:758
%     plotFCVoxelsInsideACluster(iCluster);
% end

% getCorrVoxelsInsideAClusterPetra;

% getCorrVoxelsInsideALLCluster;

% idx_cluster = 546;
% plotFCVoxelsInsideACluster(idx_cluster);

%%% FUNCTIONAL CLUSTER

% getCorrClustersDifferentClusterSameDifferentAAL;
% getCorrClustersDifferentClusterEverybodyAALPETRA;
% plotVarianceDistributionClustersDifferentClusterSameDifferentAAL;
% plotVarianceDistributionClustersDifferentClusterEverybodyAALPETRA;;

% plotVarianceDistributionClustersDifferentClusterEverybody;

% label = 'Total';
% nStartCluster = 1;
% nEndCluster = 758;
% plotMeanFCFunctionalClusters(nStartCluster,nEndCluster,label);
% 
% label = 'Total';
% nStartCluster = 1;
% nEndCluster = 758;
%plotMeanFCFunctionalClustersPETRA(nStartCluster,nEndCluster,label);
% plotMeanFCFunctionalClustersPETRAcustom(nStartCluster,nEndCluster,label);

% label = 'Total';
% nStartCluster = 1;
% nEndCluster = 758;
% plotMeanFCFunctionalClustersHCP(nStartCluster,nEndCluster,label);

% label = 'Total';
% nStartCluster = 1;
% nEndCluster = 758;
% plotMeanFCFunctionalClustersPETRAWholeSegment(nStartCluster,nEndCluster,label);
 
% label = 'Angular-L';
% nStartCluster = 546;
% nEndCluster = 550;
% plotMeanFCFunctionalClusters(nStartCluster,nEndCluster,label);
% 
% label = 'Frontal-Mid-L';
% nStartCluster = 79;
% nEndCluster = 102;
% plotMeanFCFunctionalClusters(nStartCluster,nEndCluster,label);

%%% ANATOMICAL CLUSTER

% getClustersOfVoxelsPer3DCoordinatesPerROI;
% getCorrAnatomicalDifferentClusterSameDifferentAAL;
% plotVarianceDistributionAnatomicalDifferentClusterSameDifferentAAL;

%%% AAL REGIONS

% getMeanTimeSeriesFrom90AALRunByRun;
% getCorrFromMeanTimeSeriesFrom90AALRunByRun;
% getAndPlotVarianceDistributionFromFC90AAL;
% getAndPlotVarianceDistributionFromFC90AALInsideVoxelLevel;
% plotMeanFCFrom90AAL;

% plotMeanFCFrom90AALPETRA;

% getCorrClustersDifferentClusterEverybodyAALPETRACustomParcellation;

% plotVarCluEveryAALPETRACustomParcellation;

% idx_ROI = 62;
% plotFCSubject(idx_ROI);
% plotFCSubjectPETRA(idx_ROI);
% plotFCSubjectHCP(idx_ROI);

% idx_ROI_A = 62;
% idx_ROI_B = 66;
% getFCcrossROI(idx_ROI_A,idx_ROI_B);

% getInformationRunAALVoxelLevel;
% getInformationRunAALVoxelLevelPETRA;
% getInformationRunAALVoxelLevelHCP;

% idx_ROI = 52;
% checkThreeLevelsOfFC_MD(idx_ROI);
% idx_ROI = 62;
% checkThreeLevelsOfFC_MD(idx_ROI);
% idx_ROI = 66;
% checkThreeLevelsOfFC_MD(idx_ROI);

% checkThreeLevelsOfFC_MD_OnlyDistribution_WholeBrain;

% checkDataWithWithOUTFieldMap;

% idx_ROI = 62;
% plotFCROImeanSorted(idx_ROI);
% idx_ROI = 66;
% plotFCROImeanSorted(idx_ROI);

% idx_ROI_A = 62;
% idx_ROI_B = 66;
% plotFCcrossROImeanSorted(idx_ROI_A,idx_ROI_B);

% idx_ROI_A = 62;
% idx_ROI_B = 66;
% plotFCROImeanSortedcross(idx_ROI_A,idx_ROI_B);

% for idx_ROI=1:90
%        
%     checkThreeLevelsOfFC_MD(idx_ROI);
%     
% end

% label = 'MD';
% RestingState = getRestingStateFromMD;
% getMeanSTDVoxels(RestingState,label);

% label = 'PETRA';
% RestingState = getRestingStateFromPETRA;
% getMeanSTDVoxels(RestingState,label);

% label = 'HCP';
% RestingState = getRestingStateFromHCP;
% getMeanSTDVoxels(RestingState,label);

% label = 'Martin';
% RestingState = getRestingStateFromMARTIN;
% getMeanSTDVoxels(RestingState,label);

% plotSNRRestingStateMDPETRAHCP;

% plotSNRRestingStateMDPETRAHCPv2;

% label = 'MD';
% SNRinaVolume(label);
% 
% label = 'PETRA';
% SNRinaVolume(label);
% 
% label = 'HCP';
% SNRinaVolume(label);
% 
% label = 'Martin';
% SNRinaVolume(label);

% getMeanTimeSeriesFrom90AALRunByRunMARTIN;
% getCorrFromMeanTimeSeriesFrom90AALRunByRunMARTIN;
% getAndPlotVarianceDistributionFromFC90AALMARTIN;

% getMeanMD758MARTIN;
% getCorrMD758MARTIN;
% getAndPlotVarianceMD758MARTIN;

% getAndPlotVarianceDistributionFromFC90AAL;
% plotVarianceDistributionClustersDifferentClusterEverybody;

% plotMeanFCAALMARTIN;
% plotMeanFCMD758MARTIN;

% nStartCluster = 1;
% nEndCluster = 758;
% nStartRun = 17;
% nEndRun = 20;
% label = 'SUBJ-5';
% iSubject = 5;
% plotMeanFCFunctionalClustersONESubject(nStartRun,nEndRun,nStartCluster,nEndCluster,label);
% plotVarianceDistributionClustersDifferentClusterEverybodyONESubject(nStartRun,nEndRun);
% getAndPlotVarianceDistributionFromFC90AALONESubject(nStartRun,nEndRun);
% plotMeanFCFrom90AALONESubject(nStartRun,nEndRun);
% idx_ROI = 53; %%% Occipital_Inf_L
% getFCROImeanONESubject(idx_ROI,iSubject);
% idx_ROI = 61; %%% Parietal_Inf_L
% getFCROImeanONESubject(idx_ROI,iSubject);
% idx_ROI = 53; %%% Occipital_Inf_L
% plotFCROImeanRawONESubject(idx_ROI,iSubject);
% idx_ROI = 61; %%% Parietal_Inf_L
% plotFCROImeanRawONESubject(idx_ROI,iSubject);

% for idx_ROI=1:90
%     getFCROImeanMartin(idx_ROI);
% end

% idx_ROI = 53; %%% Occipital_Inf_L
% plotFCROImeanRawMartin(idx_ROI);
% idx_ROI = 61; %%% Parietal_Inf_L
% plotFCROImeanRawMartin(idx_ROI);

% getMeanMD758SUBJECT5Again;

% [PETRA_custom, PETRA_melodic, MAGDEBURG_custom, MAGDEBURG_melodic] = loadAllDataMartinDataset;

% label = 'PETRA-Custom';
% run = PETRA_custom;
% getMeanMD758MARTINbothScans(run,label);
% getCorrMD758MARTINbothScans(label);
% plotMeanFCMD758MARTINbothScans(label);

% label = 'PETRA-Melodic';
% run = PETRA_melodic;
% getMeanMD758MARTINbothScans(run,label);
% getCorrMD758MARTINbothScans(label);
% plotMeanFCMD758MARTINbothScans(label);

% label = 'MAGDEBURG-Custom';
% run = MAGDEBURG_custom;
% getMeanMD758MARTINbothScans(run,label);
% getCorrMD758MARTINbothScans(label);
% plotMeanFCMD758MARTINbothScans(label);

% label = 'MAGDEBURG-Melodic';
% run = MAGDEBURG_melodic;
% getMeanMD758MARTINbothScans(run,label);
% getCorrMD758MARTINbothScans(label);
% plotMeanFCMD758MARTINbothScans(label);

% label = 'SUBJ5AgainCustom';
% run = loadAllDataSUBJ5Dataset;
% getMeanMD758MARTINbothScans(run,label);
% getCorrMD758MARTINbothScans(label);
% plotMeanFCMD758MARTINbothScans(label);

% label = 'SUBJ5AgainCustomPreMask';
% run = loadAllDataSUBJ5DatasetPreMask;
% getMeanMD758MARTINbothScans(run,label);
% getCorrMD758MARTINbothScans(label);
% plotMeanFCMD758MARTINbothScans(label);

% label = 'SUBJ5Original';
% run = loadAllDataSUBJ5DatasetOriginal;
% getMeanMD758MARTINbothScans(run,label);
% getCorrMD758MARTINbothScans(label);
% plotMeanFCMD758MARTINbothScans(label);

% label = 'SUBJ5OriginalFiltered';
% run = loadAllDataSUBJ5DatasetOriginal;
% getMeanMD758MARTINbothScans(run,label);
% getCorrMD758MARTINbothScans(label);
% plotMeanFCMD758MARTINbothScans(label);

% label = 'ALLLHROriginalTeste2';
% run = getRestingStateFromMD;
% getMeanMD758MARTINbothScans(run,label);
% getCorrMD758MARTINbothScans(label);
% plotMeanFCMD758MARTINbothScans(label);

end


function plotFCVoxelLevel

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-1-22-10-2015\output\FC_Voxels_AAL_ROI\ROI-corr\LHR-SUBJ1-FC-Voxels-AAL-ROI-corr-RestingState-Angular-L.mat');

FC(1).corr = FC_Voxels.run(1).rho_RestingState;
FC(2).corr = FC_Voxels.run(2).rho_RestingState;

clear FC_Voxels

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-8-14-01-2016\output\FC_voxel_AAL_ROI\ROI-corr\LHR-SUBJ8-FC-Voxels-AAL-ROI-corr-RestingState-Angular-L.mat');

FC(3).corr = FC_Voxels.run(4).rho_RestingState;

clear FC_Voxels

for i=1:3
    
    plotCorr(FC(i).corr,int2str(i));
    
end

end

function plotCorr(corr_matrix,label)

    Ncomponent = size(corr_matrix,1);
    
    rlimit = 1;
    
    rmin = -rlimit;
    rmax = rlimit;
    
    f = figure;

    clrmp = colormap('jet');
    %clrmp(32,:) = [1 1 1];

    hold on;
    caxis([rmin rmax]);
    h = imagesc(corr_matrix);
    colormap(clrmp);
    hold off;
    axis 'square'; 

    print(f,'-depsc',strcat('FC_voxel_per_ROI_AAL-',label,'.eps'));

end

function [target_order, cluster_idx_sorted, cluster_size_sorted] = TargetOrder( cluster_assignment, ncluster )

% cluster_assignment assigns nc rows/columns to ncluster clusters

% target_order assigns nc matrix locations to nc rows/columns

% the first locations are filled with rows/columns from the largest cluster

% the next locations are filled with rows/columns from the next smaller
% cluster

% and so on

% cluster_idx_sorted lists the clusters in descending order of size ...

nc = length( cluster_assignment );

kk = [];
kk = find(isnan(cluster_assignment));

if ~isempty(kk)
    
    nNaN = length(kk);
    
    ncluster = ncluster + 1;
    cluster_assignment(kk) = ncluster;
    
end


for kc = 1 : ncluster   % loop over clusters
    
    cluster_size(kc) = length( find( cluster_assignment == kc ) );
    
end

[cluster_size_sorted, cluster_idx_sorted] = sort( cluster_size, 'descend' );   % c_size_sorted = c_size( c_idx_sorted )

cluster_size_sorted;

target_order  = [];  

for kc = 1 : ncluster % loop over sorted clusters
    'target';
    size(target_order);
    'find';
    size(find( cluster_idx_sorted(kc) == cluster_assignment ));
    target_order = [target_order; find( cluster_idx_sorted(kc) == cluster_assignment )];  % append indices of matrix rows/columns in current cluster
       
end

target_order = reshape( target_order, 1, nc );

end

function [m_sorted, source_order] = SortMatrix( m_original, source_order, target_order, nc )

% typically the rows/columns of a matrix are numbered by location, ie., 1st, 2nd, 3rd, ...
% however, when we want to rearrange a matrix, we need to distinguish between
% row/column number and location in the matrix

% source_order           lists row/column numbers in the source order (1st, 2nd, 3rd, ... location)
% target_order            lists row/column numbers in the target order (1st, 2nd, 3rd, ... location)

m_sorted = m_original;

ix_must_swap = find( source_order ~= target_order );    % find positions where component numbers differ in source and target

ix = ix_must_swap(1);     % first location to be rectified

Nswap = length( ix_must_swap );

for ic = 1 : Nswap                  % loop number of necessary swaps
    
    i_from   = ix;                  % location of component to be moved away
    c_source = source_order(ix);    % number of component to be moved away
 
    c_target = target_order(ix);                 % number of component to replace it
    i_to     = find( c_target == source_order );   % location of other component to replace it


    
    if i_from ~= i_to               % locations are different are different
        
        [m_sorted, source_order] = MSwap( m_sorted, source_order, i_from, i_to ); % swap components
                                            
                                            % now location i_from has component c_to, which is correct
                                            % but location i_to has component c_from, which is incorrect
                                            
        ix = i_to;                          % next location to be rectified                                
        
              
    elseif ic < Nswap
        ix_must_swap = find( source_order ~= target_order ); 
        ix = ix_must_swap(1);
    end
        
end

end

function [m_out, order_out]  = MSwap( m_in, order_in, i, k )

temp = m_in;

temp(i,:) = m_in(k,:);

temp(k,:) = m_in(i,:);

m_out = temp;

m_out(:,i) = temp(:,k);

m_out(:,k) = temp(:,i);

order_out = order_in;

order_out(i) = order_in(k);

order_out(k) = order_in(i);

end

function plotFCVoxelLevelSorted

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-1-22-10-2015\output\FC_Voxels_AAL_ROI\ROI-corr\LHR-SUBJ1-FC-Voxels-AAL-ROI-corr-RestingState-Angular-L.mat');

FC(1).corr = FC_Voxels.run(1).rho_RestingState;
FC(2).corr = FC_Voxels.run(2).rho_RestingState;

clear FC_Voxels

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-8-14-01-2016\output\FC_voxel_AAL_ROI\ROI-corr\LHR-SUBJ8-FC-Voxels-AAL-ROI-corr-RestingState-Angular-L.mat');

FC(3).corr = FC_Voxels.run(4).rho_RestingState;

clear FC_Voxels

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\KMeans\LHR-All-Subjects-FC-Voxels-AAL-ROI-corr-RestingState-Angular-L-KMeans.mat');

cluster_assignment = FC_KMeans.all_runs.IdxClusters;

Ncluster = FC_KMeans.all_runs.Ncluster;
 
Ncomponent = size(FC(1).corr,1);
    
source_order = 1:Ncomponent;
[target_order, cluster_idx_sorted] = TargetOrder( cluster_assignment, Ncluster );

for i=1:3
 
    [rho_sorted, ~] = SortMatrix( FC(i).corr, source_order, target_order, Ncomponent );
    
    plotCorr(rho_sorted,strcat(int2str(i),'-sorted'));
    
end

end

function getClustersOfVoxelsPer3DCoordinates

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nClusters = 758;
nROI = 90;
nRuns = 32;
MNI_img = [91 109 91];
nTotalVoxels = 160990;
nTotalClusters = 758;
nClusterSize = 200;

iiVoxel = 0;

idx = zeros(3,nTotalVoxels);
all_idx_voxels = [];

for iROI=1:nROI

    nClusters = length(ROI(iROI).clusters);
    
    for iCluster=1:nClusters
       
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
        
        all_idx_voxels = [all_idx_voxels;idx_voxels];
        
        nVoxels = length(idx_voxels);
        
        for iVoxel=1:nVoxels
            
           iiVoxel = iiVoxel + 1;
            
           [idx(1,iiVoxel) idx(2,iiVoxel) idx(3,iiVoxel)] = ind2sub(MNI_img,idx_voxels(iVoxel)); 
            
        end
        
        
    end
    
end

D = pdist(idx','euclidean');

closest_voxels = zeros(nTotalVoxels,nClusterSize*5);

for iVoxel=1:nTotalVoxels
    
    %D((I-1)*(M-I/2)+J-I);
    
    distances = zeros(1,nTotalVoxels);
    
    if mod(iVoxel,20000) == 0; disp(int2str(iVoxel)); end

    distances((iVoxel+1):nTotalVoxels) = D((iVoxel-1)*(nTotalVoxels-iVoxel/2) + ((iVoxel+1):nTotalVoxels)-iVoxel);

    if iVoxel > 1

        distances((iVoxel-1):-1:1) = D((((iVoxel-1):-1:1)-1).*(nTotalVoxels-((iVoxel-1):-1:1)./2)+iVoxel-((iVoxel-1):-1:1));
    
    end

    distances(iVoxel) = 0;

    [s,i] = sort(distances,'ascend');
    
    closest_voxels(iVoxel,:) = i(1:nClusterSize*5);
    
end

idx_clusters = kmeans(closest_voxels,nTotalClusters);

for iCluster=1:nTotalClusters
   
    idx_voxels_on_cluster = find(idx_clusters==iCluster);
    
    clusters(iCluster).idx_voxels = all_idx_voxels(idx_voxels_on_cluster);
    
end

save('Clusters-Voxels-Per-3D-Coordinate.mat','clusters','idx');

end

function getCorrFromMeanTimeSeriesFromVoxelsFromClustersOf3DCoordinates

load('Clusters-Voxels-Per-3D-Coordinate.mat');

nTR = 150;
nRuns = 4;

all_settings = getAllSettings;

nSubjects = length(all_settings);

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

iiRun = 0;

for iSubject=1:nSubjects

    settings = all_settings(iSubject).settings;

    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;

    kind = 'RestingState';
    for irun=1:nRuns;
        [RestingState(irun).run, RestingState(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end

    for irun=1:nRuns

        iiRun = iiRun + 1;
        
        disp(strcat('Run:',int2str(iiRun)));
        
        sample_clusters = zeros(nTR,length(clusters));

        for iCluster=1:length(clusters)

            idx_voxels = clusters(iCluster).idx_voxels;

            sample_voxels = zeros(nTR,length(idx_voxels));

            for iVoxel=1:length(idx_voxels)

                [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

                sample_voxels(:,iVoxel) = RestingState(irun).run(idxx,idxy,idxz,:);

            end
            
            sample_clusters(:,iCluster) = mean(sample_voxels,1);

        end
        
        coordinates_clusters_corr(iiRun,:) = nonzeros(triu(corr(sample_clusters)));

    end

end

save('Coordinates_Clusters_Corr.mat','coordinates_clusters_corr');

end

function getAndPlotVarianceDistributionFromFC758Clusters

load('758_Clusters_mean_run_corr.mat');

nClusters = 758;
iPair = 0;
nRuns = 32;

pairs = zeros(nClusters*(nClusters-1),nRuns);

for iClusters=1:nClusters
    
    for iiClusters=1:nClusters
        
        if iClusters ~= iiClusters
        
            iPair = iPair + 1;
            
            if mod(iPair,20000) == 0; disp(int2str(iPair)); end

            for iRun=1:nRuns

                pairs(iRun,iPair) = corr_clusters(iRun).corr(iClusters,iiClusters);

            end
        
        end
        
    end
    
end

my_var = var(pairs);

my_cv = std(pairs) ./ abs(mean(pairs,1));

min_val = min(my_var);
max_val = max(my_var);
vector_one = my_var;
label_one = 'Variance';
title_label = '758-Clusters-Variance';
x_label = 'variance';
y_label = 'density';

plotDistributionDensityOnlyOne(min_val,max_val,vector_one,label_one,title_label,x_label,y_label);

min_val = min(my_cv);
max_val = max(my_cv);
vector_one = my_cv;
label_one = 'Coefficient-Of-Variation';
title_label = '758-Clusters-Coefficient-Of-Variation';
x_label = 'CV';
y_label = 'density';

plotDistributionDensityOnlyOne(min_val,max_val,vector_one,label_one,title_label,x_label,y_label);

plotDistributionDensityOnlyOneXInLog(min_val,max_val,vector_one,label_one,title_label,x_label,y_label);


end

function getAndPlotAllVarianceDistributions

%%% ROI LEVEL

load('AAL_ROI_mean_run_corr.mat');

nROI = 90;
nTotalRuns = 32;

iPair = 0;

pairs_aal = zeros(nTotalRuns,nROI*(nROI-1));

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).run(iROI,iiROI);
                
                pairs_aal(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            end
            
        end
        
    end
    
end

my_var_aal = var(pairs_aal);

my_cv_aal = std(pairs_aal) ./ abs(mean(pairs_aal,1));

my_mean_aal = mean(pairs_aal,1);
my_std_aal = std(pairs_aal);

%%% CLUSTER LEVEL

load('758_Clusters_mean_run_corr.mat');

nClusters = 758;
iPair = 0;
nRuns = 32;

pairs_clusters = zeros(nRuns,nClusters*(nClusters-1));

for iClusters=1:nClusters
    
    for iiClusters=1:nClusters
        
        if iClusters ~= iiClusters
        
            iPair = iPair + 1;
            
            if mod(iPair,20000) == 0; disp(int2str(iPair)); end

            for iRun=1:nRuns

                single_rho = corr_clusters(iRun).corr(iClusters,iiClusters);
                
                pairs_clusters(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

            end
        
        end
        
    end
    
end

my_var_clusters = var(pairs_clusters);

my_cv_clusters = std(pairs_clusters) ./ abs(mean(pairs_clusters,1));

my_mean_clusters = mean(pairs_clusters,1);
my_std_clusters = std(pairs_clusters);

%%% VOXEL LEVEL

load('Corr-Samples-Of-Voxels.mat');

% voxels_corr = getCorrFromFromVoxelsFromASampleOf20VoxelsPerCluster;

nVoxels = 15160;
nRuns = 32; 
           
pairs_voxels = (1/2)*log((1 + voxels_corr)./(1 - voxels_corr));
            
my_var_voxels = var(pairs_voxels);

my_cv_voxels = std(pairs_voxels) ./ abs(mean(pairs_voxels,1));

my_mean_voxels = mean(pairs_voxels,1);
my_std_voxels = std(pairs_voxels);

vector_one = my_std_aal;
vector_two = my_mean_aal;
vector_three = my_cv_aal;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Contour - AAL';

MyContourCorrelation(vector_two,vector_one,title_label);

vector_one = my_std_clusters;
vector_two = my_mean_clusters;
vector_three = my_cv_clusters;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Contour - Clusters';

MyContourCorrelation(vector_two,vector_one,title_label);

vector_one = my_std_voxels;
vector_two = my_mean_voxels;
vector_three = my_cv_voxels;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Contour - Voxels';

MyContourCorrelation(vector_two,vector_one,title_label);
 

end

function getFCcrossROI(idx_ROI_A,idx_ROI_B)

nTotalRuns = 32;
nRuns = 4;
nTR = 150;
MNI_img = [91 109 91];

all_settings = getAllSettings;
nSubjects = length(all_settings);

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

iiRun = 0;
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
    
    kind = 'RestingState';
    for irun=1:nRuns;
        iiRun = iiRun + 1;
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
end

ROI_A = AAL_ROI(idx_ROI_A).ID;
ROI_A_Label = AAL_ROI(idx_ROI_A).Nom_L;
ROI_A_Label = strrep(ROI_A_Label,'_','-');
idx_voxels_A = find(AAL_img == ROI_A);
nVoxels_A = length(idx_voxels_A);
            
ROI_B = AAL_ROI(idx_ROI_B).ID;
ROI_B_Label = AAL_ROI(idx_ROI_B).Nom_L;
ROI_B_Label = strrep(ROI_B_Label,'_','-');
idx_voxels_B = find(AAL_img == ROI_B);
nVoxels_B = length(idx_voxels_B);
            
all_voxels = [idx_voxels_A;idx_voxels_B];

for iiRun=1:nTotalRuns
    
    disp(strcat('run:',int2str(iiRun)));
   
    nVoxels = length(all_voxels);
    
    voxels_ts = zeros(nTR,nVoxels);
    
    for iVoxel=1:nVoxels
       
        [idxx,idxy,idxz] = ind2sub(MNI_img,all_voxels(iVoxel));
        
        voxels_ts(:,iVoxel) = squeeze(RestingState(iiRun).run(idxx,idxy,idxz,:));
        
    end
    
    disp(datestr(now));
    voxels_corr = corr(voxels_ts);
    
    save(strcat('LHR-FC-Voxel-',ROI_A_Label,'-',ROI_B_Label,'-','run','-',int2str(iiRun)),'idx_voxels_A','idx_voxels_B','idx_ROI_A','idx_ROI_B','ROI_A_Label','ROI_B_Label','nVoxels_A','nVoxels_B','voxels_corr');
    
end

end

function getFCcrossROImean(idx_ROI_A,idx_ROI_B)

nTotalRuns = 32;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ROI_A = AAL_ROI(idx_ROI_A).ID;
ROI_A_Label = AAL_ROI(idx_ROI_A).Nom_L;
ROI_A_Label = strrep(ROI_A_Label,'_','-');
idx_voxels_A = find(AAL_img == ROI_A);
nVoxels_A = length(idx_voxels_A);
            
ROI_B = AAL_ROI(idx_ROI_B).ID;
ROI_B_Label = AAL_ROI(idx_ROI_B).Nom_L;
ROI_B_Label = strrep(ROI_B_Label,'_','-');
idx_voxels_B = find(AAL_img == ROI_B);
nVoxels_B = length(idx_voxels_B);

for iRun=1:nTotalRuns
    
   load(strcat('LHR-FC-Voxel-',ROI_A_Label,'-',ROI_B_Label,'-','run','-',int2str(iRun),'.mat'));
   
   FC(iRun,:,:) = voxels_corr(:,:);
   
   clear voxels_corr
   
end

m_FC = mean(FC,1);

save(strcat('LHR-FC-Voxel-',ROI_A_Label,'-',ROI_B_Label,'-','Mean','.mat'),'m_FC');

end

function getFCROImean(idx_ROI)

nTotalRuns = 32;
nSubjects = 8;
nRuns = 4;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ROI = AAL_ROI(idx_ROI).ID;
ROI_Label = AAL_ROI(idx_ROI).Nom_L;
ROI_Label = strrep(ROI_Label,'_','-');
idx_voxels = find(AAL_img == ROI);
nVoxels = length(idx_voxels);

iiRun = 0;
for iSubject=1:nSubjects
    
   load(strcat('LHR-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-RestingState-',ROI_Label,'.mat'));
   
   for iRun=1:nRuns
        
        iiRun = iiRun + 1;
        
        FC(iiRun,:,:) = FC_Voxels.run(iRun).rho_RestingState(:,:);
   
   end
   
   clear FC_Voxels
   
end

m_FC = squeeze(mean(FC,1));

save(strcat('LHR-FC-Voxel-',ROI_Label,'-','Mean','.mat'),'m_FC');

end

function getFCROImeanMartin(idx_ROI)

nTotalRuns = 4;
nRuns = 4;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ROI = AAL_ROI(idx_ROI).ID;
ROI_Label = AAL_ROI(idx_ROI).Nom_L;
ROI_Label = strrep(ROI_Label,'_','-');
idx_voxels = find(AAL_img == ROI);
nVoxels = length(idx_voxels);

for iRun=1:nRuns
    
   load(strcat('Martin-SUBJ1-FC-Voxels-AAL-ROI-corr-Martin-RestingState-',ROI_Label,'-',int2str(iRun),'.mat'));

   FC(iRun,:,:) = FC_Voxels.run(iRun).rho_RestingState(:,:);

   clear FC_Voxels
   
end

m_FC = squeeze(mean(FC,1));

save(strcat('Martin-SUBJ1-FC-Voxel-',ROI_Label,'-','Mean','.mat'),'m_FC');

end

function getFCROImeanONESubject(idx_ROI,iSubject)

nTotalRuns = 32;
nRuns = 4;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ROI = AAL_ROI(idx_ROI).ID;
ROI_Label = AAL_ROI(idx_ROI).Nom_L;
ROI_Label = strrep(ROI_Label,'_','-');
idx_voxels = find(AAL_img == ROI);
nVoxels = length(idx_voxels);

iiRun = 0;
    
load(strcat('LHR-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-RestingState-',ROI_Label,'.mat'));

for iRun=1:nRuns

    iiRun = iiRun + 1;

    FC(iiRun,:,:) = FC_Voxels.run(iRun).rho_RestingState(:,:);

end

clear FC_Voxels

m_FC = squeeze(mean(FC,1));

save(strcat('LHR-FC-Voxel-',ROI_Label,'-','Mean','-','SUBJ-',int2str(iSubject),'.mat'),'m_FC');

end

function clusterFCROImean(idx_ROI)

nTotalRuns = 32;
nSubjects = 8;
nRuns = 4;
nVoxelsPerCluster = 200;

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ROI = AAL_ROI(idx_ROI).ID;
ROI_Label = AAL_ROI(idx_ROI).Nom_L;
ROI_Label = strrep(ROI_Label,'_','-');
idx_voxels = find(AAL_img == ROI);
nVoxels = length(idx_voxels);

load(strcat('LHR-FC-Voxel-',ROI_Label,'-','Mean','.mat'));

% m_FC_z(:,:) = (1/2).*log((1 + m_FC)./(1 - m_FC));
% 
% pval = abs(m_FC_z) < z_threshold; 

nClusters = floor( size(m_FC,1) / nVoxelsPerCluster );

[IdxClusters, Tidx, Ncluster] = ClusterWithKmeans_all( m_FC, nClusters );
        
save(strcat('LHR-FC-Voxel-',ROI_Label,'-','Mean-KMeans','.mat'),'IdxClusters', 'Tidx', 'Ncluster');

end

function plotFCcrossROImeanSorted(idx_ROI_A,idx_ROI_B)

nTotalRuns = 32;
nSubjects = 8;
nRuns = 4;
nVoxelsPerCluster = 200;

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ROI_A = AAL_ROI(idx_ROI_A).ID;
ROI_Label_A = AAL_ROI(idx_ROI_A).Nom_L;
ROI_Label_A = strrep(ROI_Label_A,'_','-');
idx_voxels_A = find(AAL_img == ROI_A);
nVoxels_A = length(idx_voxels_A);

ROI_B = AAL_ROI(idx_ROI_B).ID;
ROI_Label_B = AAL_ROI(idx_ROI_B).Nom_L;
ROI_Label_B = strrep(ROI_Label_B,'_','-');
idx_voxels_B = find(AAL_img == ROI_B);
nVoxels_B = length(idx_voxels_B);

load(strcat('LHR-FC-Voxel-',ROI_Label_A,'-',ROI_Label_B,'-Mean','.mat'));

m_FC_A = m_FC(1:nVoxels_A,(nVoxels_A+1):(nVoxels_A+nVoxels_B));

m_FC_B = m_FC((nVoxels_A+1):(nVoxels_A+nVoxels_B),1:nVoxels_A);

load(strcat('LHR-FC-Voxel-',ROI_Label_A,'-','Mean-KMeans','.mat'));

IdxClusters_A = IdxClusters; 
Tidx_A = Tidx; 
Ncluster_A = Ncluster;

load(strcat('LHR-FC-Voxel-',ROI_Label_B,'-','Mean-KMeans','.mat'));

IdxClusters_B = IdxClusters; 
Tidx_B = Tidx; 
Ncluster_B = Ncluster;




end

function plotFCROImeanSorted(idx_ROI)

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ROI = AAL_ROI(idx_ROI).ID;
ROI_Label = AAL_ROI(idx_ROI).Nom_L;
ROI_Label = strrep(ROI_Label,'_','-');
idx_voxels = find(AAL_img == ROI);
nVoxels = length(idx_voxels);

load(strcat('LHR-FC-Voxel-',ROI_Label,'-','Mean','.mat'));

load(strcat('LHR-FC-Voxel-',ROI_Label,'-','Mean-KMeans','.mat'));

cluster_assignment = IdxClusters;

% Ncluster = Ncluster;
 
Ncomponent = size(m_FC,1);
    
source_order = 1:Ncomponent;
[target_order, cluster_idx_sorted] = TargetOrder( cluster_assignment, Ncluster );

[rho_sorted, source_order] = SortMatrix( m_FC, source_order, target_order, Ncomponent );

rho_sorted_z(:,:) = (1/2).*log((1 + rho_sorted)./(1 - rho_sorted));

rho_sorted_z(isinf(rho_sorted_z)) = NaN;
    
plotFCgeneral2(rho_sorted_z,strcat('LHR-FC-Voxel-',ROI_Label,'-','Mean-Sorted'));

rho_sorted_z_m = zeros(Ncluster,Ncluster);

clusters = unique(cluster_idx_sorted,'stable');

for iClu=1:Ncluster
    
    nVoxels(iClu) = length(find(cluster_assignment==clusters(iClu)));
    
end

voxels_so_far_i = 0;
for iClu=1:Ncluster
    
    nVoxels_i = nVoxels(iClu);
    
    voxels_so_far_i = voxels_so_far_i + nVoxels_i;
    
    voxels_so_far_ii = 0;
    for iiClu=1:Ncluster
        
        nVoxels_ii = nVoxels(iiClu);
        
        voxels_so_far_ii = voxels_so_far_ii + nVoxels_ii;
        
        if iClu == 1; start_l = 1; else start_l = voxels_so_far_i - nVoxels_i + 1; end
        if iiClu == 1; start_c = 1; else start_c = voxels_so_far_ii - nVoxels_ii + 1; end
        
        end_l = voxels_so_far_i;
        end_c = voxels_so_far_ii;
        
        tmp_corr = rho_sorted_z(start_l:end_l,start_c:end_c);
        
        rho_sorted_z_m(iClu,iiClu) = nanmean(tmp_corr(:));
     
    end
    
end

plotFCgeneral2(rho_sorted_z_m,strcat('LHR-FC-Voxel-',ROI_Label,'-','Mean-Sorted-clusters'));

save(strcat('LHR-FC-Voxel-',ROI_Label,'-','cluster_idx_sorted'),'cluster_idx_sorted');

end

function plotFCROImeanSortedcross(idx_ROI_A,idx_ROI_B)

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ROI_A = AAL_ROI(idx_ROI_A).ID;
ROI_Label_A = AAL_ROI(idx_ROI_A).Nom_L;
ROI_Label_A = strrep(ROI_Label_A,'_','-');
idx_voxels_A = find(AAL_img == ROI_A);
all_nVoxels_A = length(idx_voxels_A);

ROI_B = AAL_ROI(idx_ROI_B).ID;
ROI_Label_B = AAL_ROI(idx_ROI_B).Nom_L;
ROI_Label_B = strrep(ROI_Label_B,'_','-');
idx_voxels_B = find(AAL_img == ROI_B);
all_nVoxels_B = length(idx_voxels_B);

load(strcat('LHR-FC-Voxel-',ROI_Label_A,'-','Mean-KMeans','.mat'));
cluster_assignment_A = IdxClusters;
load(strcat('LHR-FC-Voxel-',ROI_Label_B,'-','Mean-KMeans','.mat'));
cluster_assignment_B = IdxClusters;

load(strcat('LHR-FC-Voxel-',ROI_Label_A,'-','cluster_idx_sorted'));
cluster_idx_sorted_A = cluster_idx_sorted;
load(strcat('LHR-FC-Voxel-',ROI_Label_B,'-','cluster_idx_sorted'));
cluster_idx_sorted_B = cluster_idx_sorted;

load(strcat('LHR-FC-Voxel-',ROI_Label_A,'-',ROI_Label_B,'-','Mean','.mat'));

Ncomponent_A = nVoxels_A;
Ncomponent_B = nVoxels_B;

clusters_A = unique(cluster_assignment_A);
clusters_B = unique(cluster_assignment_B);

clusters_A(isnan(clusters_A)) = [];
clusters_B(isnan(clusters_B)) = [];

nClu_A = length(clusters_A);
nClu_B = length(clusters_B);

for iClu=1:nClu_A
    
    nVoxels_A(iClu) = length(find(cluster_assignment_A==iClu));
    
end

for iClu=1:nClu_B
    
    nVoxels_B(iClu) = length(find(cluster_assignment_B==iClu));
    
end

all_voxels_A = [];
for iClu=1:nClu_A
    
    this_cluster = cluster_idx_sorted_A(iClu);
    
    idx_voxels = find(cluster_assignment_A==this_cluster);
    
    all_voxels_A = [all_voxels_A;idx_voxels];
    
end

all_voxels_B = [];
for iClu=1:nClu_B
    
    this_cluster = cluster_idx_sorted_B(iClu);
    
    idx_voxels = find(cluster_assignment_B==this_cluster);
    
    all_voxels_B = [all_voxels_B;idx_voxels];
    
end

m_FC_A = m_FC(1:all_nVoxels_A,(all_nVoxels_A+1):(all_nVoxels_A+all_nVoxels_B));
m_FC_B = m_FC((all_nVoxels_A+1):(all_nVoxels_A+all_nVoxels_B),1:all_nVoxels_A);

m_FC_A = m_FC_A(all_voxels_A(:),:);
m_FC_A = m_FC_A(:,all_voxels_B(:));

m_FC_B = m_FC_B(all_voxels_B(:),:);
m_FC_B = m_FC_B(:,all_voxels_A(:));

rho_sorted_z_A(:,:) = (1/2).*log((1 + m_FC_A)./(1 - m_FC_A));
rho_sorted_z_B(:,:) = (1/2).*log((1 + m_FC_B)./(1 - m_FC_B));

rho_sorted_z_A(isinf(rho_sorted_z_A)) = NaN;
rho_sorted_z_B(isinf(rho_sorted_z_B)) = NaN;
    
plotFCgeneral2(rho_sorted_z_A,strcat('LHR-FC-Voxel-',ROI_Label_A,'-',ROI_Label_B,'-','Mean-Sorted'));
plotFCgeneral2(rho_sorted_z_B,strcat('LHR-FC-Voxel-',ROI_Label_B,'-',ROI_Label_A,'-','Mean-Sorted'));

rho_sorted_z_m_A = zeros(Ncluster,Ncluster);
rho_sorted_z_m_B = zeros(Ncluster,Ncluster);

voxels_so_far_i = 0;
for iClu=1:nClu_A
    
    nVoxels_i = nVoxels_A(iClu);
    
    voxels_so_far_i = voxels_so_far_i + nVoxels_i;
    
    voxels_so_far_ii = 0;
    for iiClu=1:nClu_B
        
        nVoxels_ii = nVoxels_B(iiClu);
        
        voxels_so_far_ii = voxels_so_far_ii + nVoxels_ii;
        
        if iClu == 1; start_l = 1; else start_l = voxels_so_far_i - nVoxels_i + 1; end
        if iiClu == 1; start_c = 1; else start_c = voxels_so_far_ii - nVoxels_ii + 1; end
        
        end_l = voxels_so_far_i;
        end_c = voxels_so_far_ii;
        
        tmp_corr = rho_sorted_z_A(start_l:end_l,start_c:end_c);
        
        rho_sorted_z_m_A(iClu,iiClu) = nanmean(tmp_corr(:));
     
    end
    
end

plotFCgeneral2(rho_sorted_z_m_A,strcat('LHR-FC-Voxel-',ROI_Label_A,'-',ROI_Label_B,'-''Mean-Sorted-clusters'));

voxels_so_far_i = 0;
for iClu=1:nClu_B
    
    nVoxels_i = nVoxels_B(iClu);
    
    voxels_so_far_i = voxels_so_far_i + nVoxels_i;
    
    voxels_so_far_ii = 0;
    for iiClu=1:nClu_A
        
        nVoxels_ii = nVoxels_A(iiClu);
        
        voxels_so_far_ii = voxels_so_far_ii + nVoxels_ii;
        
        if iClu == 1; start_l = 1; else start_l = voxels_so_far_i - nVoxels_i + 1; end
        if iiClu == 1; start_c = 1; else start_c = voxels_so_far_ii - nVoxels_ii + 1; end
        
        end_l = voxels_so_far_i;
        end_c = voxels_so_far_ii;
        
        tmp_corr = rho_sorted_z_B(start_l:end_l,start_c:end_c);
        
        rho_sorted_z_m_B(iClu,iiClu) = nanmean(tmp_corr(:));
     
    end
    
end

plotFCgeneral2(rho_sorted_z_m_B,strcat('LHR-FC-Voxel-',ROI_Label_B,'-',ROI_Label_A,'-','Mean-Sorted-clusters'));

end

function plotFCROImeanRawONESubject(idx_ROI,iSubject)

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ROI = AAL_ROI(idx_ROI).ID;
ROI_Label = AAL_ROI(idx_ROI).Nom_L;
ROI_Label = strrep(ROI_Label,'_','-');
idx_voxels = find(AAL_img == ROI);
nVoxels = length(idx_voxels);

load(strcat('LHR-FC-Voxel-',ROI_Label,'-','Mean','-','SUBJ-',int2str(iSubject),'.mat'));

m_FC_z(:,:) = (1/2).*log((1 + m_FC)./(1 - m_FC));

m_FC_z(isinf(m_FC_z)) = NaN;
    
plotFCgeneral2(m_FC_z,strcat('LHR-FC-Voxel-',ROI_Label,'-','Mean-SUBJ-',int2str(iSubject)));

end

function plotFCROImeanRawMartin(idx_ROI)

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ROI = AAL_ROI(idx_ROI).ID;
ROI_Label = AAL_ROI(idx_ROI).Nom_L;
ROI_Label = strrep(ROI_Label,'_','-');
idx_voxels = find(AAL_img == ROI);
nVoxels = length(idx_voxels);

load(strcat('Martin-SUBJ1-FC-Voxel-',ROI_Label,'-','Mean','.mat'));

m_FC_z(:,:) = (1/2).*log((1 + m_FC)./(1 - m_FC));

m_FC_z(isinf(m_FC_z)) = NaN;
    
plotFCgeneral2(m_FC_z,strcat('Martin-SUBJ1-FC-Voxel-',ROI_Label,'-','Mean'));

end

function checkThreeLevelsOfFC_MD(idx_ROI)

tic

nTotalRuns_MD = 32;
nTotalRuns_PETRA = 27;
nTotalRuns_HCP = 32;

nRuns_MD = 4;
nRuns_PETRA = 3;
nRuns_HCP = 4;

nSubjects_MD = 8;
nSubjects_PETRA = 9;
nSubjects_HCP = 8;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

ROI = AAL_ROI(idx_ROI).ID;
ROI_Label = AAL_ROI(idx_ROI).Nom_L;
ROI_Label = strrep(ROI_Label,'_','-');
idx_voxels = find(AAL_img == ROI);
nVoxels = length(idx_voxels);

disp(strcat('AAL:',int2str(idx_ROI),':',ROI_Label));

disp('Magdeburg');

iiRun = 0;
for iSubject=1:nSubjects_MD
    
    load(strcat('LHR-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-RestingState-',ROI_Label,'.mat'));
    
    for iRun=1:nRuns_MD
        
        iiRun = iiRun + 1;
        
        disp(strcat('Run:',int2str(iiRun)));
        
        FC_MD(iiRun,:,:) = FC_Voxels.run(iRun).rho_RestingState;
        
        sub_FC_MD(iRun,:,:) = FC_Voxels.run(iRun).rho_RestingState;
        
        if iRun==1
            
            this_corr(:,:) = squeeze(FC_MD(iiRun,:,:));
            
            this_corr_z(:,:) = (1/2).*log((1 + this_corr)./(1 - this_corr));
        
            %plotFCgeneral2(this_corr_z,strcat('LHR-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-RestingState-',ROI_Label,'-run-',int2str(iRun)));
            
        end
        
    end
    
    sub_FC_z_MD(:,:,:) = (1/2).*log((1 + sub_FC_MD)./(1 - sub_FC_MD));
    
    m_sub_FC_z = squeeze(nanmean(sub_FC_z_MD,1));
    s_sub_FC_z = squeeze(nanstd(sub_FC_z_MD));
    
    m_sub_FC_z_all_MD(iSubject,:,:) = m_sub_FC_z(:,:);
    
    %plotFCgeneral(m_sub_FC_z,s_sub_FC_z,strcat('LHR-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-RestingState-',ROI_Label,'-','mean'));
    
    clear FC_Voxels
            
end

FC_z_MD(:,:,:) = (1/2).*log((1 + FC_MD)./(1 - FC_MD));

m_FC_z_MD = squeeze(nanmean(FC_z_MD,1));
s_FC_z_MD = squeeze(nanstd(FC_z_MD));

%plotFCgeneral(m_FC_z_MD,s_FC_z_MD,strcat('LHR-Average','-FC-Voxels-AAL-ROI-corr-RestingState-',ROI_Label));

clear this_corr this_corr_z m_sub_FC_z s_sub_FC_z

disp('PETRA');

all_settings = getAllSettingsPetra;

iiRun = 0;
for iSubject=1:nSubjects_PETRA
    
    settings = all_settings(iSubject).settings;
    
    load(strcat('PETRA-',settings.codes.subject,'-FC-Voxels-AAL-ROI-corr-Petra-RestingState-',ROI_Label,'.mat'));
    
    for iRun=1:nRuns_PETRA
        
        iiRun = iiRun + 1;
        
        disp(strcat('Run:',int2str(iiRun)));
        
        FC_PETRA(iiRun,:,:) = FC_Voxels.run(iRun).rho_RestingState;
        
        sub_FC_PETRA(iRun,:,:) = FC_Voxels.run(iRun).rho_RestingState;
        
        if iRun==1
            
            this_corr(:,:) = squeeze(FC_PETRA(iiRun,:,:));
            
            this_corr_z(:,:) = (1/2).*log((1 + this_corr)./(1 - this_corr));
        
            %plotFCgeneral2(this_corr_z,strcat('PETRA-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-Petra-RestingState-',ROI_Label,'-run-',int2str(iRun)));
            
        end
        
    end
    
    sub_FC_z_PETRA(:,:,:) = (1/2).*log((1 + sub_FC_PETRA)./(1 - sub_FC_PETRA));
    
    m_sub_FC_z = squeeze(nanmean(sub_FC_z_PETRA,1));
    s_sub_FC_z = squeeze(nanstd(sub_FC_z_PETRA));
    
    m_sub_FC_z_all_PETRA(iSubject,:,:) = m_sub_FC_z(:,:);
    
    %plotFCgeneral(m_sub_FC_z,s_sub_FC_z,strcat('PETRA-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-Petra-RestingState-',ROI_Label,'-','mean'));
    
    clear FC_Voxels
            
end

FC_z_PETRA(:,:,:) = (1/2).*log((1 + FC_PETRA)./(1 - FC_PETRA));

m_FC_z_PETRA = squeeze(nanmean(FC_z_PETRA,1));
s_FC_z_PETRA = squeeze(nanstd(FC_z_PETRA));

%plotFCgeneral(m_FC_z_PETRA,s_FC_z_PETRA,strcat('PETRA-Average','-FC-Voxels-AAL-ROI-corr-Petra-RestingState-',ROI_Label));

clear this_corr this_corr_z m_sub_FC_z s_sub_FC_z

disp('HCP');

all_settings = getAllSettingsHCP;

iiRun = 0;
for iSubject=1:nSubjects_HCP
    
    settings = all_settings(iSubject).settings;
    
    for iRun=1:nRuns_HCP
        
        load(strcat('HCP-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-HCP-RestingState-',ROI_Label,'-',int2str(iRun),'.mat'));
        
        iiRun = iiRun + 1;
        
        disp(strcat('Run:',int2str(iiRun)));
        
        FC_HCP(iiRun,:,:) = FC_Voxels.run(iRun).rho_RestingState;
        
        sub_FC_HCP(iRun,:,:) = FC_Voxels.run(iRun).rho_RestingState;
        
        if iRun==1
            
            this_corr(:,:) = squeeze(FC_HCP(iiRun,:,:));
            
            this_corr_z(:,:) = (1/2).*log((1 + this_corr)./(1 - this_corr));
        
            %plotFCgeneral2(this_corr_z,strcat('HCP-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-HCP-RestingState-',ROI_Label,'-run-',int2str(iRun)));
            
        end
        
    end
    
    sub_FC_z_HCP(:,:,:) = (1/2).*log((1 + sub_FC_HCP)./(1 - sub_FC_HCP));
    
    m_sub_FC_z = squeeze(nanmean(sub_FC_z_HCP,1));
    s_sub_FC_z = squeeze(nanstd(sub_FC_z_HCP));
    
    m_sub_FC_z_all_HCP(iSubject,:,:) = m_sub_FC_z(:,:);
    
    %plotFCgeneral(m_sub_FC_z,s_sub_FC_z,strcat('HCP-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-HCP-RestingState-',ROI_Label,'-','mean'));
    
    clear FC_Voxels
            
end

FC_z_HCP(:,:,:) = (1/2).*log((1 + FC_HCP)./(1 - FC_HCP));

m_FC_z_HCP = squeeze(nanmean(FC_z_HCP,1));
s_FC_z_HCP = squeeze(nanstd(FC_z_HCP));

%plotFCgeneral(m_FC_z_HCP,s_FC_z_HCP,strcat('HCP-Average','-FC-Voxels-AAL-ROI-corr-HCP-RestingState-',ROI_Label));

clear this_corr this_corr_z m_sub_FC_z s_sub_FC_z

% distributions

label_one = 'Magdeburg';
label_two = 'PETRA';
label_three = 'HCP';

x_label = 'Correlations (Fisher-Z)';
y_label = 'Probability';

title_label = strcat('Runs','-',ROI_Label);

FC_MD_z(:,:,:) = (1/2).*log((1 + FC_MD)./(1 - FC_MD));
FC_PETRA_z(:,:,:) = (1/2).*log((1 + FC_PETRA)./(1 - FC_PETRA));
FC_HCP_z(:,:,:) = (1/2).*log((1 + FC_HCP)./(1 - FC_HCP));
min_val = min([FC_MD_z(:);FC_PETRA_z(:);FC_HCP_z(:)]);
max_val = max([FC_MD_z(:);FC_PETRA_z(:);FC_HCP_z(:)]);

saveDistributionDensity3(min_val,max_val,FC_MD_z(:),FC_PETRA_z(:),FC_HCP_z(:),label_one,label_two,label_three,title_label,x_label,y_label);

title_label = strcat('Subjects','-',ROI_Label);

min_val = min([m_sub_FC_z_all_MD(:);m_sub_FC_z_all_PETRA(:);m_sub_FC_z_all_HCP(:)]);
max_val = max([m_sub_FC_z_all_MD(:);m_sub_FC_z_all_PETRA(:);m_sub_FC_z_all_HCP(:)]);

saveDistributionDensity3(min_val,max_val,m_sub_FC_z_all_MD(:),m_sub_FC_z_all_PETRA(:),m_sub_FC_z_all_HCP(:),label_one,label_two,label_three,title_label,x_label,y_label);

title_label = strcat('Average','-',ROI_Label);

min_val = min([m_FC_z_MD(:);m_FC_z_PETRA(:);m_FC_z_HCP(:)]);
max_val = max([m_FC_z_MD(:);m_FC_z_PETRA(:);m_FC_z_HCP(:)]);

saveDistributionDensity3(min_val,max_val,m_FC_z_MD(:),m_FC_z_PETRA(:),m_FC_z_HCP(:),label_one,label_two,label_three,title_label,x_label,y_label);

close all

toc

end

function checkThreeLevelsOfFC_MD_OnlyDistribution_WholeBrain

nTotalRuns_MD = 32;
nTotalRuns_PETRA = 27;
nTotalRuns_HCP = 32;

nRuns_MD = 4;
nRuns_PETRA = 3;
nRuns_HCP = 4;

nSubjects_MD = 8;
nSubjects_PETRA = 9;
nSubjects_HCP = 8;

nROI = 90;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

% Magdeburg

for idx_ROI=1:nROI
    
    ROI = AAL_ROI(idx_ROI).ID;
    ROI_Label = AAL_ROI(idx_ROI).Nom_L;
    ROI_Label = strrep(ROI_Label,'_','-');
    idx_voxels = find(AAL_img == ROI);
    nVoxels = length(idx_voxels);
    
    iiRun = 0;
    for iSubject=1:nSubjects_MD

        load(strcat('LHR-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-RestingState-',ROI_Label,'.mat'));

        for iRun=1:nRuns_MD

            iiRun = iiRun + 1;

            FC_MD(idx_ROI).corr(iiRun,:,:) = FC_Voxels.run(iRun).rho_RestingState;

            sub_FC_MD(iRun,:,:) = FC_Voxels.run(iRun).rho_RestingState;

        end

        sub_FC_z_MD(:,:,:) = (1/2).*log((1 + sub_FC_MD)./(1 - sub_FC_MD));

        m_sub_FC_z = squeeze(nanmean(sub_FC_z_MD,1));
        
        m_sub_FC_z_all_MD(idx_ROI).corr(iSubject,:,:) = m_sub_FC_z(:,:);

        clear FC_Voxels sub_FC_MD sub_FC_z_MD m_sub_FC_z

    end

end 

m_sub_FC_z_all_MD_v = [];

for idx_ROI=1:nROI
    
    m_sub_FC_z_all_MD_v = [m_sub_FC_z_all_MD_v(:); m_sub_FC_z_all_MD(idx_ROI).corr(:)];
    
end

FC_z_MD = [];

for idx_ROI=1:nROI
    
    FC_z_MD = [FC_z_MD(:); (1/2).*log((1 + FC_MD(idx_ROI).corr(:))./(1 - FC_MD(idx_ROI).corr(:)))];
    
end

m_FC_z_MD = [];

for idx_ROI=1:nROI
    
    f = (1/2).*log((1 + FC_MD(idx_ROI).corr)./(1 - FC_MD(idx_ROI).corr));
    
    m = nanmean(f,1);
    
    m_FC_z_MD = [m_FC_z_MD(:); m(:)];
    
end

clear this_corr this_corr_z m_sub_FC_z s_sub_FC_z

% PETRA

all_settings = getAllSettingsPetra;

for idx_ROI=1:nROI
    
    ROI = AAL_ROI(idx_ROI).ID;
    ROI_Label = AAL_ROI(idx_ROI).Nom_L;
    ROI_Label = strrep(ROI_Label,'_','-');
    idx_voxels = find(AAL_img == ROI);
    nVoxels = length(idx_voxels);
    
    iiRun = 0;
    for iSubject=1:nSubjects_PETRA

        settings = all_settings(iSubject).settings;

        load(strcat('PETRA-',settings.codes.subject,'-FC-Voxels-AAL-ROI-corr-Petra-RestingState-',ROI_Label,'.mat'));

        for iRun=1:nRuns_PETRA

            iiRun = iiRun + 1;

            FC_PETRA(idx_ROI).corr(iiRun,:,:) = FC_Voxels.run(iRun).rho_RestingState;

            sub_FC_PETRA(iRun,:,:) = FC_Voxels.run(iRun).rho_RestingState;

        end

        sub_FC_z_PETRA(:,:,:) = (1/2).*log((1 + sub_FC_PETRA)./(1 - sub_FC_PETRA));

        m_sub_FC_z = squeeze(nanmean(sub_FC_z_PETRA,1));
        
        m_sub_FC_z_all_PETRA(idx_ROI).corr(iSubject,:,:) = m_sub_FC_z(:,:);

        clear FC_Voxels sub_FC_PETRA sub_FC_z_PETRA m_sub_FC_z

    end

end

m_sub_FC_z_all_PETRA_v = [];

for idx_ROI=1:nROI
    
    m_sub_FC_z_all_PETRA_v = [m_sub_FC_z_all_PETRA_v(:); m_sub_FC_z_all_PETRA(idx_ROI).corr(:)];
    
end

FC_z_PETRA = [];

for idx_ROI=1:nROI
    
    FC_z_PETRA = [FC_z_PETRA(:); (1/2).*log((1 + FC_PETRA(idx_ROI).corr(:))./(1 - FC_PETRA(idx_ROI).corr(:)))];
    
end

m_FC_z_PETRA = [];

for idx_ROI=1:nROI
    
    f = (1/2).*log((1 + FC_PETRA(idx_ROI).corr)./(1 - FC_PETRA(idx_ROI).corr));
    
    m = nanmean(f,1);
    
    m_FC_z_PETRA = [m_FC_z_PETRA(:); m(:)];
    
end

clear this_corr this_corr_z m_sub_FC_z s_sub_FC_z

% HCP

all_settings = getAllSettingsHCP;

for idx_ROI=1:nROI
    
    ROI = AAL_ROI(idx_ROI).ID;
    ROI_Label = AAL_ROI(idx_ROI).Nom_L;
    ROI_Label = strrep(ROI_Label,'_','-');
    idx_voxels = find(AAL_img == ROI);
    nVoxels = length(idx_voxels);
    
    iiRun = 0;
    for iSubject=1:nSubjects_HCP

        settings = all_settings(iSubject).settings;

        for iRun=1:nRuns_HCP

            load(strcat('HCP-SUBJ',int2str(iSubject),'-FC-Voxels-AAL-ROI-corr-HCP-RestingState-',ROI_Label,'-',int2str(iRun),'.mat'));

            iiRun = iiRun + 1;

            FC_HCP(idx_ROI).corr(iiRun,:,:) = FC_Voxels.run(iRun).rho_RestingState;

            sub_FC_HCP(iRun,:,:) = FC_Voxels.run(iRun).rho_RestingState;

        end

        sub_FC_z_HCP(:,:,:) = (1/2).*log((1 + sub_FC_HCP)./(1 - sub_FC_HCP));

        m_sub_FC_z = squeeze(nanmean(sub_FC_z_HCP,1));
        
        m_sub_FC_z_all_HCP(idx_ROI).corr(iSubject,:,:) = m_sub_FC_z(:,:);

        clear FC_Voxels sub_FC_HCP sub_FC_z_HCP m_sub_FC_z

    end

end

m_sub_FC_z_all_HCP_v = [];

for idx_ROI=1:nROI
    
    m_sub_FC_z_all_HCP_v = [m_sub_FC_z_all_HCP_v(:); m_sub_FC_z_all_HCP(idx_ROI).corr(:)];
    
end

FC_z_HCP = [];

for idx_ROI=1:nROI
    
    FC_z_HCP = [FC_z_HCP(:); (1/2).*log((1 + FC_HCP(idx_ROI).corr(:))./(1 - FC_HCP(idx_ROI).corr(:)))];
    
end

m_FC_z_HCP = [];

for idx_ROI=1:nROI
    
    f = (1/2).*log((1 + FC_HCP(idx_ROI).corr)./(1 - FC_HCP(idx_ROI).corr));
    
    m = nanmean(f,1);
    
    m_FC_z_HCP = [m_FC_z_HCP(:); m(:)];
    
end

clear this_corr this_corr_z m_sub_FC_z s_sub_FC_z

% distributions

label_one = 'Magdeburg';
label_two = 'PETRA';
label_three = 'HCP';

x_label = 'Correlations (Fisher-Z)';
y_label = 'Probability';

title_label = strcat('Runs Level','-','Whole-Brain');

min_val = min([FC_MD_z(:);FC_PETRA_z(:);FC_HCP_z(:)]);
max_val = max([FC_MD_z(:);FC_PETRA_z(:);FC_HCP_z(:)]);

plotDistributionDensity3(min_val,max_val,FC_MD_z(:),FC_PETRA_z(:),FC_HCP_z(:),label_one,label_two,label_three,title_label,x_label,y_label);

title_label = strcat('Subjects Level','-','Whole-Brain');

min_val = min([m_sub_FC_z_all_MD_v(:);m_sub_FC_z_all_PETRA_v(:);m_sub_FC_z_all_HCP_v(:)]);
max_val = max([m_sub_FC_z_all_MD_v(:);m_sub_FC_z_all_PETRA_v(:);m_sub_FC_z_all_HCP_v(:)]);

plotDistributionDensity3(min_val,max_val,m_sub_FC_z_all_MD_v(:),m_sub_FC_z_all_PETRA_v(:),m_sub_FC_z_all_HCP_v(:),label_one,label_two,label_three,title_label,x_label,y_label);

title_label = strcat('Average','-','Whole-Brain');

min_val = min([m_FC_z_MD_v(:);m_FC_z_PETRA_v(:);m_FC_z_HCP_v(:)]);
max_val = max([m_FC_z_MD_v(:);m_FC_z_PETRA_v(:);m_FC_z_HCP_v(:)]);

plotDistributionDensity3(min_val,max_val,m_FC_z_MD_v(:),m_FC_z_PETRA_v(:),m_FC_z_HCP_v(:),label_one,label_two,label_three,title_label,x_label,y_label);

close all

end

function checkDataWithWithOUTFieldMap

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nTotalClusters = 758;
nROIs = 90;
nTR = 150;
MNI_dim = [91 109 91];

folder_field = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-1-22-10-2015\preprocessed\T2-Stimulus-RestingState\Run-1-4-1\FSL\custom';
folder_NO_field = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-1-22-10-2015\preprocessed\T2-Stimulus-RestingState\Run-1-4-1\FSL\Melodic.ica';

file = 'filtered_func_data_mcf_unwarp2standard-clean-voxel-res';
file_no = 'filtered_func_data2standard-clean-voxel-res';

data = nifti(strcat(folder_field,'\',file,'.nii'));
data_no = nifti(strcat(folder_NO_field,'\',file_no,'.nii'));

data.dat.fname = strcat(folder_field,'\',file,'.nii');
data_no.dat.fname = strcat(folder_NO_field,'\',file_no,'.nii');

FIELD = data.dat(:,:,:,:);
NO_FIELD = data_no.dat(:,:,:,:);

clusters_field = zeros(nTR,nTotalClusters);
clusters_field_no = zeros(nTR,nTotalClusters);

iiCluster = 0;
for iROI=1:nROIs
    
    nClusters = length(ROI(iROI).clusters);
   
    for iCluster=1:nClusters
        
        iiCluster = iiCluster + 1;
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
        
        nVoxels = length(idx_voxels);
        
        voxels_ts = zeros(nVoxels,nTR);
        voxels_ts_no = zeros(nVoxels,nTR);
        
        for iVoxel=1:nVoxels
           
            [idxx,idxy,idxz] = ind2sub(MNI_dim,idx_voxels(iVoxel));
            
            voxels_ts(iVoxel,1:nTR) = FIELD(idxx,idxy,idxz,1:nTR);
            voxels_ts_no(iVoxel,1:nTR) = NO_FIELD(idxx,idxy,idxz,1:nTR);
            
        end
        
        m_voxels_ts = mean(voxels_ts,1);
        m_voxels_ts_no = mean(voxels_ts_no,1);
    
        clusters_field(1:nTR,iiCluster) = m_voxels_ts(:);
        clusters_field_no(1:nTR,iiCluster) = m_voxels_ts_no(:);
        
    end
    
end

FIELD_corr = corr(clusters_field);
NO_FIELD_corr = corr(clusters_field_no);

save('FieldMap-Comparison.mat','FIELD_corr','NO_FIELD_corr');

%%% AND PLOT

mu_z = (1/2).*log((1 + FIELD_corr)./(1 - FIELD_corr));
label = 'FIELD';
plotFCgeneral2(mu_z,label);

mu_z = (1/2).*log((1 + NO_FIELD_corr)./(1 - NO_FIELD_corr));
label = 'NO-FIELD';
plotFCgeneral2(mu_z,label);

end

%%% CC400 PARCELLATION

function getCorrFromC400Parcellation

nClusters = 400;
nTR = 150;
nRuns = 4;
nTotalRuns = 32;
MNI_size = [91 109 91];
MNI_size_mm = [47 56 46];
mm = 4;

all_settings = getAllSettings;
nSubjects = length(all_settings);

mean_area_RestingState_voxel = zeros(nClusters,nTR);

iiRun = 0;
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
    
    kind = 'RestingState';
    for irun=1:nRuns;
        iiRun = iiRun + 1;
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL_mm(settings,kind,irun,file,mask,get_at_this_preprocessed_step,mm);
    end
    
end

%%% C400
file = 'C400_4mm_tcorr_2level.nii';
img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\craddock_2011_parcellations\',file));
img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\craddock_2011_parcellations\',file);
C400 = img.dat(:,:,:);

for iiRun=1:nTotalRuns
    
    for iCluster=1:nClusters
    
        idx_voxels = find(C400==iCluster); 
        nVoxels = length(idx_voxels);
   
        this_cluster_voxels = zeros(nVoxels,nTR);
        
        for iVoxel=1:nVoxels
        
            [idxx,idxy,idxz] = ind2sub(MNI_size_mm,idx_voxels(iVoxel));
            
            this_cluster_voxels(iVoxel,:) = RestingState(iiRun).run(idxx,idxy,idxz,:);
            
        end
        
        m_this_cluster_voxels = mean(this_cluster_voxels,1);
        
        all_clusters(iiRun).run(:,iCluster) = m_this_cluster_voxels(:);
   
    end

end

for iiRun=1:nTotalRuns
   
    all_corr(iiRun,:,:) = corr(all_clusters(iiRun).run);
    
end

save('C400-4mm-tcorr-2level-Corr.mat','all_corr');

end

function getCorrFromC400ParcellationInsideParcelVoxelLevel

nClusters = 400;
nTR = 150;
nRuns = 4;
nTotalRuns = 32;
MNI_size = [91 109 91];
MNI_size_mm = [47 56 46];
mm = 4;

all_settings = getAllSettings;
nSubjects = length(all_settings);

mean_area_RestingState_voxel = zeros(nClusters,nTR);

iiRun = 0;
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
    
    kind = 'RestingState';
    for irun=1:nRuns;
        iiRun = iiRun + 1;
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL_mm(settings,kind,irun,file,mask,get_at_this_preprocessed_step,mm);
    end
    
end

%%% C400
file = 'C400_4mm_tcorr_2level.nii';
img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\craddock_2011_parcellations\',file));
img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\craddock_2011_parcellations\',file);
C400 = img.dat(:,:,:);

for iiRun=1:nTotalRuns
    
    disp(int2str(iiRun));
    
    all_corr(iiRun).corr = [];
    
    for iCluster=1:nClusters
       
        idx_voxels = find(C400==iCluster); 
        nVoxels = length(idx_voxels);
   
        this_cluster_voxels = zeros(nTR,nVoxels);
        
        for iVoxel=1:nVoxels
        
            [idxx,idxy,idxz] = ind2sub(MNI_size_mm,idx_voxels(iVoxel));
            
            this_cluster_voxels(:,iVoxel) = RestingState(iiRun).run(idxx,idxy,idxz,:);
            
        end
        
        if ~isempty(idx_voxels); 
            
            corr_this_cluster_voxels = corr(this_cluster_voxels);
        
            all_corr(iiRun).corr = [all_corr(iiRun).corr, corr_this_cluster_voxels(:)'];
   
        end
        
    end

end

save('C400-4mm-2level-Inside-Corr.mat','all_corr','-v7.3');

end

function getCorrFromC400ParcellationInsideParcelVoxelLevelsavepercluster

nClusters = 400;
nTR = 150;
nRuns = 4;
nTotalRuns = 32;
MNI_size = [91 109 91];
MNI_size_4mm = [47 56 46];
mm = 4;

all_settings = getAllSettings;
nSubjects = length(all_settings);

mean_area_RestingState_voxel = zeros(nClusters,nTR);

iiRun = 0;
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
    
    kind = 'RestingState';
    for irun=1:nRuns;
        iiRun = iiRun + 1;
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL_mm(settings,kind,irun,file,mask,get_at_this_preprocessed_step,mm);
    end
    
end

%%% C400
file = 'C400_4mm_tcorr_2level.nii';
img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\craddock_2011_parcellations\',file));
img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\craddock_2011_parcellations\',file);
C400 = img.dat(:,:,:);
    
for iCluster=1:nClusters

    idx_voxels = find(C400==iCluster); 
    nVoxels = length(idx_voxels);

    clusters_corr(iCluster).corr = zeros(nTotalRuns,nVoxels,nVoxels);
    
    if ~isempty(idx_voxels)

        all_idx_zeros = [];
        for iiRun=1:nTotalRuns

            this_cluster_voxels = zeros(nTR,nVoxels);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_size_4mm,idx_voxels(iVoxel));

                this_cluster_voxels(:,iVoxel) = RestingState(iiRun).run(idxx,idxy,idxz,:);

            end

            s_voxels = sum(this_cluster_voxels,1);
            idx_zeros = find(s_voxels==0);
            all_idx_zeros = [all_idx_zeros,idx_zeros];

            clusters_corr(iCluster).corr(iiRun,:,:) = corr(this_cluster_voxels);

        end

        if ~isempty(all_idx_zeros)

            clusters_corr(iCluster).corr(:,unique(all_idx_zeros),:) = [];
            clusters_corr(iCluster).corr(:,:,unique(all_idx_zeros)) = [];

        end

        if iCluster < 5; plotFCVoxelsInsideAClusterC400(clusters_corr(iCluster).corr,iCluster); end
    
    end
    
end

save('C400-4mm-tcorr-2level-Inside-Corr.mat','clusters_corr','-v7.3');

end

function getAndPlotVarianceDistributionFromC400Parcellation

load('C400-4mm-tcorr-2level-Corr.mat');

nClusters = 400;
nTotalRuns = 32;

all_corr_z = (1/2).*log((1 + all_corr)./(1 - all_corr));
            
my_std = squeeze(std(all_corr_z,0,1));
my_mean = squeeze(mean(all_corr_z,1));

vector_one = my_std;
vector_two = my_mean;
label_one = 'STD';
label_two = 'Mean';
title_label = 'C400-tcorr-2level';

MyContourCorrelation_v5(vector_two(:),vector_one(:),title_label,1);

end

function getAndPlotVarianceDistributionFromC400ParcellationInsideParcelVoxelLevel

load('C400-4mm-tcorr-2level-Inside-Corr.mat');

nClusters = 400;
nTotalRuns = 32;

for iRun=1:nTotalRuns
   
    iPair = 0;
    for iCluster=1:nClusters
        
        nVoxels = size(clusters_corr(iCluster).corr,2);
    
        for iVoxel=1:nVoxels
            
            for iiVoxel=1:nVoxels
                
                iPair = iPair + 1;
                
                all_pairs(iRun,iPair) = clusters_corr(iCluster).corr(iRun,iVoxel,iiVoxel);
                
            end
            
        end
        
    end
    
end

all_corr_z = (1/2).*log((1 + all_pairs)./(1 - all_pairs));
            
my_std = squeeze(std(all_corr_z,0,1));
my_mean = squeeze(mean(all_corr_z,1));

vector_one = my_std;
vector_two = my_mean;
label_one = 'STD';
label_two = 'Mean';
title_label = 'C400-4mm-tcorr-2level-Inside';

MyContourCorrelation_v5(vector_two(:),vector_one(:),title_label,1);

end

function plotFCVoxelsInsideAClusterC400(voxels_corr,idx_cluster)

nRuns = 32;
nTR = 150;
nVoxels = size(voxels_corr,2);
nTotalRuns = 32;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;
for iVoxel=1:nVoxels
    
    for iiVoxel=1:nVoxels
    
        if iVoxel ~= iiVoxel
            
            iPair = iPair + 1;
            for iRun=1:nTotalRuns

                single_rho = squeeze(voxels_corr(iRun,iVoxel,iiVoxel));

                pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

                pairs_rest_rho(iRun,iPair) = single_rho;

            end
        
        end
    
    end

end
            
m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ mean(pairs_rest_z,1);


mean_corr = zeros(nVoxels,nVoxels);
s_corr = zeros(nVoxels,nVoxels);

iPair = 0;
for iVoxel=1:nVoxels
    
    for iiVoxel=1:nVoxels
        
        if iVoxel ~= iiVoxel
            
            iPair = iPair + 1;
        
            %if my_cv_aal(iPair) < cv_criterion & my_cv_aal(iPair) > 0 & p_pairs(iPair) < p_criterion
                
                mean_corr(iVoxel,iiVoxel) =  m_pairs_z(iPair);
                s_corr(iVoxel,iiVoxel) =  s_pairs_z(iPair);
                
            %end
            
        end
        
    end
    
end

new_label = strcat('FC-C400-Voxels-Cluster-',int2str(idx_cluster),'-','Voxels','-',int2str(nVoxels));
plotFCgeneral(mean_corr,s_corr,new_label);

close all;

end

function plotMeanFCC400

load('C400-4mm-tcorr-2level-Corr.mat');

all_corr_z = (1/2).*log((1 + all_corr)./(1 - all_corr));
            
my_std = squeeze(std(all_corr_z,0,1));
my_mean = squeeze(mean(all_corr_z,1));


new_label = strcat('FC-C400');
plotFCgeneral(my_mean,my_std,new_label);

end

function getInformationC400

load('C400-4mm-tcorr-2level-Corr.mat');

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

all_corr_z = (1/2).*log((1 + all_corr)./(1 - all_corr));
            
my_std = squeeze(nanstd(all_corr_z,0,1));
my_mean = squeeze(nanmean(all_corr_z,1));

new_label = strcat('Information-C400');
biInfo = TotalInformationContent( my_mean, my_std , z_threshold, cv_threshold );

save('C400-biInfo.mat','biInfo');

end


%%% VAN ESSEN PARCELLATION

function convert180to360VanEssen

nClusters = 180;
MNI_size = [91 109 91];

file = 'HCP-MMP1_on_MNI152_ICBM2009a_nlin-2mm.nii';

img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file));
img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file);
vanEssen = img.dat(:,:,:);

for iCluster=1:nClusters
    
    idx_voxels = find(vanEssen==iCluster);
    
    nVoxels = length(idx_voxels);
    
    for iVoxel=1:nVoxels
        
        [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
        
        if idxx > MNI_size(1)/2
            
            vanEssen(idxx,idxy,idxz) = vanEssen(idxx,idxy,idxz) + nClusters;
            
        end
        
    end
    
end

nifti_file = img;
offset = img.dat.offset;
scl_slope = img.dat.scl_slope;
scl_inter = img.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = img.dat.dim;

descrip = 'VanEssen-LR';

fname = strcat('VanEssen-LR-2mm','.nii');
input_data = vanEssen; 
real_save_image;

end

function getCorrFromVanEssenParcellation

%nClusters = 180;
nClusters = 360;
nTR = 150;
nRuns = 4;
nTotalRuns = 32;
MNI_size = [91 109 91];

all_settings = getAllSettings;
nSubjects = length(all_settings);

mean_area_RestingState_voxel = zeros(nClusters,nTR);

iiRun = 0;
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
    
    kind = 'RestingState';
    for irun=1:nRuns;
        iiRun = iiRun + 1;
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
end

%%% 180
% file = 'HCP-MMP1_on_MNI152_ICBM2009a_nlin-2mm.nii';
% img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file));
% img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file);
% vanEssen = img.dat(:,:,:);

%%% 360
file = 'VanEssen-LR-2mm.nii';
img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file));
img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file);
vanEssen = img.dat(:,:,:);

for iiRun=1:nTotalRuns
    
    for iCluster=1:nClusters
    
        idx_voxels = find(vanEssen==iCluster); 
        nVoxels = length(idx_voxels);
   
        this_cluster_voxels = zeros(nVoxels,nTR);
        
        for iVoxel=1:nVoxels
        
            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
            
            this_cluster_voxels(iVoxel,:) = RestingState(iiRun).run(idxx,idxy,idxz,:);
            
        end
        
        m_this_cluster_voxels = mean(this_cluster_voxels,1);
        
        all_clusters(iiRun).run(:,iCluster) = m_this_cluster_voxels(:);
   
    end

end

for iiRun=1:nTotalRuns
   
    all_corr(iiRun,:,:) = corr(all_clusters(iiRun).run);
    
end

save('VanEssen-Corr.mat','all_corr');

end

function getCorrFromVanEssenParcellationInsideParcelVoxelLevel

%nClusters = 180;
nClusters = 360;
nTR = 150;
nRuns = 4;
nTotalRuns = 32;
MNI_size = [91 109 91];

all_settings = getAllSettings;
nSubjects = length(all_settings);

mean_area_RestingState_voxel = zeros(nClusters,nTR);

iiRun = 0;
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
    
    kind = 'RestingState';
    for irun=1:nRuns;
        iiRun = iiRun + 1;
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
end

%%% 180
% file = 'HCP-MMP1_on_MNI152_ICBM2009a_nlin-2mm.nii';
% img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file));
% img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file);
% vanEssen = img.dat(:,:,:);

%%% 360
file = 'VanEssen-LR-2mm.nii';
img = nifti(strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file));
img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\Glasser_et_al._2016-parcellation\',file);
vanEssen = img.dat(:,:,:);

for iiRun=1:nTotalRuns
    
    all_corr(iiRun).corr = [];
    
    for iCluster=1:nClusters
        
        disp(int2str(iCluster));
    
        idx_voxels = find(vanEssen==iCluster); 
        nVoxels = length(idx_voxels);
   
        this_cluster_voxels = zeros(nTR,nVoxels);
        
        for iVoxel=1:nVoxels
        
            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
            
            this_cluster_voxels(:,iVoxel) = RestingState(iiRun).run(idxx,idxy,idxz,:);
            
        end
        
        corr_this_cluster_voxels = corr(this_cluster_voxels);
        
        all_corr(iiRun).corr = [all_corr(iiRun).corr, corr_this_cluster_voxels(:)'];
   
    end

end

save('VanEssen-Inside-Corr.mat','all_corr','-v7.3');

end

function getAndPlotVarianceDistributionFromVanEssenParcellation

load('VanEssen-360-Corr.mat');

% nClusters = 180;
nClusters = 360;
nTotalRuns = 32;

all_corr_z = (1/2).*log((1 + all_corr)./(1 - all_corr));
            
my_std = squeeze(std(all_corr_z,0,1));
my_mean = squeeze(mean(all_corr_z,1));

vector_one = my_std;
vector_two = my_mean;
label_one = 'STD';
label_two = 'Mean';
title_label = 'VanEssen';

MyContourCorrelation_v5(vector_two(:),vector_one(:),title_label,1);

end

function getAndPlotVarianceDistributionFromVanEssenParcellationInsideParcelVoxelLevel

load('VanEssen-Inside-Corr.mat');

% nClusters = 180;
nClusters = 360;
nTotalRuns = 32;

all_corr_struct = all_corr;
clear all_corr;

all_corr = zeros(nTotalRuns,length(all_corr_struct(1).corr));

for iRun=1:nTotalRuns
    
    all_corr(iRun,:) = all_corr_struct(iRun).corr(:);
    
end

all_corr_z = (1/2).*log((1 + all_corr)./(1 - all_corr));
            
my_std = squeeze(std(all_corr_z,0,1));
my_mean = squeeze(mean(all_corr_z,1));

vector_one = my_std;
vector_two = my_mean;
label_one = 'STD';
label_two = 'Mean';
title_label = 'VanEssen-Inside';

MyContourCorrelation_v5(vector_two(:),vector_one(:),title_label,1);

end

function plotMeanFCVanEssen

load('VanEssen-360-Corr.mat');

% nClusters = 180;
nClusters = 360;
nTotalRuns = 32;

all_corr_z = (1/2).*log((1 + all_corr)./(1 - all_corr));
            
my_std = squeeze(std(all_corr_z,0,1));
my_mean = squeeze(mean(all_corr_z,1));


new_label = strcat('FC-VanEssen');
plotFCgeneral(my_mean,my_std,new_label);

end

function getInformationVanEssen

load('VanEssen-360-Corr.mat');

% nClusters = 180;
nClusters = 360;
nTotalRuns = 32;

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

all_corr_z = (1/2).*log((1 + all_corr)./(1 - all_corr));
            
my_std = squeeze(nanstd(all_corr_z,0,1));
my_mean = squeeze(nanmean(all_corr_z,1));

new_label = strcat('Information-VanEssen');
biInfo = TotalInformationContent( my_mean, my_std , z_threshold, cv_threshold );

save('VanEssen-biInfo.mat','biInfo');

end

%%% INDIVIDUAL VOXELS

function getCorrFromFromVoxelsFromASampleOf20VoxelsPerCluster

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nClusters = 758;
nROI = 90;
nRuns = 4;
MNI_img = [91 109 91];
nTotalVoxels = 160990;
nTotalClusters = 758;
nClusterSize = 200;
nSamples = 20;
nTR = 150;

all_idx_voxels = [];

for iROI=1:nROI
   
    nClusters = length(ROI(iROI).clusters);
    
    for iCluster=1:nClusters
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels(1:nSamples);
        
        all_idx_voxels = [all_idx_voxels;idx_voxels];
        
    end
    
end

all_settings = getAllSettings;

nSubjects = length(all_settings);

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

%%% REST

voxels_corr_rest = zeros(nRuns*nSubjects,length(all_idx_voxels),length(all_idx_voxels));

iiRun = 0;

for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
    
    kind = 'RestingState';
    for irun=1:nRuns;
        [RestingState(irun).run, RestingState(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
    for irun=1:nRuns
        
        iiRun = iiRun + 1;
        
        rest_sample_voxels = zeros(nTR,nSamples*nTotalClusters);
        
        disp(strcat('Run:',int2str(iiRun)));
        
        for iVoxel=1:length(all_idx_voxels)

            [idxx,idxy,idxz] = ind2sub(size(AAL_img),all_idx_voxels(iVoxel));
            
            rest_sample_voxels(:,iVoxel) = RestingState(irun).run(idxx,idxy,idxz,:);
            
        end
        
        voxels_corr_rest(iiRun,:,:) = corr(rest_sample_voxels);
        
    end

end

disp('I will save right now...');

save('Samples-Of-Voxels-Corr-Rest.mat','voxels_corr_rest','-v7.3');

%%% PASSIVE

voxels_corr_passive = zeros(nRuns*nSubjects,length(all_idx_voxels),length(all_idx_voxels));

iiRun = 0;

for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
        
    kind = 'Passive';
    for irun=1:nRuns;
        [Passive(irun).run, Passive(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
    for irun=1:nRuns
        
        iiRun = iiRun + 1;
        
        passive_sample_voxels = zeros(nTR,nSamples*nTotalClusters);
        
        disp(strcat('Run:',int2str(iiRun)));
        
        for iVoxel=1:length(all_idx_voxels)

            [idxx,idxy,idxz] = ind2sub(size(AAL_img),all_idx_voxels(iVoxel));
            
            passive_sample_voxels(:,iVoxel) = Passive(irun).run(idxx,idxy,idxz,:);
            
        end
        
        voxels_corr_passive(iiRun,:,:) = corr(passive_sample_voxels);
        
    end

end

disp('I will save right now...');

save('Samples-Of-Voxels-Corr-Passive.mat','voxels_corr_passive','-v7.3');

%%% TRACK

voxels_corr_track = zeros(nRuns*nSubjects,length(all_idx_voxels),length(all_idx_voxels));

iiRun = 0;

for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
        
    kind = 'Track';
    for irun=1:nRuns;
        [Track(irun).run, Track(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
    for irun=1:nRuns
        
        iiRun = iiRun + 1;
        
        track_sample_voxels = zeros(nTR,nSamples*nTotalClusters);
        
        disp(strcat('Run:',int2str(iiRun)));
        
        for iVoxel=1:length(all_idx_voxels)

            [idxx,idxy,idxz] = ind2sub(size(AAL_img),all_idx_voxels(iVoxel));
            
            track_sample_voxels(:,iVoxel) = Track(irun).run(idxx,idxy,idxz,:);
            
        end
        
        voxels_corr_track(iiRun,:,:) = corr(track_sample_voxels);
        
    end

end

disp('I will save right now...');

save('Samples-Of-Voxels-Corr-Track.mat','voxels_corr_track','-v7.3');

end

function plotVarianceDistributionVoxelsFromASampleOf20VoxelsPerCluster(condition_label)

load(strcat('Samples-Of-Voxels-Corr','-',condition_label,'.mat'));
load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

eval(strcat('voxels_corr = voxels_corr','_',lower(condition_label),';'));

nROI = 90;
nRuns = 32;
nSamplesOfVoxels = 20;

nTotalVoxels = 758*20;

iPair_same_cluster = 0;
iPair_same_aal = 0;
iPair_diff_aal = 0;
IVoxel = 0;
for iROI=1:nROI
    nClusters = length(ROI(iROI).clusters);
    for iCluster=1:nClusters
       for iVoxel=1:nSamplesOfVoxels
          IVoxel = IVoxel + 1;
          if mod(IVoxel,1000) == 0; disp(int2str(IVoxel)); end
          IIVoxel = 0;
          for iiROI=1:nROI
              nnClusters = length(ROI(iiROI).clusters);
              for iiCluster=1:nnClusters
                  for iiVoxel=1:nSamplesOfVoxels
                      IIVoxel = IIVoxel + 1;
                      if IIVoxel >= IVoxel
                          if iROI == iiROI & iCluster == iiCluster 
                              iPair_same_cluster = iPair_same_cluster + 1;
                          end
                          if iROI == iiROI & iCluster ~= iiCluster
                              iPair_same_aal = iPair_same_aal + 1;
                          end
                          if iROI ~= iiROI
                              iPair_diff_aal = iPair_diff_aal + 1;
                          end
                      end
                  end
              end
          end
       end
    end
end

pairs_same_cluster = zeros(nRuns,iPair_same_cluster);
pairs_same_aal = zeros(nRuns,iPair_same_aal);
pairs_diff_aal = zeros(nRuns,iPair_diff_aal);

iPair_same_cluster = 0;
iPair_same_aal = 0;
iPair_diff_aal = 0;
IVoxel = 0;
for iROI=1:nROI
    
    nClusters = length(ROI(iROI).clusters);
   
    for iCluster=1:nClusters
       
       for iVoxel=1:nSamplesOfVoxels
           
          IVoxel = IVoxel + 1;
          
          if mod(IVoxel,1000) == 0; disp(int2str(IVoxel)); end
          
          IIVoxel = 0;
          
          for iiROI=1:nROI
              
              nnClusters = length(ROI(iiROI).clusters);
              
              for iiCluster=1:nnClusters

                  for iiVoxel=1:nSamplesOfVoxels
                      
                      IIVoxel = IIVoxel + 1;
                      
                      if IIVoxel >= IVoxel
                    
                          if iROI == iiROI & iCluster == iiCluster 

                              iPair_same_cluster = iPair_same_cluster + 1;

                               single_rho = squeeze(voxels_corr(:,IVoxel,IIVoxel));
            
                               pairs_same_cluster(:,iPair_same_cluster) = (1/2)*log((1 + single_rho)./(1 - single_rho));

                          end

                          if iROI == iiROI & iCluster ~= iiCluster

                              iPair_same_aal = iPair_same_aal + 1;

                              single_rho = squeeze(voxels_corr(:,IVoxel,IIVoxel));
                              
                              pairs_same_aal(:,iPair_same_aal) = (1/2)*log((1 + single_rho)./(1 - single_rho));

                          end

                          if iROI ~= iiROI

                              iPair_diff_aal = iPair_diff_aal + 1;

                              single_rho = squeeze(voxels_corr(:,IVoxel,IIVoxel));
                              
                              pairs_diff_aal(:,iPair_diff_aal) = (1/2)*log((1 + single_rho)./(1 - single_rho));

                          end
                      
                      end

                  end
                  
              end
              
          end
           
       end
        
    end
    
end

cv_threshold = 0.5;
% r_threshold = 0.186728;  % corresponds to p-value = 0.01 for N=150
r_threshold = 0.3;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );

my_var_same_cluster = nanvar(pairs_same_cluster);

my_cv_same_cluster = nanstd(pairs_same_cluster) ./ abs(nanmean(pairs_same_cluster,1));

my_mean_same_cluster = nanmean(pairs_same_cluster,1);
my_std_same_cluster = nanstd(pairs_same_cluster);

s = sum(pairs_same_cluster,1);
idx_zeros = find(s==0);

my_var_same_cluster(idx_zeros) = [];
my_cv_same_cluster(idx_zeros) = [];
my_mean_same_cluster(idx_zeros) = [];
my_std_same_cluster(idx_zeros) = [];

vector_one = my_std_same_cluster;
vector_two = my_mean_same_cluster;
vector_three = my_cv_same_cluster;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = strcat('Voxels-Same-Cluster','-',condition_label);

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);
fraction = FractionSignificantConsistent( vector_two, vector_one, z_threshold, cv_threshold );

my_var_same_aal = var(pairs_same_aal);

my_cv_same_aal = std(pairs_same_aal) ./ abs(mean(pairs_same_aal,1));

my_mean_same_aal = mean(pairs_same_aal,1);
my_std_same_aal = std(pairs_same_aal);

s = sum(pairs_same_aal,1);
idx_zeros = find(s==0);

my_var_same_aal(idx_zeros) = [];
my_cv_same_aal(idx_zeros) = [];
my_mean_same_aal(idx_zeros) = [];
my_std_same_aal(idx_zeros) = [];

vector_one = my_std_same_aal;
vector_two = my_mean_same_aal;
vector_three = my_cv_same_aal;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = strcat('Voxels-Same-AAL','-',condition_label);

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

my_var_diff_aal = var(pairs_diff_aal);

my_cv_diff_aal = std(pairs_diff_aal) ./ abs(mean(pairs_diff_aal,1));

my_mean_diff_aal = mean(pairs_diff_aal,1);
my_std_diff_aal = std(pairs_diff_aal);

s = sum(pairs_diff_aal,1);
idx_zeros = find(s==0);

my_var_diff_aal(idx_zeros) = [];
my_cv_diff_aal(idx_zeros) = [];
my_mean_diff_aal(idx_zeros) = [];
my_std_diff_aal(idx_zeros) = [];

vector_one = my_std_diff_aal;
vector_two = my_mean_diff_aal;
vector_three = my_cv_diff_aal;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = strcat('Voxels-Diff-AAL','-',condition_label);

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function getCorrVoxelsInsideACluster(idx_cluster)

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nROI = 90;
iiClusters = 0;
nTR = 150;
MNI_img = [91 109 91];
nTotalRuns = 32;
nRuns = 4;

all_settings = getAllSettings;
nSubjects = length(all_settings);

iiRun = 0;
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
        
    kind = 'RestingState';
    for irun=1:nRuns;
        iiRun = iiRun + 1;
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
end

for iROI=1:nROI
    
    nClusters = length(ROI(iROI).clusters);
    
    for iCluster=1:nClusters
        
        iiClusters = iiClusters + 1;
        
        if iiClusters == idx_cluster
            
            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
            
            nVoxels = length(idx_voxels);
            
            for iRun=1:nTotalRuns
            
                voxels_ts = zeros(nTR,nVoxels);
            
                for iVoxel=1:nVoxels
               
                    [idxx,idxy,idxz] = ind2sub(MNI_img,idx_voxels(iVoxel));
                
                    voxels_ts(:,iVoxel) = squeeze(RestingState(iRun).run(idxx,idxy,idxz,:));
                    
                end
                
                voxels_corr(iRun,:,:) = corr(voxels_ts);
                
            end
            
        end
        
    end

end

save(strcat('Voxels-Inside-Cluster-Corr-',int2str(idx_cluster),'.mat'),'voxels_corr');

end

function getCorrVoxelsInsideAClusterPetra

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nROI = 90;
iiClusters = 0;
nTR = 150;
MNI_img = [91 109 91];
nTotalRuns = 9;
nSegments = 3;

all_settings = getAllSettingsPetra;
nSubjects = length(all_settings);

for iSubject=1:nSubjects
    
    run = 1;
    
    settings = all_settings(iSubject).settings;
    
    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;

    [RestingState(iSubject).run, settings] = real_get_data_FSL_Petra(settings,run,file,get_at_this_preprocessed_step);
    
end

iiClusters = 0;
for iROI=1:nROI
    
    nClusters = length(ROI(iROI).clusters);
    
    for iCluster=1:nClusters
        
        iiClusters = iiClusters + 1;
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

        nVoxels = length(idx_voxels);

        iiRunSeg = 0;
        for iRun=1:nTotalRuns

            voxels_ts = zeros(nTR,nVoxels);

            for iSeg=1:nSegments
                
                iiRunSeg = iiRunSeg + 1;
                
                for iVoxel=1:nVoxels
                    
                    [idxx,idxy,idxz] = ind2sub(MNI_img,idx_voxels(iVoxel));

                    voxels_ts(:,iVoxel) = squeeze(RestingState(iRun).run(idxx,idxy,idxz,(1+(iSeg-1)*nTR):(nTR+(iSeg-1)*nTR)));

                end
                
                voxels_corr(iiRunSeg,:,:) = corr(voxels_ts);

            end
  

        end

        save(strcat('Voxels-Inside-Cluster-Petra-Corr-',int2str(iiClusters),'.mat'),'voxels_corr');
            
        clear voxels_corr
        
   end
        
end

end

function getCorrVoxelsInsideALLCluster

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nROI = 90;
iiClusters = 0;
nTR = 150;
MNI_img = [91 109 91];
nTotalRuns = 32;
nRuns = 4;

all_settings = getAllSettings;
nSubjects = length(all_settings);

iiRun = 0;
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
        
    kind = 'RestingState';
    for irun=1:nRuns;
        iiRun = iiRun + 1;
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
end

for iROI=1:nROI
    
    nClusters = length(ROI(iROI).clusters);
    
    for iCluster=1:nClusters
        
        iiClusters = iiClusters + 1;
        
        iiClusters
        
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

        nVoxels = length(idx_voxels);

        for iRun=1:nTotalRuns

            voxels_ts = zeros(nTR,nVoxels);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_img,idx_voxels(iVoxel));

                voxels_ts(:,iVoxel) = squeeze(RestingState(iRun).run(idxx,idxy,idxz,:));

            end

            voxels_corr(iRun,:,:) = corr(voxels_ts);

        end
        
        save(strcat('Voxels-Inside-Cluster-Corr-',int2str(iiClusters),'.mat'),'voxels_corr');
        
        clear voxels_corr

    end

end



end

function plotVarianceDistributionClustersVoxelsInsideALLCluster

nTotalClusters = 758;
nRuns = 32;

pairs = [];
for iCluster=1:nTotalClusters
    
    load(strcat('Voxels-Inside-Cluster-Corr-',int2str(iCluster)));
    nVoxels = size(voxels_corr,2);
    v = reshape(voxels_corr,[nRuns nVoxels*nVoxels]);
    clear voxels_corr
    pairs = [pairs,v];

end

pairs_z = (1/2).*log((1 + pairs)./(1 - pairs));

%%% EVERYBODY - VOXELS INSIDE - REST

my_var = var(pairs_z);

my_cv = std(pairs_z) ./ abs(mean(pairs_z,1));

my_mean = mean(pairs_z,1);
my_std = std(pairs_z);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Everybody-Voxels-Inside-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotVarianceDistributionClustersVoxelsInsideALLClusterPetra

nTotalClusters = 758;
nRuns = 27;

pairs = [];
for iCluster=1:nTotalClusters
    
    load(strcat('Voxels-Inside-Cluster-Petra-Corr-',int2str(iCluster)));
    nVoxels = size(voxels_corr,2);
    v = reshape(voxels_corr,[nRuns nVoxels*nVoxels]);
    clear voxels_corr
    pairs = [pairs,v];

end

pairs_z = (1/2).*log((1 + pairs)./(1 - pairs));

%%% EVERYBODY - VOXELS INSIDE - REST

my_var = var(pairs_z);

my_cv = std(pairs_z) ./ abs(mean(pairs_z,1));

my_mean = mean(pairs_z,1);
my_std = std(pairs_z);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Everybody-Voxels-Inside-Petra-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotFCVoxelsInsideACluster(idx_cluster)

load(strcat('Voxels-Inside-Cluster-Corr-',int2str(idx_cluster),'.mat'));

nRuns = 32;
nTR = 150;
nVoxels = size(voxels_corr,2);
nTotalRuns = 32;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;
for iVoxel=1:nVoxels
    
    for iiVoxel=1:nVoxels
    
        if iVoxel ~= iiVoxel
            
            iPair = iPair + 1;
            for iRun=1:nTotalRuns

                single_rho = squeeze(voxels_corr(iRun,iVoxel,iiVoxel));

                pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

                pairs_rest_rho(iRun,iPair) = single_rho;

            end
        
        end
    
    end

end
            
m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ mean(pairs_rest_z,1);


mean_corr = zeros(nVoxels,nVoxels);
s_corr = zeros(nVoxels,nVoxels);

iPair = 0;
for iVoxel=1:nVoxels
    
    for iiVoxel=1:nVoxels
        
        if iVoxel ~= iiVoxel
            
            iPair = iPair + 1;
        
            %if my_cv_aal(iPair) < cv_criterion & my_cv_aal(iPair) > 0 & p_pairs(iPair) < p_criterion
                
                mean_corr(iVoxel,iiVoxel) =  m_pairs_z(iPair);
                s_corr(iVoxel,iiVoxel) =  s_pairs_z(iPair);
                
            %end
            
        end
        
    end
    
end

new_label = strcat('FC-Voxels-Cluster-',int2str(idx_cluster),'-','Voxels','-',int2str(nVoxels));
plotFCgeneral(mean_corr,s_corr,new_label);

close all;

end

function getInformationRunAALVoxelLevel


nROIs = 90;
nSubjects = 8;
nRuns = 4;
nTotalRuns = 32;

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

all_settings = getAllSettings;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

biInfo = zeros(nSubjects,nROIs);
for iROI=1:nROIs
   
    area_Label = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    
    disp(strcat(area_Label,'-',datestr(now)));
    
    iiRun = 0;
    for iSubject=1:nSubjects
       
        settings = all_settings(iSubject).settings;
        
        load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-RestingState-',area_Label,'.mat'));
        
        for iRun=1:nRuns
            
            iiRun = iiRun + 1;
            
            single_rho = FC_Voxels.run(iRun).rho_RestingState;
            
            fisher(iRun,:,:) = (1/2).*log((1 + single_rho)./(1 - single_rho));
            
        end
        
        [l,c,z] = size(fisher);
        
        fisher_pairs = reshape(fisher,[l c*z]);
        
        mu_z = mean(fisher_pairs,1);
        sigma_z = std(fisher_pairs);
   
        biInfo(iSubject,iROI) = TotalInformationContent( mu_z, sigma_z, z_threshold, cv_threshold );
            
        clear FC_Voxels subject file area_label fisher fisher_pairs mu_z sigma_z
        
    end
    
end

save('LHR-biInfo.mat','biInfo');

end

function getInformationRunAALVoxelLevelPETRA


nROIs = 90;
nSubjects = 9;
nRuns = 3;
nTotalRuns = 27;

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

all_settings = getAllSettingsPetra;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

biInfo = zeros(nSubjects,nROIs);
for iROI=1:nROIs
   
    area_Label = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    
    disp(strcat(area_Label,'-',datestr(now)));
    
    iiRun = 0;
    for iSubject=1:nSubjects
       
        settings = all_settings(iSubject).settings;
        
        load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-Petra-RestingState-',area_Label,'.mat'));
        
        for iRun=1:nRuns
            
            iiRun = iiRun + 1;
            
            single_rho = FC_Voxels.run(iRun).rho_RestingState;
            
            fisher(iRun,:,:) = (1/2).*log((1 + single_rho)./(1 - single_rho));
            
        end
        
        [l,c,z] = size(fisher);
        
        fisher_pairs = reshape(fisher,[l c*z]);
        
        mu_z = mean(fisher_pairs,1);
        sigma_z = std(fisher_pairs);
   
        biInfo(iSubject,iROI) = TotalInformationContent( mu_z, sigma_z, z_threshold, cv_threshold );
            
        clear FC_Voxels subject file area_label fisher fisher_pairs mu_z sigma_z
        
    end
    
end

save('PETRA-biInfo.mat','biInfo');

end

function getInformationRunAALVoxelLevelHCP


nROIs = 90;
nSubjects = 8;
nRuns = 4;
nTotalRuns = 32;

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

biInfo = zeros(nSubjects,nROIs);
for iROI=1:nROIs
   
    area_Label = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    
    disp(strcat(area_Label,'-',datestr(now)));
    
    iiRun = 0;
    for iSubject=1:nSubjects

        for iRun=1:nRuns
            
            load(strcat('HCP','-','SUBJ',int2str(iSubject),'-','FC-Voxels-AAL-ROI-corr-HCP-RestingState-',area_Label,'-',int2str(iRun),'.mat'));
            
            iiRun = iiRun + 1;
            
            single_rho = FC_Voxels.run(iRun).rho_RestingState;
            
            fisher(iRun,:,:) = (1/2).*log((1 + single_rho)./(1 - single_rho));
            
        end
        
        [l,c,z] = size(fisher);
        
        fisher_pairs = reshape(fisher,[l c*z]);
        
        mu_z = mean(fisher_pairs,1);
        sigma_z = std(fisher_pairs);
   
        biInfo(iSubject,iROI) = TotalInformationContent( mu_z, sigma_z, z_threshold, cv_threshold );
            
        clear FC_Voxels subject file area_label fisher fisher_pairs mu_z sigma_z
        
    end
    
end

save('HCP-biInfo.mat','biInfo');

end

function plotInfoLHRPETRA

nROI = 90;

LHR = load('LHR-biInfo.mat');
PETRA = load('PETRA-biInfo.mat');

biInfo = [LHR.biInfo;PETRA.biInfo];

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for iROI=1:nROI
    
   idx_ROI  = AAL_ROI(iROI).ID; 
   idx_voxels = find(AAL_img == idx_ROI);
   
   nVoxels = length(idx_voxels);
   
   %biInfo(:,iROI) = biInfo(:,iROI) ./ ( nVoxels * ( nVoxels - 1 ) / 2 );
   
   LHR.biInfo(:,iROI) = LHR.biInfo(:,iROI) ./ ( nVoxels * ( nVoxels - 1 ) / 2 );
   PETRA.biInfo(:,iROI) = PETRA.biInfo(:,iROI) ./ ( nVoxels * ( nVoxels - 1 ) / 2 );
   
   
end

[H P] = ttest2(LHR.biInfo,PETRA.biInfo);

m_LHR = mean(LHR.biInfo,1);
m_PETRA = mean(PETRA.biInfo,1);

most_LHR = (m_LHR > m_PETRA) & H;
most_PETRA = (m_PETRA > m_LHR) & H;

for iROI=1:nROI
    
    if most_LHR(iROI)
        
        disp(strrep(AAL_ROI(iROI).Nom_L,'_','-'));
        
    end
    
end

imagesc(biInfo);

hold on
plot([0 90],[8.5 8.5],'w-');

for i=1:90
    my_labels{i} = AAL_ROI(i).Nom_L;
end

set(gca,'XTick',1:90);
set(gca,'XTickLabel',my_labels);

xticklabel_rotate;

colorbar;

end

function plotInfoLHRHCP

nROI = 90;

LHR = load('LHR-biInfo.mat');
HCP = load('HCP-biInfo.mat');

biInfo = [LHR.biInfo;HCP.biInfo];

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

for iROI=1:nROI
    
   idx_ROI  = AAL_ROI(iROI).ID; 
   idx_voxels = find(AAL_img == idx_ROI);
   
   nVoxels = length(idx_voxels);
   
   %biInfo(:,iROI) = biInfo(:,iROI) ./ ( nVoxels * ( nVoxels - 1 ) / 2 );
   
   LHR.biInfo(:,iROI) = LHR.biInfo(:,iROI) ./ ( nVoxels * ( nVoxels - 1 ) / 2 );
   HCP.biInfo(:,iROI) = HCP.biInfo(:,iROI) ./ ( nVoxels * ( nVoxels - 1 ) / 2 );
   
   
end

[H P] = ttest2(LHR.biInfo,HCP.biInfo);

m_LHR = mean(LHR.biInfo,1);
m_HCP = mean(HCP.biInfo,1);

most_LHR = (m_LHR > m_HCP) & H;
most_HCP = (m_HCP > m_LHR) & H;

for iROI=1:nROI
    
    if most_LHR(iROI)
        
        disp(strrep(AAL_ROI(iROI).Nom_L,'_','-'));
        
    end
    
end

imagesc(biInfo);

hold on
plot([0 90],[8.5 8.5],'w-');

for i=1:90
    my_labels{i} = AAL_ROI(i).Nom_L;
end

set(gca,'XTick',1:90);
set(gca,'XTickLabel',my_labels);

xticklabel_rotate;

colorbar;

end


%%% FUNCTIONAL CLUSTER

function getCorrClustersDifferentClusterSameDifferentAAL

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nClusters = 758;
nROI = 90;
nRuns = 32;

for iRun=1:nRuns
    
    iiCluster = 0;
    
    for iROI=1:nROI
    
        nClusters = length(ROI(iROI).clusters);
        
        for iCluster=1:nClusters
      
            iiCluster = iiCluster + 1;
            
            func_clusters(iRun).rest_run(iiCluster,:) = ROI(iROI).clusters(iCluster).rest_mean(iRun,:);
            func_clusters(iRun).passive_run(iiCluster,:) = ROI(iROI).clusters(iCluster).passive_mean(iRun,:);
            func_clusters(iRun).track_run(iiCluster,:) = ROI(iROI).clusters(iCluster).track_mean(iRun,:);
            
        end
        
    end
    
end

for iRun=1:nRuns
    
   func_clusters(iRun).rest_corr = corr(func_clusters(iRun).rest_run');
   func_clusters(iRun).passive_corr = corr(func_clusters(iRun).passive_run');
   func_clusters(iRun).track_corr = corr(func_clusters(iRun).track_run');
    
end
       
save('758-Functional-Clusters-Corr.mat','func_clusters');

end

function getCorrClustersPARCELSpetraDATAmd

load('Z:\Dropbox (Uni Magdeburg)\_DATA\PETRA-RITTER\all_subjects\FC-Voxels-AAL-Kmeans-Volume\PETRA-FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nTotalClusters = 758;
nROI = 90;
nTotalRuns = 32;
nRuns = 4;
nTR = 150;
MNI_size = [91 109 91];

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

for iRun=1:nTotalRuns
    
    disp(strcat('iRun:',int2str(iRun)));
    
    func_clusters(iRun).corr = zeros(nTotalClusters);
    
    this_run_clusters = zeros(nTR,nTotalClusters);
    
    iiCluster = 0;
    
    for iROI=1:nROI
    
        nClusters = length(ROI(iROI).clusters);
        
        for iCluster=1:nClusters
      
            iiCluster = iiCluster + 1;
            
            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
            
            nVoxels = length(idx_voxels);
            
            voxels_this_cluster = zeros(nTR,nVoxels);
            
            for iVoxel=1:nVoxels
                
                [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
                
                voxels_this_cluster(:,iVoxel) = RestingState(iRun).run(idxx,idxy,idxz,1:nTR);
                
            end
            
            m_voxels_this_cluster = mean(voxels_this_cluster,2);
            
            this_run_clusters(:,iiCluster) = m_voxels_this_cluster(:);
            
        end
        
    end
    
    func_clusters(iRun).corr = corr(this_run_clusters);
    
end
       
save('758-Functional-Clusters-Corr-PARCELS-PETRA-DATA-MD.mat','func_clusters');

end

function getCorrClustersDifferentClusterEverybodyAALPETRA

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nClusters = 758;
nROI = 90;
nRuns = 9*3;
MNI_size = [91 109 91];
nTR = 150;
nSegments = 3;

run = 1;
all_settings = getAllSettingsPetra;

iRun = 0;
for iSet=1:length(all_settings)
    
    settings = all_settings(iSet).settings;

    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;

    [residual, settings] = real_get_data_FSL_Petra(settings,run,file,get_at_this_preprocessed_step);
    
    for iSeg=1:nSegments
        
        iRun = iRun + 1;
        
        Subject(iRun).residual(:,:,:,:) = residual(:,:,:,(1+(iSeg-1)*nTR):(nTR+(iSeg-1)*nTR));
        
    end

end

for iRun=1:nRuns
               
    func_clusters(iRun).rest_run = zeros(nTR,nClusters);
    
    iiCluster = 0;
    for iROI=1:nROI

        nCluster = length(ROI(iROI).clusters);

        for iCluster=1:nCluster
            
            iiCluster = iiCluster + 1;

            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);
            
            this_cluster_voxels = zeros(nVoxels,nTR);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

                this_cluster_voxels(iVoxel,:) = squeeze(Subject(iRun).residual(idxx,idxy,idxz,:));

            end
            
            func_clusters(iRun).rest_run(:,iiCluster) = mean(this_cluster_voxels,1);

        end

    end

end


for iRun=1:nRuns
    
   func_clusters(iRun).rest_corr = corr(func_clusters(iRun).rest_run);
    
end
       
save('758-Functional-Clusters-Petra-Corr.mat','func_clusters','-v7.3');

end

function getCorrClustersDifferentClusterEverybodyAALPETRACustomParcellation

load('PETRA-FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nClusters = 758;
nROI = 90;
nRuns = 9*3;
MNI_size = [91 109 91];
nTR = 150;
nSegments = 3;

run = 1;
all_settings = getAllSettingsPetra;

iRun = 0;
for iSet=1:length(all_settings)
    
    settings = all_settings(iSet).settings;

    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;

    [residual, settings] = real_get_data_FSL_Petra(settings,run,file,get_at_this_preprocessed_step);
    
    for iSeg=1:nSegments
        
        iRun = iRun + 1;
        
        Subject(iRun).residual(:,:,:,:) = residual(:,:,:,(1+(iSeg-1)*nTR):(nTR+(iSeg-1)*nTR));
        
    end

end

for iRun=1:nRuns
               
    func_clusters(iRun).rest_run = zeros(nTR,nClusters);
    
    iiCluster = 0;
    for iROI=1:nROI

        nCluster = length(ROI(iROI).clusters);

        for iCluster=1:nCluster
            
            iiCluster = iiCluster + 1;

            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);
            
            this_cluster_voxels = zeros(nVoxels,nTR);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

                this_cluster_voxels(iVoxel,:) = squeeze(Subject(iRun).residual(idxx,idxy,idxz,:));

            end
            
            func_clusters(iRun).rest_run(:,iiCluster) = mean(this_cluster_voxels,1);

        end

    end

end


for iRun=1:nRuns
    
   func_clusters(iRun).rest_corr = corr(func_clusters(iRun).rest_run);
    
end
       
save('758-Functional-Clusters-Petra-Custom-Corr.mat','func_clusters','-v7.3');

end

function getCorrClustersDifferentClusterEverybodyAALHCP

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nClusters = 758;
nROI = 90;
nTotalRuns = 32;
nRuns = 4;
MNI_size = [91 109 91];
nTR = 1200;

all_settings = getAllSettingsHCP;
nSubjects = length(all_settings);

iiRun = 0;
for iSet=1:5
    
    settings = all_settings(iSet).settings;
    
    disp(strcat('Subj:',int2str(iSet)));

    for iRun=1:nRuns
        
        iiRun = iiRun + 1;
        
        [all_runs(iiRun).run, settings] = real_get_data_HCP(settings,iRun);
    
    end
    
end

i_global_Run = 0;

for iRun=1:20
    
    i_global_Run = i_global_Run + 1;
               
    func_clusters(i_global_Run).rest_run = zeros(nTR,nClusters);
    
    iiCluster = 0;
    for iROI=1:nROI

        nCluster = length(ROI(iROI).clusters);

        for iCluster=1:nCluster
            
            iiCluster = iiCluster + 1;

            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);
            
            this_cluster_voxels = zeros(nVoxels,nTR);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

                this_cluster_voxels(iVoxel,:) = squeeze(all_runs(iRun).run(idxx,idxy,idxz,:));

            end
            
            func_clusters(i_global_Run).rest_run(:,iiCluster) = mean(this_cluster_voxels,1);

        end

    end

end

clear all_runs

iiRun = 0;
for iSet=6:8
    
    settings = all_settings(iSet).settings;
    
    disp(strcat('Subj:',int2str(iSet)));

    for iRun=1:nRuns
        
        iiRun = iiRun + 1;
        
        [all_runs(iiRun).run, settings] = real_get_data_HCP(settings,iRun);
    
    end
    
end

for iRun=1:12
    
    i_global_Run = i_global_Run + 1;
               
    func_clusters(i_global_Run).rest_run = zeros(nTR,nClusters);
    
    iiCluster = 0;
    for iROI=1:nROI

        nCluster = length(ROI(iROI).clusters);

        for iCluster=1:nCluster
            
            iiCluster = iiCluster + 1;

            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);
            
            this_cluster_voxels = zeros(nVoxels,nTR);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

                this_cluster_voxels(iVoxel,:) = squeeze(all_runs(iRun).run(idxx,idxy,idxz,:));

            end
            
            func_clusters(i_global_Run).rest_run(:,iiCluster) = mean(this_cluster_voxels,1);

        end

    end

end

for iRun=1:nTotalRuns
    
   func_clusters(iRun).rest_corr = corr(func_clusters(iRun).rest_run);
    
end
       
save('758-Functional-Clusters-HCP-Corr.mat','func_clusters','-v7.3');

end

function getCorrClustersDifferentClusterEverybodyAALPETRAWholeSegment

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nClusters = 758;
nROI = 90;
nRuns = 9;
MNI_size = [91 109 91];
nTR = 450;

run = 1;
all_settings = getAllSettingsPetra;

iRun = 0;
for iSet=1:length(all_settings)
    
    settings = all_settings(iSet).settings;

    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;

    [residual, settings] = real_get_data_FSL_Petra(settings,run,file,get_at_this_preprocessed_step);

    iRun = iRun + 1;

    Subject(iRun).residual(:,:,:,:) = residual(:,:,:,:);
 
end

for iRun=1:nRuns
               
    func_clusters(iRun).rest_run = zeros(nTR,nClusters);
    
    iiCluster = 0;
    for iROI=1:nROI

        nCluster = length(ROI(iROI).clusters);

        for iCluster=1:nCluster
            
            iiCluster = iiCluster + 1;

            idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);
            
            this_cluster_voxels = zeros(nVoxels,nTR);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

                this_cluster_voxels(iVoxel,:) = squeeze(Subject(iRun).residual(idxx,idxy,idxz,:));

            end
            
            func_clusters(iRun).rest_run(:,iiCluster) = mean(this_cluster_voxels,1);

        end

    end

end


for iRun=1:nRuns
    
   func_clusters(iRun).rest_corr = corr(func_clusters(iRun).rest_run);
    
end
       
save('758-Functional-Clusters-Petra-Whole-Segment-Corr.mat','func_clusters','-v7.3');

end

function plotVarianceDistributionClustersDifferentClusterSameDifferentAAL

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Corr.mat');

nROI = 90;
nRuns = 32;

%%% SAME ALL - REST

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        for iiCluster=1:nClusters
            
            iPair = iPair + 1;
        
            for iRun=1:nRuns

                single_rho = func_clusters(iRun).rest_corr(iiiCluster,iiiCluster-iCluster+iiCluster);
                
                pairs_inside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

            end
        
        end

    end
    
end

my_var_inside = var(pairs_inside);

my_cv_inside = std(pairs_inside) ./ abs(mean(pairs_inside,1));

my_mean_inside = mean(pairs_inside,1);
my_std_inside = std(pairs_inside);

vector_one = my_std_inside;
vector_two = my_mean_inside;
vector_three = my_cv_inside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Same-AAL-Rest';

MyContourCorrelation_v2(vector_two,vector_one,title_label,1);

nROI = 90;
nRuns = 32;

%%% SAME ALL - PASSIVE

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        for iiCluster=1:nClusters
            
            iPair = iPair + 1;
        
            for iRun=1:nRuns

                single_rho = func_clusters(iRun).passive_corr(iiiCluster,iiiCluster-iCluster+iiCluster);
                
                pairs_inside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

            end
        
        end

    end
    
end

my_var_inside = var(pairs_inside);

my_cv_inside = std(pairs_inside) ./ abs(mean(pairs_inside,1));

my_mean_inside = mean(pairs_inside,1);
my_std_inside = std(pairs_inside);

vector_one = my_std_inside;
vector_two = my_mean_inside;
vector_three = my_cv_inside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Same-AAL-Passive';

MyContourCorrelation_v2(vector_two,vector_one,title_label,1);

nROI = 90;
nRuns = 32;

%%% SAME ALL - TRACK

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        for iiCluster=1:nClusters
            
            iPair = iPair + 1;
        
            for iRun=1:nRuns

                single_rho = func_clusters(iRun).track_corr(iiiCluster,iiiCluster-iCluster+iiCluster);
                
                pairs_inside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

            end
        
        end

    end
    
end

my_var_inside = var(pairs_inside);

my_cv_inside = std(pairs_inside) ./ abs(mean(pairs_inside,1));

my_mean_inside = mean(pairs_inside,1);
my_std_inside = std(pairs_inside);

vector_one = my_std_inside;
vector_two = my_mean_inside;
vector_three = my_cv_inside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Same-AAL-Track';

MyContourCorrelation_v2(vector_two,vector_one,title_label,1);

%%% DIFFERENT AAL - REST

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        ivCluster = 0;
        for iiROI=1:nROI
            
            nnClusters = length(ROI(iiROI).clusters);
            
            for iiCluster=1:nnClusters
                
                ivCluster = ivCluster + 1;
                
                if iROI ~= iiROI
                    
                    iPair = iPair + 1;
                    
                    for iRun=1:nRuns

                        single_rho = func_clusters(iRun).rest_corr(iiiCluster,ivCluster);
                        
                        pairs_outside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

my_var_outside = var(pairs_outside);

my_cv_outside = std(pairs_outside) ./ abs(mean(pairs_outside,1));

my_mean_outside = mean(pairs_outside,1);
my_std_outside = std(pairs_outside);

vector_one = my_std_outside;
vector_two = my_mean_outside;
vector_three = my_cv_outside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Different-AAL-Rest';

MyContourCorrelation_v2(vector_two,vector_one,title_label,1);

%%% DIFFERENT AAL - PASSIVE

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        ivCluster = 0;
        for iiROI=1:nROI
            
            nnClusters = length(ROI(iiROI).clusters);
            
            for iiCluster=1:nnClusters
                
                ivCluster = ivCluster + 1;
                
                if iROI ~= iiROI
                    
                    iPair = iPair + 1;
                    
                    for iRun=1:nRuns
                        
                        single_rho = func_clusters(iRun).passive_corr(iiiCluster,ivCluster);
                        
                        pairs_outside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

my_var_outside = var(pairs_outside);

my_cv_outside = std(pairs_outside) ./ abs(mean(pairs_outside,1));

my_mean_outside = mean(pairs_outside,1);
my_std_outside = std(pairs_outside);

vector_one = my_std_outside;
vector_two = my_mean_outside;
vector_three = my_cv_outside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Different-AAL-Passive';

MyContourCorrelation_v2(vector_two,vector_one,title_label,1);

%%% DIFFERENT AAL- TRACK

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        ivCluster = 0;
        for iiROI=1:nROI
            
            nnClusters = length(ROI(iiROI).clusters);
            
            for iiCluster=1:nnClusters
                
                ivCluster = ivCluster + 1;
                
                if iROI ~= iiROI
                    
                    iPair = iPair + 1;
                    
                    for iRun=1:nRuns
                        
                        single_rho = func_clusters(iRun).track_corr(iiiCluster,ivCluster);
                        
                        pairs_outside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

my_var_outside = var(pairs_outside);

my_cv_outside = std(pairs_outside) ./ abs(mean(pairs_outside,1));

my_mean_outside = mean(pairs_outside,1);
my_std_outside = std(pairs_outside);

vector_one = my_std_outside;
vector_two = my_mean_outside;
vector_three = my_cv_outside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Different-AAL-Track';

MyContourCorrelation_v2(vector_two,vector_one,title_label,1);

end

function plotVarianceDistributionClustersDifferentClusterEverybodyAALPETRA

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Petra-Corr.mat');

nROI = 90;
% nRuns = 27;
nRuns = 24;
nClusters = 758;

iPair = 0;
for iCluster=1:nClusters

    for iiCluster=1:nClusters

        iPair = iPair + 1;

        for iRun=1:nRuns

            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end

end

my_var = nanvar(pairs);

my_cv = nanstd(pairs) ./ abs(nanmean(pairs,1));

my_mean = nanmean(pairs,1);
my_std = nanstd(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Petra-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotVarCluPARCELSpetraDATAmd

load('758-Functional-Clusters-Corr-PARCELS-PETRA-DATA-MD.mat');

nTotalRuns = 32;
nTotalClusters = 758;

iPair = 0;
for iCluster=1:nTotalClusters

    for iiCluster=1:nTotalClusters

        iPair = iPair + 1;

        for iRun=1:nTotalRuns

            single_rho = func_clusters(iRun).corr(iCluster,iiCluster);

            pairs(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end

end

my_var = nanvar(pairs);

my_cv = nanstd(pairs) ./ abs(nanmean(pairs,1));

my_mean = nanmean(pairs,1);
my_std = nanstd(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Func-Parcel-Petra-Data-MD';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotFCCluPARCELSpetraDATAmd

load('758-Functional-Clusters-Corr-PARCELS-PETRA-DATA-MD.mat');

nTotalRuns = 32;
nTotalClusters = 758;

iPair = 0;
for iCluster=1:nTotalClusters

    for iiCluster=1:nTotalClusters

        iPair = iPair + 1;

        for iRun=1:nTotalRuns

            single_rho = func_clusters(iRun).corr(iCluster,iiCluster);

            pairs(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end

end

my_var = nanvar(pairs);

my_cv = nanstd(pairs) ./ abs(nanmean(pairs,1));

my_mean = nanmean(pairs,1);
my_std = nanstd(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Func-Parcel-Petra-Data-MD';

iPair = 0;
jCluster = 0;
for iCluster=1:nTotalClusters
    
    jCluster = jCluster + 1;
    
    jjCluster = 0;
    for iiCluster=1:nTotalClusters
        
        jjCluster = jjCluster + 1;
        
            iPair = iPair + 1;
            
            if iCluster ~= iiCluster
        
                %if my_cv_aal(iPair) < cv_criterion & p_pairs(iPair) < p_criterion

                    mean_corr(jCluster,jjCluster) =  my_mean(iPair);
                    s_corr(jCluster,jjCluster) =  my_std(iPair);

                %end
            
            end
            
    end
    
end

plotFCgeneral(mean_corr,s_corr,title_label);

end

function plotVarCluEveryAALPETRACustomParcellation

load('PETRA-FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');
load('758-Functional-Clusters-Petra-Custom-Corr.mat');

nROI = 90;
nRuns = 27;
nClusters = 758;

iPair = 0;
for iCluster=1:nClusters

    for iiCluster=1:nClusters

        iPair = iPair + 1;

        for iRun=1:nRuns

            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end

end

my_var = nanvar(pairs);

my_cv = nanstd(pairs) ./ abs(nanmean(pairs,1));

my_mean = nanmean(pairs,1);
my_std = nanstd(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Petra-Custom-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotVarianceDistributionClustersDifferentClusterEverybodyAALHCP

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-HCP-Corr.mat');

nROI = 90;
nRuns = 32;
nClusters = 758;

iPair = 0;
for iCluster=1:nClusters

    for iiCluster=1:nClusters

        iPair = iPair + 1;

        for iRun=1:nRuns

            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end

end

my_var = nanvar(pairs);

my_cv = nanstd(pairs) ./ abs(nanmean(pairs,1));

my_mean = nanmean(pairs,1);
my_std = nanstd(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-HCP-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotVarianceDistributionClustersDifferentClusterEverybodyAALPETRAWholeSegment

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Petra-Whole-Segment-Corr.mat');

nROI = 90;
nRuns = 9;
nClusters = 758;

iPair = 0;
for iCluster=1:nClusters

    for iiCluster=1:nClusters

        iPair = iPair + 1;

        for iRun=1:nRuns

            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end

end

my_var = nanvar(pairs);

my_cv = nanstd(pairs) ./ abs(nanmean(pairs,1));

my_mean = nanmean(pairs,1);
my_std = nanstd(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Petra-Whole-Segment-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotVarianceDistributionClustersDifferentClusterEverybody

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Corr.mat');

nTotalCluster = 758;
nRuns = 32;

%%% EVERYBODY - REST

iPair = 0;
for iCluster=1:nTotalCluster
   
    for iiCluster=1:nTotalCluster

        iPair = iPair + 1;

        for iRun=1:nRuns

            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end
    
end

my_var = var(pairs);

my_cv = std(pairs) ./ abs(mean(pairs,1));

my_mean = mean(pairs,1);
my_std = std(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Functional-Everybody-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

fraction = FractionSignificantConsistent( vector_two, vector_one, z_threshold, cv_threshold )

end

function plotVarianceDistributionClustersDifferentClusterEverybodyONESubject(nStartRun,nEndRun)

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Corr.mat');

nTotalCluster = 758;
nRuns = 4;

%%% EVERYBODY - REST

iPair = 0;
for iCluster=1:nTotalCluster
   
    for iiCluster=1:nTotalCluster

        iPair = iPair + 1;

        iiRun = 0;
        for iRun=nStartRun:nEndRun
            
            iiRun = iiRun + 1;

            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs(iiRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

        end

    end
    
end

my_var = var(pairs);

my_cv = std(pairs) ./ abs(mean(pairs,1));

my_mean = mean(pairs,1);
my_std = std(pairs);

vector_one = my_std;
vector_two = my_mean;
vector_three = my_cv;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = strcat('Functional-Everybody-Rest','-Runs-',int2str(nStartRun),'-',int2str(nEndRun));

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

fraction = FractionSignificantConsistent( vector_two, vector_one, z_threshold, cv_threshold )

end

function plotMeanFCFunctionalClusters(nStartCluster,nEndCluster,label)

showAAL = 0;

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Corr.mat');

nROI = 90;
nRuns = 32;
nTotalClusters = 758;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

%%% ALL PAIRS - REST

iPair = 0;
for iCluster=nStartCluster:nEndCluster
    
    for iiCluster=nStartCluster:nEndCluster
        
        iPair = iPair + 1;
        
        for iRun=1:nRuns
            
            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            pairs_rest_rho(iRun,iPair) = single_rho;
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ abs(mean(pairs_rest_z,1));

% mean_corr = zeros(nEndCluster-nStartCluster + 1,nEndCluster-nStartCluster + 1);

iPair = 0;
jCluster = 0;
for iCluster=nStartCluster:nEndCluster
    
    jCluster = jCluster + 1;
    
    jjCluster = 0;
    for iiCluster=nStartCluster:nEndCluster
        
        jjCluster = jjCluster + 1;
        
            iPair = iPair + 1;
            
            if iCluster ~= iiCluster
        
                %if my_cv_aal(iPair) < cv_criterion & p_pairs(iPair) < p_criterion

                    mean_corr(jCluster,jjCluster) =  m_pairs_z(iPair);
                    s_corr(jCluster,jjCluster) =  s_pairs_z(iPair);

                %end
            
            end
            
    end
    
end

% f = figure;
% 
% hold on
% 
% clmap = colormap('jet');
% clmap(1,:) = [1 1 1];
% 
% min_C = 0;
% max_C = 1;
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% if showAAL
%     
%     before_Angular_L = 545;
%     after_Angular_L = 6;
%     
%     before_Frontal_Mid_L = 78;
%     after_Frontal_Mid_L = 25;
%     
%     rectangle('Position',[before_Angular_L before_Angular_L after_Angular_L after_Angular_L],'Curvature',0.2);
%     rectangle('Position',[before_Frontal_Mid_L before_Frontal_Mid_L after_Frontal_Mid_L after_Frontal_Mid_L],'Curvature',0.2);
% 
% end
% 
% print(f,'-depsc',strcat('FC-Functional-Clusters','-',label,'.eps'));

new_label = strcat('FC-Functional-Clusters','-',label);
plotFCgeneral(mean_corr,s_corr,new_label);


end

function plotMeanFCFunctionalClustersONESubject(nStartRun,nEndRun,nStartCluster,nEndCluster,label)

showAAL = 0;

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Corr.mat');

nROI = 90;
nRuns = 32;
nTotalClusters = 758;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

%%% ALL PAIRS - REST

iPair = 0;
for iCluster=nStartCluster:nEndCluster
    
    for iiCluster=nStartCluster:nEndCluster
        
        iPair = iPair + 1;
        
        iiRun = 0;
        for iRun=nStartRun:nEndRun
            
            iiRun = iiRun + 1;
            
            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs_rest_z(iiRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            pairs_rest_rho(iiRun,iPair) = single_rho;
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ abs(mean(pairs_rest_z,1));

% mean_corr = zeros(nEndCluster-nStartCluster + 1,nEndCluster-nStartCluster + 1);

iPair = 0;
jCluster = 0;
for iCluster=nStartCluster:nEndCluster
    
    jCluster = jCluster + 1;
    
    jjCluster = 0;
    for iiCluster=nStartCluster:nEndCluster
        
        jjCluster = jjCluster + 1;
        
            iPair = iPair + 1;
            
            if iCluster ~= iiCluster
        
                %if my_cv_aal(iPair) < cv_criterion & p_pairs(iPair) < p_criterion

                    mean_corr(jCluster,jjCluster) =  m_pairs_z(iPair);
                    s_corr(jCluster,jjCluster) =  s_pairs_z(iPair);

                %end
            
            end
            
    end
    
end

% f = figure;
% 
% hold on
% 
% clmap = colormap('jet');
% clmap(1,:) = [1 1 1];
% 
% min_C = 0;
% max_C = 1;
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% if showAAL
%     
%     before_Angular_L = 545;
%     after_Angular_L = 6;
%     
%     before_Frontal_Mid_L = 78;
%     after_Frontal_Mid_L = 25;
%     
%     rectangle('Position',[before_Angular_L before_Angular_L after_Angular_L after_Angular_L],'Curvature',0.2);
%     rectangle('Position',[before_Frontal_Mid_L before_Frontal_Mid_L after_Frontal_Mid_L after_Frontal_Mid_L],'Curvature',0.2);
% 
% end
% 
% print(f,'-depsc',strcat('FC-Functional-Clusters','-',label,'.eps'));

new_label = strcat('FC-Functional-Clusters','-',label,'-Runs-',int2str(nStartRun),'-',int2str(nEndRun));
plotFCgeneral(mean_corr,s_corr,new_label);


end

function plotMeanFCFunctionalClustersPETRA(nStartCluster,nEndCluster,label)

showAAL = 0;

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Petra-Corr.mat');

nROI = 90;
nRuns = 27;
nTotalClusters = 758;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

%%% ALL PAIRS - REST

iPair = 0;
for iCluster=nStartCluster:nEndCluster
    
    for iiCluster=nStartCluster:nEndCluster
        
        iPair = iPair + 1;
        
        for iRun=1:nRuns
            
            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            pairs_rest_rho(iRun,iPair) = single_rho;
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ abs(mean(pairs_rest_z,1));

% mean_corr = zeros(nEndCluster-nStartCluster + 1,nEndCluster-nStartCluster + 1);

iPair = 0;
jCluster = 0;
for iCluster=nStartCluster:nEndCluster
    
    jCluster = jCluster + 1;
    
    jjCluster = 0;
    for iiCluster=nStartCluster:nEndCluster
        
        jjCluster = jjCluster + 1;
        
            iPair = iPair + 1;
            
            if iCluster ~= iiCluster
        
                %if my_cv_aal(iPair) < cv_criterion & p_pairs(iPair) < p_criterion

                    mean_corr(jCluster,jjCluster) =  m_pairs_z(iPair);
                    s_corr(jCluster,jjCluster) =  s_pairs_z(iPair);

                %end
            
            end
            
    end
    
end

% f = figure;
% 
% hold on
% 
% clmap = colormap('jet');
% clmap(1,:) = [1 1 1];
% 
% min_C = 0;
% max_C = 1;
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% if showAAL
%     
%     before_Angular_L = 545;
%     after_Angular_L = 6;
%     
%     before_Frontal_Mid_L = 78;
%     after_Frontal_Mid_L = 25;
%     
%     rectangle('Position',[before_Angular_L before_Angular_L after_Angular_L after_Angular_L],'Curvature',0.2);
%     rectangle('Position',[before_Frontal_Mid_L before_Frontal_Mid_L after_Frontal_Mid_L after_Frontal_Mid_L],'Curvature',0.2);
% 
% end
% 
% print(f,'-depsc',strcat('FC-Functional-Clusters','-',label,'.eps'));

new_label = strcat('FC-Functional-Clusters-Petra','-',label);
plotFCgeneral(mean_corr,s_corr,new_label);


end

function plotMeanFCFunctionalClustersPETRAcustom(nStartCluster,nEndCluster,label)

showAAL = 0;

load('PETRA-FC-Voxels-AAL-ROI-corr-custom-KMeans-Info.mat');
load('758-Functional-Clusters-Petra-Custom-Corr.mat');

nROI = 90;
nRuns = 27;
nTotalClusters = 758;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

%%% ALL PAIRS - REST

iPair = 0;
for iCluster=nStartCluster:nEndCluster
    
    for iiCluster=nStartCluster:nEndCluster
        
        iPair = iPair + 1;
        
        for iRun=1:nRuns
            
            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            pairs_rest_rho(iRun,iPair) = single_rho;
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ abs(mean(pairs_rest_z,1));

% mean_corr = zeros(nEndCluster-nStartCluster + 1,nEndCluster-nStartCluster + 1);

iPair = 0;
jCluster = 0;
for iCluster=nStartCluster:nEndCluster
    
    jCluster = jCluster + 1;
    
    jjCluster = 0;
    for iiCluster=nStartCluster:nEndCluster
        
        jjCluster = jjCluster + 1;
        
            iPair = iPair + 1;
            
            if iCluster ~= iiCluster
        
                %if my_cv_aal(iPair) < cv_criterion & p_pairs(iPair) < p_criterion

                    mean_corr(jCluster,jjCluster) =  m_pairs_z(iPair);
                    s_corr(jCluster,jjCluster) =  s_pairs_z(iPair);

                %end
            
            end
            
    end
    
end

% f = figure;
% 
% hold on
% 
% clmap = colormap('jet');
% clmap(1,:) = [1 1 1];
% 
% min_C = 0;
% max_C = 1;
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% if showAAL
%     
%     before_Angular_L = 545;
%     after_Angular_L = 6;
%     
%     before_Frontal_Mid_L = 78;
%     after_Frontal_Mid_L = 25;
%     
%     rectangle('Position',[before_Angular_L before_Angular_L after_Angular_L after_Angular_L],'Curvature',0.2);
%     rectangle('Position',[before_Frontal_Mid_L before_Frontal_Mid_L after_Frontal_Mid_L after_Frontal_Mid_L],'Curvature',0.2);
% 
% end
% 
% print(f,'-depsc',strcat('FC-Functional-Clusters','-',label,'.eps'));

new_label = strcat('FC-Functional-Clusters-Petra-custom','-',label);
plotFCgeneral(mean_corr,s_corr,new_label);


end

function plotMeanFCFunctionalClustersHCP(nStartCluster,nEndCluster,label)

showAAL = 0;

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-HCP-Corr.mat');

nROI = 90;
nRuns = 32;
nTotalClusters = 758;
nTR = 1200;
cv_criterion = 0.75;
p_criterion = 0.01;

%%% ALL PAIRS - REST

iPair = 0;
for iCluster=nStartCluster:nEndCluster
    
    for iiCluster=nStartCluster:nEndCluster
        
        iPair = iPair + 1;
        
        for iRun=1:nRuns
            
            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            pairs_rest_rho(iRun,iPair) = single_rho;
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ abs(mean(pairs_rest_z,1));

% mean_corr = zeros(nEndCluster-nStartCluster + 1,nEndCluster-nStartCluster + 1);

iPair = 0;
jCluster = 0;
for iCluster=nStartCluster:nEndCluster
    
    jCluster = jCluster + 1;
    
    jjCluster = 0;
    for iiCluster=nStartCluster:nEndCluster
        
        jjCluster = jjCluster + 1;
        
            iPair = iPair + 1;
            
            if iCluster ~= iiCluster
        
                %if my_cv_aal(iPair) < cv_criterion & p_pairs(iPair) < p_criterion

                    mean_corr(jCluster,jjCluster) =  m_pairs_z(iPair);
                    s_corr(jCluster,jjCluster) =  s_pairs_z(iPair);

                %end
            
            end
            
    end
    
end

% f = figure;
% 
% hold on
% 
% clmap = colormap('jet');
% clmap(1,:) = [1 1 1];
% 
% min_C = 0;
% max_C = 1;
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% if showAAL
%     
%     before_Angular_L = 545;
%     after_Angular_L = 6;
%     
%     before_Frontal_Mid_L = 78;
%     after_Frontal_Mid_L = 25;
%     
%     rectangle('Position',[before_Angular_L before_Angular_L after_Angular_L after_Angular_L],'Curvature',0.2);
%     rectangle('Position',[before_Frontal_Mid_L before_Frontal_Mid_L after_Frontal_Mid_L after_Frontal_Mid_L],'Curvature',0.2);
% 
% end
% 
% print(f,'-depsc',strcat('FC-Functional-Clusters','-',label,'.eps'));

new_label = strcat('FC-Functional-Clusters-HCP','-',label);
plotFCgeneral(mean_corr,s_corr,new_label);


end

function plotMeanFCFunctionalClustersPETRAWholeSegment(nStartCluster,nEndCluster,label)

showAAL = 0;

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('758-Functional-Clusters-Petra-Whole-Segment-Corr.mat');

nROI = 90;
nRuns = 9;
nTotalClusters = 758;
nTR = 450;
cv_criterion = 0.75;
p_criterion = 0.01;

%%% ALL PAIRS - REST

iPair = 0;
for iCluster=nStartCluster:nEndCluster
    
    for iiCluster=nStartCluster:nEndCluster
        
        iPair = iPair + 1;
        
        for iRun=1:nRuns
            
            single_rho = func_clusters(iRun).rest_corr(iCluster,iiCluster);

            pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            pairs_rest_rho(iRun,iPair) = single_rho;
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ abs(mean(pairs_rest_z,1));

% mean_corr = zeros(nEndCluster-nStartCluster + 1,nEndCluster-nStartCluster + 1);

iPair = 0;
jCluster = 0;
for iCluster=nStartCluster:nEndCluster
    
    jCluster = jCluster + 1;
    
    jjCluster = 0;
    for iiCluster=nStartCluster:nEndCluster
        
        jjCluster = jjCluster + 1;
        
            iPair = iPair + 1;
            
            if iCluster ~= iiCluster
        
                %if my_cv_aal(iPair) < cv_criterion & p_pairs(iPair) < p_criterion

                    mean_corr(jCluster,jjCluster) =  m_pairs_z(iPair);
                    s_corr(jCluster,jjCluster) =  s_pairs_z(iPair);

                %end
            
            end
            
    end
    
end

% f = figure;
% 
% hold on
% 
% clmap = colormap('jet');
% clmap(1,:) = [1 1 1];
% 
% min_C = 0;
% max_C = 1;
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% if showAAL
%     
%     before_Angular_L = 545;
%     after_Angular_L = 6;
%     
%     before_Frontal_Mid_L = 78;
%     after_Frontal_Mid_L = 25;
%     
%     rectangle('Position',[before_Angular_L before_Angular_L after_Angular_L after_Angular_L],'Curvature',0.2);
%     rectangle('Position',[before_Frontal_Mid_L before_Frontal_Mid_L after_Frontal_Mid_L after_Frontal_Mid_L],'Curvature',0.2);
% 
% end
% 
% print(f,'-depsc',strcat('FC-Functional-Clusters','-',label,'.eps'));

new_label = strcat('FC-Functional-Clusters-Petra-Whole-Segment','-',label);
plotFCgeneral(mean_corr,s_corr,new_label);


end

function totalInformationContentPerSubject

load('758-Functional-Clusters-Corr.mat');

nSegments = 4;
nSubjects = 8;
r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

iCorr = 0;
for iSubject=1:nSubjects

    for iSegment=1:nSegments
        
        iCorr = iCorr + 1;
        
        subject(iSubject).corr(iSegment,:,:) = func_clusters(iCorr).rest_corr;
        
    end
    
end

for iSubject=1:nSubjects

    single_rho = subject(iSubject).corr;

    subject(iSubject).fisher = (1/2).*log((1 + single_rho)./(1 - single_rho));
    
end
         
for iSubject=1:nSubjects
    
    my(iSubject).mean = squeeze(mean(subject(iSubject).fisher,1));
    my(iSubject).std = squeeze(std(subject(iSubject).fisher));      
    
end

for iSubject=1:nSubjects
    
    info(iSubject).IM_total = TotalInformationContent( my(iSubject).mean(:), my(iSubject).std(:), z_threshold, cv_threshold );
    
end

end

function totalInformationContentPerSubjectPetra

load('758-Functional-Clusters-Petra-Corr.mat');

nSegments = 3;
nSubjects = 9;
r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

iCorr = 0;
for iSubject=1:nSubjects

    for iSegment=1:nSegments
        
        iCorr = iCorr + 1;
        
        subject(iSubject).corr(iSegment,:,:) = func_clusters(iCorr).rest_corr;
        
    end
    
end

for iSubject=1:nSubjects

    single_rho = subject(iSubject).corr;

    subject(iSubject).fisher = (1/2).*log((1 + single_rho)./(1 - single_rho));
    
end
         
for iSubject=1:nSubjects
    
    my(iSubject).mean = squeeze(mean(subject(iSubject).fisher,1));
    my(iSubject).std = squeeze(std(subject(iSubject).fisher));      
    
end

for iSubject=1:nSubjects
    
    info(iSubject).IM_total = TotalInformationContent( my(iSubject).mean(:), my(iSubject).std(:), z_threshold, cv_threshold );
    
end

end


%%% ANATOMICAL CLUSTER

function getClustersOfVoxelsPer3DCoordinatesPerROI

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nROI = 90;
MNI_img = [91 109 91];
nTotalVoxels = 160990;
nTotalClusters = 758;
nClusterSize = 200;

iiVoxel = 0;

idx = zeros(3,nTotalVoxels);
all_idx_voxels = [];

for iROI=1:nROI
    
    disp(strcat('ROI:',int2str(iROI)));

    nClusters = length(ROI(iROI).clusters);
    
    nVoxels = 0;
    
    all_idx_voxels = [];
    for iCluster=1:nClusters
       
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
        
        nVoxels = nVoxels + length(idx_voxels);
        
        all_idx_voxels = [all_idx_voxels; idx_voxels];
        
    end
    
    if nClusters == 1
        
        ROI_clusters(iROI).clusters(1).idx_voxels = ROI(iROI).clusters(1).idx_voxels;
        
    else
        
        idx = zeros(3,nVoxels);
        
        for iVoxel=1:nVoxels
            
           iiVoxel = iiVoxel + 1;
            
           [idx(1,iiVoxel) idx(2,iiVoxel) idx(3,iiVoxel)] = ind2sub(MNI_img,all_idx_voxels(iVoxel)); 
            
        end
        
        D = pdist(idx','euclidean');

        closest_voxels = zeros(nVoxels,nClusterSize);

        for iVoxel=1:nVoxels
    
            %D((I-1)*(M-I/2)+J-I);
    
            distances = zeros(1,nVoxels);
    
            if mod(iVoxel,20000) == 0; disp(int2str(iVoxel)); end

            distances((iVoxel+1):nVoxels) = D((iVoxel-1)*(nVoxels-iVoxel/2) + ((iVoxel+1):nVoxels)-iVoxel);

            if iVoxel > 1

                distances((iVoxel-1):-1:1) = D((((iVoxel-1):-1:1)-1).*(nVoxels-((iVoxel-1):-1:1)./2)+iVoxel-((iVoxel-1):-1:1));
    
            end

            distances(iVoxel) = 0;

            [s,i] = sort(distances,'ascend');

            closest_voxels(iVoxel,:) = i(1:nClusterSize);
    
        end
        
        idx_clusters = kmeans(closest_voxels,nClusters);
        
        for iCluster=1:nClusters
   
            idx_voxels_on_cluster = find(idx_clusters==iCluster);
    
            ROI_clusters(iROI).clusters(iCluster).idx_voxels = all_idx_voxels(idx_voxels_on_cluster);
    
        end

    end
    
end

save('Anatomical-Spatial-Clusters.mat','ROI_clusters');

end

function getCorrAnatomicalDifferentClusterSameDifferentAAL

load('Anatomical-Spatial-Clusters.mat');

nRuns = 4;
nTotalRuns = 32;
nROI = 90;
nTR = 150;
nSubjects = 8;
nTotalClusters = 758;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

iiRun = 0;
all_settings = getAllSettings;
for iSubject=1:nSubjects
    settings = all_settings(iSubject).settings;
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
    
    kind = 'RestingState';
    for irun=1:nRuns;
        
        iiRun = iiRun + 1;
        
        [RestingState(iiRun).run, RestingState(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    
    end
    
end

iiRun = 0;
all_settings = getAllSettings;
for iSubject=1:nSubjects
    settings = all_settings(iSubject).settings;
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;

    kind = 'Passive';
    for irun=1:nRuns;
        
        iiRun = iiRun + 1;
        
        [Passive(iiRun).run, Passive(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    
    end

end

iiRun = 0;
all_settings = getAllSettings;
for iSubject=1:nSubjects
    settings = all_settings(iSubject).settings;
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;

    kind = 'Track';
    for irun=1:nRuns;
        
        iiRun = iiRun + 1;
        
        [Track(iiRun).run, Track(iiRun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    
    end
    
end

resting_state_clusters = zeros(nTotalClusters,nTR);
passive_clusters = zeros(nTotalClusters,nTR);
track_clusters = zeros(nTotalClusters,nTR);

for iRun=1:nTotalRuns
    
    iiCluster = 0;
    for iROI=1:nROI

        nClusters = length(ROI_clusters(iROI).clusters);

        for iCluster=1:nClusters

            iiCluster = iiCluster + 1;

            idx_voxels = ROI_clusters(iROI).clusters(iCluster).idx_voxels;

            nVoxels = length(idx_voxels);

            resting_state_voxels = zeros(nTR,nVoxels);
            passive_voxels = zeros(nTR,nVoxels);
            track_voxels = zeros(nTR,nVoxels);

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

                resting_state_voxels(:,iVoxel) = RestingState(iRun).run(idxx,idxy,idxz,:);
                passive_voxels(:,iVoxel) = Passive(iRun).run(idxx,idxy,idxz,:);
                track_voxels(:,iVoxel) = Track(iRun).run(idxx,idxy,idxz,:);

            end

            resting_state_clusters(iiCluster,:) = mean(resting_state_voxels',1);
            passive_clusters(iiCluster,:) = mean(passive_voxels',1);
            track_clusters(iiCluster,:) = mean(track_voxels',1);

        end

    end
    
    Spatial_Clusters(iRun).rest_corr = corr(resting_state_clusters');
    Spatial_Clusters(iRun).passive_corr = corr(passive_clusters');
    Spatial_Clusters(iRun).track_corr = corr(track_clusters');
    
end

save('Anatomical-Spatial-Clusters-Corr.mat','Spatial_Clusters');

end

function plotVarianceDistributionAnatomicalDifferentClusterSameDifferentAAL

load('Anatomical-Spatial-Clusters.mat');
load('Anatomical-Spatial-Clusters-Corr.mat');

nROI = 90;
nRuns = 32;

%%% SAME ALL - RESTING STATE

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI_clusters(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        for iiCluster=1:nClusters
            
            iPair = iPair + 1;
        
            for iRun=1:nRuns

                single_rho = Spatial_Clusters(iRun).rest_corr(iiiCluster,iiiCluster-iCluster+iiCluster);
                
                pairs_inside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
    
            end
        
        end

    end
    
end

my_var_inside = var(pairs_inside);

my_cv_inside = std(pairs_inside) ./ abs(mean(pairs_inside,1));

my_mean_inside = mean(pairs_inside,1);
my_std_inside = std(pairs_inside);

vector_one = my_std_inside;
vector_two = my_mean_inside;
vector_three = my_cv_inside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Anatomical-Same-AAL-Rest';

MyContourCorrelation_v2(vector_two,vector_one,title_label,1);

%%% SAME ALL - PASSIVE

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI_clusters(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        for iiCluster=1:nClusters
            
            iPair = iPair + 1;
        
            for iRun=1:nRuns

                single_rho = Spatial_Clusters(iRun).passive_corr(iiiCluster,iiiCluster-iCluster+iiCluster);
                
                pairs_inside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
            end
        
        end

    end
    
end

my_var_inside = var(pairs_inside);

my_cv_inside = std(pairs_inside) ./ abs(mean(pairs_inside,1));

my_mean_inside = mean(pairs_inside,1);
my_std_inside = std(pairs_inside);

vector_one = my_std_inside;
vector_two = my_mean_inside;
vector_three = my_cv_inside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Anatomical-Same-AAL-Passive';

MyContourCorrelation_v2(vector_two,vector_one,title_label,1);

%%% SAME ALL - TRACK

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI_clusters(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        for iiCluster=1:nClusters
            
            iPair = iPair + 1;
        
            for iRun=1:nRuns

                single_rho = Spatial_Clusters(iRun).track_corr(iiiCluster,iiiCluster-iCluster+iiCluster);
                
                pairs_inside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));

            end
        
        end

    end
    
end

my_var_inside = var(pairs_inside);

my_cv_inside = std(pairs_inside) ./ abs(mean(pairs_inside,1));

my_mean_inside = mean(pairs_inside,1);
my_std_inside = std(pairs_inside);

vector_one = my_std_inside;
vector_two = my_mean_inside;
vector_three = my_cv_inside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Anatomical-Same-AAL-Track';

MyContourCorrelation_v2(vector_two,vector_one,title_label,1);

%%% DIFFERENT AAL - RESTING STATE

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI_clusters(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        ivCluster = 0;
        for iiROI=1:nROI
            
            nnClusters = length(ROI_clusters(iiROI).clusters);
            
            for iiCluster=1:nnClusters
                
                ivCluster = ivCluster + 1;
                
                if iROI ~= iiROI
                    
                    iPair = iPair + 1;
                    
                    for iRun=1:nRuns
                        
                        single_rho = Spatial_Clusters(iRun).rest_corr(iiiCluster,ivCluster);
                        
                        pairs_outside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

my_var_outside = var(pairs_outside);

my_cv_outside = std(pairs_outside) ./ abs(mean(pairs_outside,1));

my_mean_outside = mean(pairs_outside,1);
my_std_outside = std(pairs_outside);

vector_one = my_std_outside;
vector_two = my_mean_outside;
vector_three = my_cv_outside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Anatomical-Different-AAL-Rest';

MyContourCorrelation_v2(vector_two,vector_one,title_label,1);

%%% DIFFERENT AAL - PASSIVE

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI_clusters(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        ivCluster = 0;
        for iiROI=1:nROI
            
            nnClusters = length(ROI_clusters(iiROI).clusters);
            
            for iiCluster=1:nnClusters
                
                ivCluster = ivCluster + 1;
                
                if iROI ~= iiROI
                    
                    iPair = iPair + 1;
                    
                    for iRun=1:nRuns
                        
                        single_rho = Spatial_Clusters(iRun).passive_corr(iiiCluster,ivCluster);
                        
                        pairs_outside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

my_var_outside = var(pairs_outside);

my_cv_outside = std(pairs_outside) ./ abs(mean(pairs_outside,1));

my_mean_outside = mean(pairs_outside,1);
my_std_outside = std(pairs_outside);

vector_one = my_std_outside;
vector_two = my_mean_outside;
vector_three = my_cv_outside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Anatomical-Different-AAL-Passive';

MyContourCorrelation_v2(vector_two,vector_one,title_label,1);

%%% DIFFERENT AAL - TRACK

iPair = 0;
iiiCluster = 0;
for iROI=1:nROI
   
    nClusters = length(ROI_clusters(iROI).clusters);

    for iCluster=1:nClusters
        
        iiiCluster = iiiCluster + 1;
        
        ivCluster = 0;
        for iiROI=1:nROI
            
            nnClusters = length(ROI_clusters(iiROI).clusters);
            
            for iiCluster=1:nnClusters
                
                ivCluster = ivCluster + 1;
                
                if iROI ~= iiROI
                    
                    iPair = iPair + 1;
                    
                    for iRun=1:nRuns
                        
                        single_rho = Spatial_Clusters(iRun).track_corr(iiiCluster,ivCluster);
                        
                        pairs_outside(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

my_var_outside = var(pairs_outside);

my_cv_outside = std(pairs_outside) ./ abs(mean(pairs_outside,1));

my_mean_outside = mean(pairs_outside,1);
my_std_outside = std(pairs_outside);

vector_one = my_std_outside;
vector_two = my_mean_outside;
vector_three = my_cv_outside;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'Anatomical-Different-AAL-Track';

MyContourCorrelation_v2(vector_two,vector_one,title_label,1);


end

%%% AAL REGIONS

function getMeanTimeSeriesFrom90AALRunByRun

nROI = 90;
nTotalRuns = 32;
nRuns = 4;
nTR = 150;

all_settings = getAllSettings;

nSubjects = length(all_settings);

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

mean_area_RestingState_voxel = zeros(nROI,nTR);
mean_area_Passive_voxel = zeros(nROI,nTR);
mean_area_Track_voxel = zeros(nROI,nTR);

iiRun = 0;

for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
    
    kind = 'RestingState';
    for irun=1:nRuns;
        [RestingState(irun).run, RestingState(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
    kind = 'Passive';
    for irun=1:nRuns;
        [Passive(irun).run, Passive(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
    kind = 'Track';
    for irun=1:nRuns;
        [Track(irun).run, Track(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
    for irun=1:nRuns
        
        iiRun = iiRun + 1;
        
        disp(strcat('Run:',int2str(iiRun)));
       
        for iROI=1:nROI

            idx_ROI = AAL_ROI(iROI).ID;
            idx_ROI_Label = AAL_ROI(iROI).Nom_L;
            idx_ROI_Label = strrep(idx_ROI_Label,'_','-');

            idx_voxels = find(AAL_img == idx_ROI);

            nVoxels = length(idx_voxels);

            area_RestingState = zeros(nVoxels,nTR);
            area_Passive = zeros(nVoxels,nTR);
            area_Track = zeros(nVoxels,nTR);

            izr = 0;
            izp = 0;
            izt = 0;

            zeros_RestingState = [];
            zeros_Passive = [];
            zeros_Track = [];

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

                area_RestingState(iVoxel,:) = RestingState(irun).run(idxx,idxy,idxz,:);
                area_Passive(iVoxel,:) = Passive(irun).run(idxx,idxy,idxz,:);
                area_Track(iVoxel,:) = Track(irun).run(idxx,idxy,idxz,:);

                if sum(area_RestingState(iVoxel,:)) == 0

                    izr = izr + 1;
                    zeros_RestingState(izr) = iVoxel;

                end
                
                if sum(area_Passive(iVoxel,:)) == 0

                    izp = izp + 1;
                    zeros_Passive(izp) = iVoxel;

                end
                
                if sum(area_Track(iVoxel,:)) == 0

                    izt = izt + 1;
                    zeros_Track(izt) = iVoxel;

                end


            end

            area_RestingState(zeros_RestingState,:) = [];
            area_Passive(zeros_Passive,:) = [];
            area_Track(zeros_Track,:) = [];

            mn_RestingState = mean(area_RestingState,1);
            mn_Passive = mean(area_Passive,1);
            mn_Track = mean(area_Track,1);

            mean_area_RestingState_voxel(iROI,:) = mn_RestingState(:);
            mean_area_Passive_voxel(iROI,:) = mn_Passive(:);
            mean_area_Track_voxel(iROI,:) = mn_Track(:);

        end
        
        mean_run(iiRun).rest = mean_area_RestingState_voxel;
        mean_run(iiRun).passive = mean_area_Passive_voxel;
        mean_run(iiRun).track = mean_area_Track_voxel;
        
        clear mean_area_RestingState_voxel
        clear mean_area_Passive_voxel
        clear mean_area_Track_voxel
        
    end 

end

save('AAL_ROI_mean_run.mat','mean_run');

end

function getMeanTimeSeriesFrom90AALRunByRunPETRA

nROI = 90;
nTotalRuns = 27;
nSegments = 3;
nTR = 150;

all_settings = getAllSettingsPetra;
nSubjects = length(all_settings);

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

run = 1;
all_settings = getAllSettingsPetra;

iRun = 0;
for iSet=1:length(all_settings)
    
    settings = all_settings(iSet).settings;

    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;

    [residual, settings] = real_get_data_FSL_Petra(settings,run,file,get_at_this_preprocessed_step);
    
    for iSeg=1:nSegments
        
        iRun = iRun + 1;
        
        Subject(iRun).residual(:,:,:,:) = residual(:,:,:,(1+(iSeg-1)*nTR):(nTR+(iSeg-1)*nTR));
        
    end

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

            area_RestingState(iVoxel,:) = Subject(iRun).residual(idxx,idxy,idxz,:);

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

save('AAL_ROI_mean_run-Petra.mat','mean_run');

end

function getMeanTimeSeriesFrom90AALRunByRunHCP

nROI = 90;
nTotalRuns = 38*4;
nTotalRunsQuarter = 38*4/4;
nRuns = 4;
nTR = 1200;

all_settings = getAllSettingsHCP;
nSubjects = length(all_settings);

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

iiRun = 0;
for iSet=1:5
    
    settings = all_settings(iSet).settings;
    
    disp(strcat('Subj:',int2str(iSet)));

    for iRun=1:nRuns
        
        iiRun = iiRun + 1;
        
        [all_runs(iiRun).run, settings] = real_get_data_HCP(settings,iRun);
    
    end
    
end

i_global_Run = 0;

for iRun=1:20
    
    i_global_Run = i_global_Run + 1;

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

            area_RestingState(iVoxel,:) = all_runs(iRun).run(idxx,idxy,idxz,:);

            if sum(area_RestingState(iVoxel,:)) == 0

                izr = izr + 1;
                zeros_RestingState(izr) = iVoxel;

            end

        end

        area_RestingState(zeros_RestingState,:) = [];

        mn_RestingState = mean(area_RestingState,1);

        mean_area_RestingState_voxel(iROI,:) = mn_RestingState(:);

    end

    mean_run(i_global_Run).rest = mean_area_RestingState_voxel;

    clear mean_area_RestingState_voxel

end 

iiRun = 0;
for iSet=6:8
    
    settings = all_settings(iSet).settings;
    
    disp(strcat('Subj:',int2str(iSet)));

    for iRun=1:nRuns
        
        iiRun = iiRun + 1;
        
        [all_runs(iiRun).run, settings] = real_get_data_HCP(settings,iRun);
    
    end
    
end

for iRun=1:12
    
    i_global_Run = i_global_Run + 1;

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

            area_RestingState(iVoxel,:) = all_runs(iRun).run(idxx,idxy,idxz,:);

            if sum(area_RestingState(iVoxel,:)) == 0

                izr = izr + 1;
                zeros_RestingState(izr) = iVoxel;

            end

        end

        area_RestingState(zeros_RestingState,:) = [];

        mn_RestingState = mean(area_RestingState,1);

        mean_area_RestingState_voxel(iROI,:) = mn_RestingState(:);

    end

    mean_run(i_global_Run).rest = mean_area_RestingState_voxel;

    clear mean_area_RestingState_voxel

end 

save('AAL_ROI_mean_run-HCP.mat','mean_run');

end

function getCorrFromMeanTimeSeriesFrom90AALRunByRun

load('AAL_ROI_mean_run.mat');

nTotalRuns = 32;

for iRun=1:nTotalRuns
    
   run = mean_run(iRun).rest;
   run = run';
   rho(iRun).rest_corr = corr(run);
   
   run = mean_run(iRun).passive;
   run = run';
   rho(iRun).passive_corr = corr(run);
   
   run = mean_run(iRun).track;
   run = run';
   rho(iRun).track_corr = corr(run);
    
end

save('AAL_ROI_mean_run_corr.mat','rho');

end

function getCorrFromMeanTimeSeriesFrom90AALRunByRunPETRA

load('AAL_ROI_mean_run-Petra.mat');

nTotalRuns = 27;

for iRun=1:nTotalRuns
    
   run = mean_run(iRun).rest;
   run = run';
   rho(iRun).rest_corr = corr(run);
    
end

save('AAL_ROI_mean_run_corr-Petra.mat','rho');

end

function getCorrFromMeanTimeSeriesFrom90AALRunByRunHCP

load('AAL_ROI_mean_run-HCP.mat');

nTotalRuns = 32;

for iRun=1:nTotalRuns
    
   run = mean_run(iRun).rest;
   run = run';
   rho(iRun).rest_corr = corr(run);
    
end

save('AAL_ROI_mean_run_corr-HCP.mat','rho');

end

function getAndPlotVarianceDistributionFromFC90AAL

load('AAL_ROI_mean_run_corr.mat');

nROI = 90;
nTotalRuns = 32;

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

% iPair = 0;
% 
% for iROI=1:nROI
%     
%     for iiROI=1:nROI
%         
%         if iROI ~= iiROI
%             
%             iPair = iPair + 1;
%         
%             for iRun=1:nTotalRuns
%             
%                 single_rho = rho(iRun).passive_corr(iROI,iiROI);
%             
%                 pairs_passive(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
%                 
%             end
%             
%         end
%         
%     end
%     
% end

% iPair = 0;
% 
% for iROI=1:nROI
%     
%     for iiROI=1:nROI
%         
%         if iROI ~= iiROI
%             
%             iPair = iPair + 1;
%         
%             for iRun=1:nTotalRuns
%             
%                 single_rho = rho(iRun).track_corr(iROI,iiROI);
%             
%                 pairs_track(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
%             
%             end
%             
%         end
%         
%     end
%     
% end

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
title_label = 'AAL-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

fraction = FractionSignificantConsistent( vector_two, vector_one, z_threshold, cv_threshold )

% my_var_aal = var(pairs_passive);
% my_cv_aal = std(pairs_passive) ./ abs(mean(pairs_passive,1));
% my_mean_aal = mean(pairs_passive,1);
% my_std_aal = std(pairs_passive);
% 
% vector_one = my_std_aal;
% vector_two = my_mean_aal;
% vector_three = my_cv_aal;
% label_one = 'STD';
% label_two = 'Mean';
% label_three = 'CV';
% title_label = 'AAL-Passive';
% 
% MyContourCorrelation_v5(vector_two,vector_one,title_label,1);
% 
% my_var_aal = var(pairs_track);
% my_cv_aal = std(pairs_track) ./ abs(mean(pairs_track,1));
% my_mean_aal = mean(pairs_track,1);
% my_std_aal = std(pairs_track);
% 
% vector_one = my_std_aal;
% vector_two = my_mean_aal;
% vector_three = my_cv_aal;
% label_one = 'STD';
% label_two = 'Mean';
% label_three = 'CV';
% title_label = 'AAL-Track';
% 
% MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function getAndPlotVarianceDistributionFromFC90AALONESubject(nStartRun,nEndRun)

load('AAL_ROI_mean_run_corr.mat');

nROI = 90;
nTotalRuns = 4;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            iiRun = 0;
            for iRun=nStartRun:nEndRun
                
                iiRun = iiRun + 1;
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest(iiRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
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
title_label = strcat('AAL-Rest','-Runs-',int2str(nStartRun),'-',int2str(nEndRun));

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

fraction = FractionSignificantConsistent( vector_two, vector_one, z_threshold, cv_threshold )

end

function getAndPlotVarianceDistributionFromFC90AALPETRA

load('AAL_ROI_mean_run_corr-Petra.mat');

nROI = 90;
nTotalRuns = 27;

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
title_label = 'AAL-Petra-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function getAndPlotVarianceDistributionFromFC90AALHCP

load('AAL_ROI_mean_run_corr-HCP.mat');

nROI = 90;
nTotalRuns = 32;

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
title_label = 'AAL-HCP-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function getAndPlotVarianceDistributionFromFC90AALInsideVoxelLevel

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_ROI = 1:90;
nRuns = 4;
nROI = length(idx_ROI);
nTotalRuns = 32;

all_settings = getAllSettings;

nSubjects = length(all_settings);

for iiRun=1:nTotalRuns
    
    aal_rhos(iiRun).rho = [];
    
end

for iROI=1:nROI

    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    area_label = strrep(label_ROI{iROI},'_','-');       
    disp(area_label);
        
    iiRun = 0;
    for iSubject=1:length(all_settings)
        
        settings = all_settings(iSubject).settings;
  
        RestingStateRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr','-','RestingState','-',area_label,'.mat'));
        
        for irun=1:nRuns
            
            iiRun = iiRun + 1;
       
            rest_rho = RestingStateRHOs.FC_Voxels.run(irun).rho_RestingState;
       
            aal_rhos(iiRun).rho = [aal_rhos(iiRun).rho(:)',rest_rho(:)'];
            
        end
     
    end

end

nRhos = length(aal_rhos(1).rho);

AAL_rhos = zeros(nTotalRuns,nRhos);

for iRun=1:nTotalRuns
    
    AAL_rhos(iRun,:) = aal_rhos(iRun).rho(:);

end

clear aal_rhos
            
AAL_rhos_z = (1/2).*log((1 + AAL_rhos)./(1 - AAL_rhos));

my_mean_aal = mean(AAL_rhos_z,1);
my_std_aal = std(AAL_rhos_z,0,1);

vector_one = my_std_aal;
vector_two = my_mean_aal;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'AAL-Inside-Rest';

r_threshold = 0.186728;  % corresponds to p-value = 0.01 for N=150
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1;

IM_total = TotalInformationContent( vector_two, vector_one, z_threshold, cv_threshold );

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotMeanFCFrom90AAL

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

label = 'FC-AAL';
plotFCgeneral(mean_corr,s_corr,label);

% f = figure;
% 
% clmap = colormap('jet');
% clmap(65,:) = clmap(64,:);
% clmap(33,:) = [1 1 1];
% 
% min_C = -1;
% max_C = 1;
% 
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% 
% print(f,'-depsc','FC-AAL.eps');
% export_fig('FC-AAL','-pdf');

end

function plotMeanFCFrom90AALONESubject(nStartRun,nEndRun)

load('AAL_ROI_mean_run_corr.mat');

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
        
            iiRun = 0;
            for iRun=nStartRun:nEndRun
                
                iiRun = iiRun + 1;
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest_z(iiRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                pairs_rest_rho(iiRun,iPair) = single_rho;
            
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

label = strcat('FC-AAL','-Runs-',int2str(nStartRun),'-',int2str(nEndRun));
plotFCgeneral(mean_corr,s_corr,label);

% f = figure;
% 
% clmap = colormap('jet');
% clmap(65,:) = clmap(64,:);
% clmap(33,:) = [1 1 1];
% 
% min_C = -1;
% max_C = 1;
% 
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% 
% print(f,'-depsc','FC-AAL.eps');
% export_fig('FC-AAL','-pdf');

end

function plotMeanFCFrom90AALPETRA

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

label = 'FC-AAL-Petra';
plotFCgeneral(mean_corr,s_corr,label);

% f = figure;
% 
% clmap = colormap('jet');
% clmap(65,:) = clmap(64,:);
% clmap(33,:) = [1 1 1];
% 
% min_C = -1;
% max_C = 1;
% 
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% 
% print(f,'-depsc','FC-AAL.eps');
% export_fig('FC-AAL','-pdf');

end

function plotMeanFCFrom90AALHCP

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

label = 'FC-AAL-HCP';
plotFCgeneral(mean_corr,s_corr,label);

% f = figure;
% 
% clmap = colormap('jet');
% clmap(65,:) = clmap(64,:);
% clmap(33,:) = [1 1 1];
% 
% min_C = -1;
% max_C = 1;
% 
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% 
% print(f,'-depsc','FC-AAL.eps');
% export_fig('FC-AAL','-pdf');

end

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










  