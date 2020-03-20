function real_methods_paper_v1_top


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

function getClustersOfVoxelsPer3DCoordinatesV2

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nClusters = 758;
nROI = 90;
nRuns = 32;
MNI_img = [91 109 91];
nTotalVoxels = 160990;
nTotalClusters = 758;
nClusterSize = 200;

for iROI=1:nROI
    
    disp(ROI(iROI).label);
    
    nTotalVoxels = ROI(iROI).nVoxels;
    
    iiVoxel = 0;
    all_idx_voxels = [];
    idx = [];
    iiCluster = 0;

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
    
    disp('...getting Euclidean');
    D = pdist(idx','euclidean');
    
    min_size = min(nClusterSize*5,nTotalVoxels);
    
    closest_voxels = zeros(nTotalVoxels,min_size);

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
        
        closest_voxels(iVoxel,:) = i(1:min_size);

    end
    
    disp('...doing kmeans');
    k = nClusters;
    idx_clusters = kmeans(closest_voxels,k);
    
    for iCluster=1:k
        
        iiCluster = iiCluster + 1;
   
        idx_voxels_on_cluster = find(idx_clusters==iCluster);
    
        clusters(iiCluster).idx_voxels = all_idx_voxels(idx_voxels_on_cluster);
    
    end
    
    save(strcat('Clusters-Voxels-Per-3D-Coordinate-v2-',ROI(iROI).label,'.mat'),'clusters','idx');
    
    clear clusters idx
    
end

end

function save3DVolumeWithSpatialClustersInOneROIV2(idx_ROI)

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

load(strcat('Clusters-Voxels-Per-3D-Coordinate-v2-',ROI(idx_ROI).label,'.mat'));

MNI_size = [91 109 91];

my_volume = zeros(MNI_size);

nClusters = length(clusters);

for iCluster=1:nClusters
    
    idx_voxels = clusters(iCluster).idx_voxels;
    
    nVoxels = length(idx_voxels);
    
    for iVoxel=1:nVoxels
        
       [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
       
       my_volume(idxx,idxy,idxz) = iCluster;
        
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

fname = strcat('Clusters-Voxels-Per-3D-Coordinate-v2-',ROI(idx_ROI).label,'.nii');
input_data = my_volume; 
real_save_image;

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


end

