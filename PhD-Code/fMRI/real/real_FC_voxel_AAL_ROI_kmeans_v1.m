function real_FC_voxel_AAL_ROI_kmeans_v1

%settings_subj1_2210;
%settings_subj2_2610;
%settings_subj3_0311;
%settings_subj4_0211;
%settings_subj5_0211;
%settings_subj6_2411;

% settings_subj2_2610;
% doTheMath(settings);
% clear settings

% settings_subj3_0311;
% doTheMath(settings);
% clear settings

% settings_subj4_0211;
% doTheMath(settings);
% clear settings

% settings_subj5_0211;
% doTheMath(settings);
% clear settings

% settings_subj6_2411;
% doTheMath(settings);
% clear settings

% nROI = 90;
% idx_ROI = 8; %% Frontal Mid R - 5104 voxels
% idx_ROI = 1:nROI;
% 
% for iROI=[55 56 87 88]
%     
%     if iROI ~= 8
%     
%         clusterConcatenatedAllSubjects(iROI);
%         
%     end
%     
% end

% for iROI=88:-1:64
%  
%     clusterConcatenatedAllSubjectsPETRA(iROI);
%  
% end

%plotDistributionVoxelsConcatenatedAllSubjects(idx_ROI);
%saveSortedFirstClusterAllRuns(idx_ROI);

% for iROI=1:90
%     save3DimageWithClusterInformation(iROI);
% end

% giveMeAmountOfVoxelsAndClusters;
% getClusterParcellationInfoAndMeanTimeSeries;
% saveVolumeWithClusterInfo;

% getCorrelationsOfMeansInsideClusters;
% 
% getDensitiesOfCorrelationsOfMeansInsideClusters;
% 
% doTTestOnDensitiesOfCorrelationsOFMeansInsideClusters;

% doFCOnAllClusters;
% plotFCOnAllClustersMeanRuns;
% plotFCOnAllClustersPerRun;

% doContrastsFCOnAllClusters;
% plotContrastFCOnAllClusters;

% computeGrangerAmongClusters;

% computeGrangerAtGroupLevel;

% computeGrangerAtGroupLevelCorbettasMethod;

% plotGrangerAtGroupLevelMean;

% plotGrangerAtGroupLevelContrast;

% plotGrangerAtGroupLevelContrastDTI;

% computeFunctionalParcellationForClusters;

% plotFCAndGrangerAtGroupLevelContrastJB;

% plotGrangerAtGroupLevelContrastJBSubcortical;
% plotFCClusterContrastJBSubcortical;

% plotGrangerAtGroupLevelContrastJBAllLobes;
% plotFCClusterContrastJBAllLobes;

% getMeaningfulClustersWithChangeBasedOnFC;

% getMeaningfulClustersWithChangeBasedOnGranger;

% getMeaningfulClustersWithChangeBasedOnFCGCAttentionStimulusOnly;

% plotMeaningfulClustersWithChangeBasedOnFCGCAttentionOnly;

% plotFCAndGrangerAtGroupLevelContrastOnlySeeds;

% plotFCAndGrangerAtGroupLevelContrastOnlySeedsSortedPerNetwork;

% plotFCAndGrangerAtGroupLevelMeanOnlySeedsSortedPerNetwork;

% whoIsTheStrongestConnection;

% plotSignOfChangeAndCommonOrigin;

% giveMeAmountOfVoxelsAndClustersPETRA;

% saveVolumeWithClusterInfoPETRA;

% getClusterParcellationInfoAndMeanTimeSeries;

saveVolumeWithClusterInfo;

end

function doTheMath(settings)
        
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

for iROI=1:nROI
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    TrackRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'.mat'));
    PassiveRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'.mat'));
    RestingStateRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'.mat'));
        
     for irun=1:nRuns
       
        disp(strcat('Track:','irun=',int2str(irun)));
        
        track_rho = TrackRHOs.FC_Voxels.run(irun).rho_Track;
        track_pval = TrackRHOs.FC_Voxels.run(irun).pval_Track;

        
        [FC_KMeans.run(irun).IdxClusters, FC_KMeans.run(irun).Tidx, FC_KMeans.run(irun).Ncluster] = ClusterWithKmeans( track_rho, track_pval );
        
     end
    
     save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-KMeans','.mat'),'area_label','FC_KMeans');
     clear FC_KMeans
     
     for irun=1:nRuns
       
        disp(strcat('Passive:','irun=',int2str(irun)));
        
        passive_rho = PassiveRHOs.FC_Voxels.run(irun).rho_Passive;
        passive_pval = PassiveRHOs.FC_Voxels.run(irun).pval_Passive;
        
        [FC_KMeans.run(irun).IdxClusters, FC_KMeans.run(irun).Tidx, FC_KMeans.run(irun).Ncluster] = ClusterWithKmeans( passive_rho, passive_pval );
        
     end
     
     save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-KMeans','.mat'),'area_label','FC_KMeans');
     clear FC_KMeans
      
     for irun=1:nRuns
       
        disp(strcat('RestingState:','irun=',int2str(irun)));

        rest_rho = RestingStateRHOs.FC_Voxels.run(irun).rho_RestingState;
        rest_pval = RestingStateRHOs.FC_Voxels.run(irun).pval_RestingState;

        [FC_KMeans.run(irun).IdxClusters, FC_KMeans.run(irun).Tidx, FC_KMeans.run(irun).Ncluster] = ClusterWithKmeans( rest_rho, rest_pval );
        
     end
     
     save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-KMeans','.mat'),'area_label','FC_KMeans');
     clear FC_KMeans
     
end
     
end

function clusterConcatenatedAllSubjects(idx_ROI)

nVoxels_per_Clusters = 200;
        
disp(datestr(now));

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_voxels = find(AAL_img==AAL_ROI(idx_ROI).ID);

nVoxels = length(idx_voxels);

nRuns = 4;

all_settings = getAllSettings;

nSubjects = length(all_settings);

label_ROI = AAL_ROI(idx_ROI).Nom_L;
    
area_label = strrep(label_ROI,'_','-');       

disp(area_label);

all_rhos = [];
all_pvals = [];
    
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    disp(settings.codes.subject);
    
    RestingStateRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr','-','RestingState','-',area_label,'.mat'));
    
    for irun=1:nRuns
       
        disp(strcat('RestingState:','irun=',int2str(irun)));

        rest_rho = RestingStateRHOs.FC_Voxels.run(irun).rho_RestingState;
        rest_pval = RestingStateRHOs.FC_Voxels.run(irun).pval_RestingState;
        
        all_rhos = [all_rhos,rest_rho];
        all_pvals = [all_pvals,rest_pval];

    end
     
    clear RestingStateRHOs

end

[x,y] = size(all_rhos);
[l,c] = size(all_pvals);

% ii = 0;
% for i=1:c
%     gett = squeeze(all_rhos(:,i));
%     if length(find(isnan(gett))) == l
%         ii = ii + 1;
%         idx_c_nan(ii) = i;
%     end
% end
% 
% all_rhos(:,idx_c_nan(:)) = [];
% all_pvals(:,idx_c_nan(:)) = [];
%     
% ii = 0;
% for i=1:l
%     gett = squeeze(all_rhos(i,:));
%     if length(find(isnan(gett))) >= 1
%         ii = ii + 1;
%         idx_l_nan(ii) = i;
%     end
% end
% 
% all_rhos(idx_l_nan(:),:) = [];
% all_pvals(idx_l_nan(:),:) = [];

nClusters = floor(nVoxels/nVoxels_per_Clusters);

disp(strcat('nVoxels:',int2str(nVoxels)));
disp(strcat('rhos:',int2str(x),'-',int2str(y)));
disp(strcat('pvals:',int2str(l),'-',int2str(c)));

[FC_KMeans.all_runs.IdxClusters, FC_KMeans.all_runs.Ncluster, distance] = ClusterWithKmeans_no_symmetrical( all_rhos, all_pvals, nClusters );
        
save(strcat(settings.codes.experiment,'-','All-Subjects','-','FC-Voxels-AAL-ROI-corr','-','RestingState','-',area_label,'-KMeans','.mat'),'area_label','nVoxels','FC_KMeans','distance');
        
end

function clusterConcatenatedAllSubjectsPETRA(idx_ROI)

nVoxels_per_Clusters = 200;
        
disp(datestr(now));

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_voxels = find(AAL_img==AAL_ROI(idx_ROI).ID);

nVoxels = length(idx_voxels);

nRuns = 3;

all_settings = getAllSettingsPetra;

nSubjects = length(all_settings);

label_ROI = AAL_ROI(idx_ROI).Nom_L;
    
area_label = strrep(label_ROI,'_','-');       

disp(area_label);

all_rhos = [];
all_pvals = [];
    
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    disp(settings.codes.subject);
    
    RestingStateRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-Petra','-','RestingState','-',area_label,'.mat'));
    
    for irun=1:nRuns
       
        disp(strcat('RestingState:','irun=',int2str(irun)));

        rest_rho = RestingStateRHOs.FC_Voxels.run(irun).rho_RestingState;
        rest_pval = RestingStateRHOs.FC_Voxels.run(irun).pval_RestingState;
        
        all_rhos = [all_rhos,rest_rho];
        all_pvals = [all_pvals,rest_pval];

    end
     
    clear RestingStateRHOs

end

[x,y] = size(all_rhos);
[l,c] = size(all_pvals);

% ii = 0;
% for i=1:c
%     gett = squeeze(all_rhos(:,i));
%     if length(find(isnan(gett))) == l
%         ii = ii + 1;
%         idx_c_nan(ii) = i;
%     end
% end
% 
% all_rhos(:,idx_c_nan(:)) = [];
% all_pvals(:,idx_c_nan(:)) = [];
%     
% ii = 0;
% for i=1:l
%     gett = squeeze(all_rhos(i,:));
%     if length(find(isnan(gett))) >= 1
%         ii = ii + 1;
%         idx_l_nan(ii) = i;
%     end
% end
% 
% all_rhos(idx_l_nan(:),:) = [];
% all_pvals(idx_l_nan(:),:) = [];

nClusters = floor(nVoxels/nVoxels_per_Clusters);

disp(strcat('nVoxels:',int2str(nVoxels)));
disp(strcat('rhos:',int2str(x),'-',int2str(y)));
disp(strcat('pvals:',int2str(l),'-',int2str(c)));

[FC_KMeans.all_runs.IdxClusters, FC_KMeans.all_runs.Ncluster, distance] = ClusterWithKmeans_no_symmetrical( all_rhos, all_pvals, nClusters );
        
save(strcat(settings.codes.experiment,'-','All-Subjects','-','FC-Voxels-AAL-ROI-corr','-','RestingState','-',area_label,'-KMeans','.mat'),'area_label','nVoxels','FC_KMeans','distance');
        
end

function plotDistributionVoxelsConcatenatedAllSubjects(idx_ROI)

load(strcat('LHR','-','All-Subjects','-','FC-Voxels-AAL-ROI-corr','-','RestingState','-',area_label,'-KMeans','.mat'));
 
nRuns = 32;
n_c_nans = length(idx_c_nan);

idx_voxels_concatenated = 1:(nVoxels*nRuns);

idx_runs_concatenated = [];

for iRun=1:nRuns
    
    new = repmat(iRun,[1 nVoxels]);
    
    idx_runs_concatenated = [idx_runs_concatenated, new];
    
end

idx_voxels_concatenated_nonans = 1:((nVoxels*nRuns)-n_c_nans);


end

function saveSortedFirstClusterAllRuns(idx_ROI)

settings.codes.experiment = 'LHR';

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_ROI_voxels = find(AAL_img==AAL_ROI(idx_ROI).ID);

nVoxels = length(idx_ROI_voxels);

nRuns = 4;

all_settings = getAllSettings;

nSubjects = length(all_settings);

label_ROI = AAL_ROI(idx_ROI).Nom_L;
    
area_label = strrep(label_ROI,'_','-');       

disp(area_label);

%%% LOAD CLUSTER INFO

load(strcat(settings.codes.experiment,'-','All-Subjects','-','FC-Voxels-AAL-ROI-corr','-','RestingState','-',area_label,'-KMeans','.mat'));

idx_cluster_voxels = FC_KMeans.all_runs.IdxClusters;
nCluster = FC_KMeans.all_runs.Ncluster;

disp(strcat('TargetOrder',datestr(now)));
source_order = 1:nVoxels;
[target_order, cluster_idx_sorted, cluster_size_sorted] = TargetOrder( idx_cluster_voxels, nCluster );

rest_rho_1st_cluster = [];

iirun = 0;
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    disp(settings.codes.subject);
    
    RestingStateRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr','-','RestingState','-',area_label,'.mat'));
    
    for irun=1:nRuns
        
        iirun = iirun + 1;
       
        disp(strcat('RestingState:','irun=',int2str(irun)));

        rest_rho = RestingStateRHOs.FC_Voxels.run(irun).rho_RestingState;
        %rest_pval = RestingStateRHOs.FC_Voxels.run(irun).pval_RestingState;
        
        plot_label = strcat('RestingState-',int2str(iirun),'-',area_label);

        rest_rho_1st_cluster = [rest_rho_1st_cluster,rest_rho(:,target_order(1:cluster_size_sorted(1)))];
        
    end
     
    clear RestingStateRHOs
    clear rest_rho
    clear rest_pval

end

settings.codes.subject = 'All-Subjects';

save(strcat(settings.codes.experiment,'-','All-Subjects','-','FC-Voxels-AAL-ROI-corr','-','RestingState','-',area_label,'-KMeans-1st','.mat'),'rest_rho_1st_cluster');

%plotFCClusteredOnly1stCluster(settings,rest_rho_1st_cluster,plot_label);

end

function plotFCClusteredOnly1stCluster(settings,rho_sorted,plot_label)
    
    rlimit = 0.5;

    rmin = -rlimit + rlimit/32;
    rmax =  rlimit;

    Tick = [1 Ncomponent];
    TickLabel = { '1' ; num2str(Ncomponent) };

    f = figure;

    clrmp = colormap('jet');
    clrmp(32,:) = [1 1 1];

    hold on;
    caxis([rmin rmax]);
    h = pcolor( 0.5:Ncomponent-0.5, 0.5:Ncomponent-0.5, rho_sorted );
    set( h, 'EdgeColor', 'none');
    colormap( clrmp);
    hold off;
    axis 'square'; 
    axis([0 Ncomponent 0 Ncomponent]);
    xlabel(plot_label);
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

    h = colorbar;
    set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1]);

    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'.jpg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'.eps'));
    print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC_voxel_per_ROI_AAL-',plot_label,'.pdf'));
    
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

function save3DimageWithClusterInformation(idx_ROI)

settings.codes.experiment = 'LHR';

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_ROI_voxels = find(AAL_img==AAL_ROI(idx_ROI).ID);

nVoxels = length(idx_ROI_voxels);

nRuns = 4;

all_settings = getAllSettings;

nSubjects = length(all_settings);

label_ROI = AAL_ROI(idx_ROI).Nom_L;
    
area_label = strrep(label_ROI,'_','-');       

disp(area_label);

%%% LOAD CLUSTER INFO

load(strcat(settings.codes.experiment,'-','All-Subjects','-','FC-Voxels-AAL-ROI-corr','-','RestingState','-',area_label,'-KMeans','.mat'));

idx_cluster_voxels = FC_KMeans.all_runs.IdxClusters;
nCluster = FC_KMeans.all_runs.Ncluster;

area_img = zeros(size(AAL_img));

for iVoxel=1:nVoxels
    
    [x,y,z] = ind2sub(size(area_img),idx_ROI_voxels(iVoxel));
    
    area_img(x,y,z) = idx_cluster_voxels(iVoxel);
    
end

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = load_aal.dat.dim;

descrip = 'KMeans';

fname = strcat('LHR','-','All-Subjects','-','FC-Voxels-AAL-ROI-corr','-','RestingState','-',area_label,'-KMeans','.nii');
input_data = area_img; 
real_save_image;

end

function giveMeAmountOfVoxelsAndClusters

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nROI = 90;

for iROI=1:nROI
    
    idx_voxels = find(AAL_img==AAL_ROI(iROI).ID);
    
    label_ROI = AAL_ROI(iROI).Nom_L;
    
    area_label = strrep(label_ROI,'_','-');       

    disp(area_label);
    
    ROI(iROI).label = area_label;
    ROI(iROI).nVoxels = length(idx_voxels);
    
    load(strcat('LHR','-','All-Subjects','-','FC-Voxels-AAL-ROI-corr','-','RestingState','-',area_label,'-KMeans','.mat'));
   
    idx_clusters = FC_KMeans.all_runs.IdxClusters;
    
    ROI(iROI).nClusters = FC_KMeans.all_runs.Ncluster;

    for iClu=1:ROI(iROI).nClusters
       
        idx_this_cluster = find(idx_clusters == iClu);
        
        idx_this_cluster_voxels = idx_voxels(idx_this_cluster);
        
        ROI(iROI).clusters(iClu).nVoxels = length(idx_this_cluster);
        ROI(iROI).clusters(iClu).idx_voxels = idx_this_cluster_voxels;
        
    end
    
end

save('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat','ROI');
    
end

function giveMeAmountOfVoxelsAndClustersPETRA

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nROI = 90;

for iROI=1:nROI
    
    idx_voxels = find(AAL_img==AAL_ROI(iROI).ID);
    
    label_ROI = AAL_ROI(iROI).Nom_L;
    
    area_label = strrep(label_ROI,'_','-');       

    disp(area_label);
    
    ROI(iROI).label = area_label;
    ROI(iROI).nVoxels = length(idx_voxels);
    
    load(strcat('PETRA','-','All-Subjects','-','FC-Voxels-AAL-ROI-corr','-','RestingState','-',area_label,'-KMeans','.mat'));
   
    idx_clusters = FC_KMeans.all_runs.IdxClusters;
    
    ROI(iROI).nClusters = FC_KMeans.all_runs.Ncluster;

    for iClu=1:ROI(iROI).nClusters
       
        idx_this_cluster = find(idx_clusters == iClu);
        
        idx_this_cluster_voxels = idx_voxels(idx_this_cluster);
        
        ROI(iROI).clusters(iClu).nVoxels = length(idx_this_cluster);
        ROI(iROI).clusters(iClu).idx_voxels = idx_this_cluster_voxels;
        
    end
    
end

save('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat','ROI');
    
end

function getClusterParcellationInfoAndMeanTimeSeries

load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

all_settings = getAllSettings;
nSet = length(all_settings);

nRuns = 4;

for iSet=1:nSet
    
    %%% LOAD DATA
    
    settings = all_settings(iSet).settings;

    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    %% file = settings.FSL.files.functional.custom.filtered; BUG FOUND ON 18/01/2017
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;

    real_load_all_data_FSL;
    
    for iRun=1:nRuns
        
        All_Track(iSet).runs(iRun).run = Track(iRun).run;
        All_Passive(iSet).runs(iRun).run = Passive(iRun).run;
        All_Resting(iSet).runs(iRun).run = RestingState(iRun).run;
        
    end
    
    clear Track
    clear Passive
    clear RestingState
    
end

nROI = 90;
nTR = 150;
MNI_size = [91 109 91];

for iROI=1:nROI
    
    disp(ROI(iROI).label);
    
    nClusters = ROI(iROI).nClusters;
    
    disp(strcat('nClusters:',int2str(nClusters)));
    
    for iClu=1:nClusters

        idx_voxels = ROI(iROI).clusters(iClu).idx_voxels;
        nVoxels = length(idx_voxels);
        
        iirun = 0;
        
        for iSet=1:nSet
            
            for iRun=1:nRuns
                
                iirun = iirun + 1;
                
                cluster_voxels_track = zeros(nVoxels,nTR);
                cluster_voxels_passive = zeros(nVoxels,nTR);
                cluster_voxels_rest = zeros(nVoxels,nTR);
        
                for iVoxel=1:nVoxels
            
                    [x,y,z] = ind2sub(MNI_size,idx_voxels(iVoxel));
            
                    cluster_voxels_track(iVoxel,:) = squeeze(All_Track(iSet).runs(iRun).run(x,y,z,:));
                    cluster_voxels_passive(iVoxel,:) = squeeze(All_Passive(iSet).runs(iRun).run(x,y,z,:));
                    cluster_voxels_rest(iVoxel,:) = squeeze(All_Resting(iSet).runs(iRun).run(x,y,z,:));
        
                end
                
                ROI(iROI).clusters(iClu).track_mean(iirun,:) = mean(cluster_voxels_track,1);
                ROI(iROI).clusters(iClu).passive_mean(iirun,:) = mean(cluster_voxels_passive,1);
                ROI(iROI).clusters(iClu).rest_mean(iirun,:) = mean(cluster_voxels_rest,1);
        
                clear cluster_voxels_track
                clear cluster_voxels_passive
                clear cluster_voxels_rest
                
            end
            
        end
        
    end
        
end

save('FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat','ROI');

end

function saveVolumeWithClusterInfo

load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nROI = 90;
MNI_size = [91 109 91];

FunctionalClusterParcellationVolume = zeros(MNI_size);
FunctionalClusterParcellationVolume_For_Plot = zeros(MNI_size);

iiclu = 0;

for iROI=1:nROI
    
    disp(ROI(iROI).label);
    
    nClusters = ROI(iROI).nClusters;
    
    disp(strcat('nClusters:',int2str(nClusters)));
    
    for iClu=1:nClusters

        iiclu = iiclu + 1;
        
        idx_voxels = ROI(iROI).clusters(iClu).idx_voxels;
        nVoxels = length(idx_voxels);
        
         for iVoxel=1:nVoxels
            
            [x,y,z] = ind2sub(MNI_size,idx_voxels(iVoxel));
            
            FunctionalClusterParcellationVolume(x,y,z) = iiclu;
            
            FunctionalClusterParcellationVolume_For_Plot(x,y,z) = iClu;
                    
         end
        
    end
    
end

load_aal = nifti('ROI_MNI_V4.nii');

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = load_aal.dat.dim;

descrip = 'FunctionalClustersParcellation';

fname = strcat('LHR','-','All-Subjects','-','FC-Voxel-AAL-ROI-KMeans-Parcellation','.nii');
input_data = FunctionalClusterParcellationVolume; 
% real_save_image;

load_aal = nifti('ROI_MNI_V4.nii');

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = load_aal.dat.dim;

descrip = 'FunctionalClustersParcellation';

fname = strcat('LHR','-','All-Subjects','-','FC-Voxel-AAL-ROI-KMeans-Parcellation-For-Plot','.nii');
input_data = FunctionalClusterParcellationVolume_For_Plot; 
real_save_image;

end

function saveVolumeWithClusterInfoPETRA

load('PETRA-FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nROI = 90;
MNI_size = [91 109 91];

FunctionalClusterParcellationVolume = zeros(MNI_size);
FunctionalClusterParcellationVolume_For_Plot = zeros(MNI_size);

iiclu = 0;

for iROI=1:nROI
    
    disp(ROI(iROI).label);
    
    nClusters = ROI(iROI).nClusters;
    
    disp(strcat('nClusters:',int2str(nClusters)));
    
    for iClu=1:nClusters

        iiclu = iiclu + 1;
        
        idx_voxels = ROI(iROI).clusters(iClu).idx_voxels;
        nVoxels = length(idx_voxels);
        
         for iVoxel=1:nVoxels
            
            [x,y,z] = ind2sub(MNI_size,idx_voxels(iVoxel));
            
            FunctionalClusterParcellationVolume(x,y,z) = iiclu;
            
            FunctionalClusterParcellationVolume_For_Plot(x,y,z) = iClu;
                    
         end
        
    end
    
end

load_aal = nifti('ROI_MNI_V4.nii');

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = load_aal.dat.dim;

descrip = 'FunctionalClustersParcellation';

fname = strcat('PETRA','-','All-Subjects','-','FC-Voxel-AAL-ROI-KMeans-Parcellation','.nii');
input_data = FunctionalClusterParcellationVolume; 
real_save_image;

load_aal = nifti('ROI_MNI_V4.nii');

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;

dtype = 'FLOAT32';
offset = 0;

dim = load_aal.dat.dim;

descrip = 'FunctionalClustersParcellation';

fname = strcat('PETRA','-','All-Subjects','-','FC-Voxel-AAL-ROI-KMeans-Parcellation-For-Plot','.nii');
input_data = FunctionalClusterParcellationVolume_For_Plot; 
real_save_image;

end

function getCorrelationsOfMeansInsideClusters

load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nRuns = 32;
nROI = 90;

for iRun=1:nRuns
   
   for iROI=1:nROI
   
        nClusters = length(ROI(iROI).clusters);
    
        for iCluster=1:nClusters
    
            track_mean(iCluster,:) = squeeze(ROI(iROI).clusters(iCluster).track_mean(iRun,:));
            passive_mean(iCluster,:) = squeeze(ROI(iROI).clusters(iCluster).passive_mean(iRun,:));
            rest_mean(iCluster,:) = squeeze(ROI(iROI).clusters(iCluster).rest_mean(iRun,:));
            
        end
        
        [FC_means_clusters(iROI).run(iRun).track_corr FC_means_clusters(iROI).run(iRun).track_pval] = corr(track_mean');
        [FC_means_clusters(iROI).run(iRun).passive_corr FC_means_clusters(iROI).run(iRun).passive_pval] = corr(passive_mean');
        [FC_means_clusters(iROI).run(iRun).rest_corr FC_means_clusters(iROI).run(iRun).rest_pval]= corr(rest_mean');
        
        clear track_mean
        clear passive_mean
        clear rest_mean
        
    end
    
end

save('FC-Voxels-AAL-ROI-corr-KMeans-FC-Inside-Clusters.mat','FC_means_clusters');

end

function getDensitiesOfCorrelationsOfMeansInsideClusters

load('FC-Voxels-AAL-ROI-corr-KMeans-FC-Inside-Clusters.mat');

nRuns = 32;
nROI = 90;

pcriterion = 0.01;

for iRun=1:nRuns
   
   for iROI=1:nROI
       
       track_corr = FC_means_clusters(iROI).run(iRun).track_corr;
       track_pval = FC_means_clusters(iROI).run(iRun).track_pval;
       
       track_corr(track_pval > pcriterion) = 0;
       
       [l,c] = size(track_corr);
       nClusters = l;
       
       for iCluster=1:nClusters
           
           lines = squeeze(track_corr(:,iCluster));
           
           nPositive = length(find(lines > 0))/nClusters;
           nNegative = length(find(lines < 0))/nClusters;
       
           FC_densities_clusters(iROI).track_positive(iRun,iCluster) = nPositive;
           FC_densities_clusters(iROI).track_negative(iRun,iCluster) = nNegative;
           
       end
       
       passive_corr = FC_means_clusters(iROI).run(iRun).passive_corr;
       passive_pval = FC_means_clusters(iROI).run(iRun).passive_pval;
       
       passive_corr(passive_pval > pcriterion) = 0;
       
       [l,c] = size(passive_corr);
       nClusters = l;
       
       for iCluster=1:nClusters
           
           lines = squeeze(passive_corr(:,iCluster));
           
           nPositive = length(find(lines > 0))/nClusters;
           nNegative = length(find(lines < 0))/nClusters;
       
           FC_densities_clusters(iROI).passive_positive(iRun,iCluster) = nPositive;
           FC_densities_clusters(iROI).passive_negative(iRun,iCluster) = nNegative;
           
       end
       
       rest_corr = FC_means_clusters(iROI).run(iRun).rest_corr;
       rest_pval = FC_means_clusters(iROI).run(iRun).rest_pval;
       
       rest_corr(rest_pval > pcriterion) = 0;
       
       [l,c] = size(rest_corr);
       nClusters = l;
       
       for iCluster=1:nClusters
           
           lines = squeeze(rest_corr(:,iCluster));
           
           nPositive = length(find(lines > 0))/nClusters;
           nNegative = length(find(lines < 0))/nClusters;
       
           FC_densities_clusters(iROI).rest_positive(iRun,iCluster) = nPositive;
           FC_densities_clusters(iROI).rest_negative(iRun,iCluster) = nNegative;
           
       end
       
   end
   
end

save('FC-Voxels-AAL-ROI-corr-KMeans-FC-Densities-Clusters.mat','FC_densities_clusters');

end

function doTTestOnDensitiesOfCorrelationsOFMeansInsideClusters

Bonferroni = false;

if Bonferroni; nTotalClustersBonferroni = 758; else nTotalClustersBonferroni = 1; end
if Bonferroni; label = 'Bonf'; else label = 'noBonf'; end

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

img_pos_stimulus = zeros(size(AAL_img));
img_neg_stimulus = zeros(size(AAL_img));

img_pos_attention = zeros(size(AAL_img));
img_neg_attention = zeros(size(AAL_img));

load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-FC-Densities-Clusters.mat');

nROI = 90;
nTotalClusters = 758;

for iROI=1:nROI
   
    nClusters = size(FC_densities_clusters(iROI).track_positive,2);
    
    for iCluster=1:nClusters
        
       track_pos_sample = squeeze(FC_densities_clusters(iROI).track_positive(:,iCluster));
       track_neg_sample = squeeze(FC_densities_clusters(iROI).track_negative(:,iCluster));
       
       passive_pos_sample = squeeze(FC_densities_clusters(iROI).passive_positive(:,iCluster));
       passive_neg_sample = squeeze(FC_densities_clusters(iROI).passive_negative(:,iCluster));
       
       rest_pos_sample = squeeze(FC_densities_clusters(iROI).rest_positive(:,iCluster));
       rest_neg_sample = squeeze(FC_densities_clusters(iROI).rest_negative(:,iCluster));
       
       [H, p_pos_stimulus] = ttest(rest_pos_sample(:),passive_pos_sample(:));
       if (mean(passive_pos_sample(:)) - mean(rest_pos_sample(:))) > 0; z_pos_stimulus = +1 * norminv(p_pos_stimulus / nTotalClustersBonferroni); else z_pos_stimulus = -1 * norminv(p_pos_stimulus / nTotalClustersBonferroni); end
       
       [H, p_neg_stimulus] = ttest(rest_neg_sample(:),passive_neg_sample(:));
       if (mean(passive_neg_sample(:)) - mean(rest_neg_sample(:))) > 0; z_neg_stimulus = +1 * norminv(p_neg_stimulus / nTotalClustersBonferroni); else z_neg_stimulus = -1 * norminv(p_neg_stimulus / nTotalClustersBonferroni); end
       
       [H, p_pos_attention] = ttest(passive_pos_sample(:),track_pos_sample(:));
       if (mean(track_pos_sample(:)) - mean(passive_pos_sample(:))) > 0; z_pos_attention = +1 * norminv(p_pos_attention / nTotalClustersBonferroni); else z_pos_attention = -1 * norminv(p_pos_attention / nTotalClustersBonferroni); end
       
       [H, p_neg_attention] = ttest(passive_neg_sample(:),track_neg_sample(:));
       if (mean(track_neg_sample(:)) - mean(passive_neg_sample(:))) > 0; z_neg_attention = +1 * norminv(p_neg_attention / nTotalClustersBonferroni); else z_neg_attention = -1 * norminv(p_neg_attention / nTotalClustersBonferroni); end
      
       idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
       
       nVoxels = length(idx_voxels);
       
       for iVoxel=1:nVoxels
           
          [x,y,z] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
          
          img_pos_stimulus(x,y,z) = z_pos_stimulus;
          img_neg_stimulus(x,y,z) = z_neg_stimulus;
          
          img_pos_attention(x,y,z) = z_pos_attention;
          img_neg_attention(x,y,z) = z_neg_attention;
           
       end
       
    end
    
end

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;
dtype = 'FLOAT32';
offset = 0;
dim = load_aal.dat.dim;
descrip = 'img_pos_stimulus';
fname = strcat('LHR','-','All-Subjects','-','FC-Voxel-AAL-ROI-KMeans-Densities-Pos-Stimulus','-',label,'.nii');
input_data = img_pos_stimulus; 
real_save_image;

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;
dtype = 'FLOAT32';
offset = 0;
dim = load_aal.dat.dim;
descrip = 'img_neg_stimulus';
fname = strcat('LHR','-','All-Subjects','-','FC-Voxel-AAL-ROI-KMeans-Densities-Neg-Stimulus','-',label,'.nii');
input_data = img_neg_stimulus; 
real_save_image;

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;
dtype = 'FLOAT32';
offset = 0;
dim = load_aal.dat.dim;
descrip = 'img_pos_attention';
fname = strcat('LHR','-','All-Subjects','-','FC-Voxel-AAL-ROI-KMeans-Densities-Pos-Attention','-',label,'.nii');
input_data = img_pos_attention; 
real_save_image;

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;
dtype = 'FLOAT32';
offset = 0;
dim = load_aal.dat.dim;
descrip = 'img_neg_attention';
fname = strcat('LHR','-','All-Subjects','-','FC-Voxel-AAL-ROI-KMeans-Densities-Neg-Attention','-',label,'.nii');
input_data = img_neg_attention; 
real_save_image;

end

function doFCOnAllClusters

load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nROI = 90;
nRuns = 32;
nTR = 150;
nTotalClusters = 758;

for iRun=1:nRuns

    iiCluster = 0;
    
    track_ts = zeros(nTR,nTotalClusters);
    passive_ts = zeros(nTR,nTotalClusters);
    rest_ts = zeros(nTR,nTotalClusters);
    
    for iROI=1:nROI

        nClusters = length(ROI(iROI).clusters);
        
        for iCluster=1:nClusters
           
            iiCluster = iiCluster + 1;
            
            track_ts(:,iiCluster) = squeeze(ROI(iROI).clusters(iCluster).track_mean(iRun,:));
            passive_ts(:,iiCluster) = squeeze(ROI(iROI).clusters(iCluster).passive_mean(iRun,:));
            rest_ts(:,iiCluster) = squeeze(ROI(iROI).clusters(iCluster).rest_mean(iRun,:));
            
        end

    end
    
    [FC_clusters(iRun).track_rho, FC_clusters(iRun).track_pval] = corr(track_ts);
    [FC_clusters(iRun).passive_rho, FC_clusters(iRun).passive_pval] = corr(passive_ts);
    [FC_clusters(iRun).rest_rho, FC_clusters(iRun).rest_pval] = corr(rest_ts);

end

save('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters.mat','FC_clusters');

end

function plotFCOnAllClustersMeanRuns

load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nRun = 32;
nTotalClusters = 758;
nTR = 150;
pcriterion = 0.01;
nROI = 90;

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

x_idx(1) = 29;
x_idx(2) = 43;
x_idx(3) = 55;
x_idx(4) = 56;
x_idx(5) = 68;
x_idx(6) = 70;
x_idx(7) = 78;
x_idx(8) = 90;

x_label{1} = 'frontal';
x_label{2} = 'subcortical';
x_label{3} = 'occipital';
x_label{4} = 'temporal';
x_label{5} = 'parietal';
x_label{6} = 'frontal';
x_label{7} = 'subcortical';
x_label{8} = 'temporal';

mean_track = zeros(nTotalClusters,nTotalClusters);
mean_passive = zeros(nTotalClusters,nTotalClusters);
mean_rest = zeros(nTotalClusters,nTotalClusters);

for iRun=1:nRun

    mean_track = mean_track + FC_clusters(iRun).track_rho;
    mean_passive = mean_passive + FC_clusters(iRun).passive_rho;
    mean_rest = mean_rest + FC_clusters(iRun).rest_rho;
    
end

mean_track = mean_track ./ nRun;
mean_passive = mean_passive ./ nRun;
mean_rest = mean_rest ./ nRun;

mean_track_pval = erfc((abs(mean_track).*sqrt(nTR))./sqrt(2));
mean_passive_pval = erfc((abs(mean_passive).*sqrt(nTR))./sqrt(2));
mean_rest_pval = erfc((abs(mean_rest).*sqrt(nTR))./sqrt(2));

mean_track(mean_track_pval > pcriterion) = 0;
mean_passive(mean_passive_pval > pcriterion) = 0;
mean_rest(mean_rest_pval > pcriterion) = 0;


for iCondition=1:3
    
    f = figure;
    
    if iCondition == 1; rho = mean_track; label = 'Track'; end
    if iCondition == 2; rho = mean_passive; label = 'Passive'; end
    if iCondition == 3; rho = mean_rest; label = 'Rest'; end

    imagesc(rho);
    hold on

    title(label);

    clrmp = colormap('jet');
    clrmp(1,:) = [1 1 1];

    min_C = 0;
    max_C = 1;

    colorbar;
    caxis([min_C max_C]);
    %colorbar('Ticks',[min_C max_C/2 max_C]);
    colormap(clrmp);
    
    jump = 0;
    last_clusters_so_far = 0;
    for iROI=1:nROI

       if ~isempty(find(ismember(x_idx,iROI)))

           jump = jump + 1;

           if jump == 1; roi_distance = 0; else roi_distance = round( x_idx(jump) - x_idx(jump-1) )/2; end

           if jump == 1; limit_roi_clusters = round(x_idx(jump)/2); else limit_roi_clusters = x_idx(jump-1); end

           clusters_so_far = 0;
           for iiROI=1:(limit_roi_clusters+roi_distance)

               clusters_so_far = clusters_so_far + ROI(iiROI).nClusters;

           end

           x_idx_jump(jump) = clusters_so_far;

       end

    end

    clusters_so_far = 0;
    nClusters = 758;
    for iROI=1:nROI

       clusters_so_far = ROI(iROI).nClusters + clusters_so_far;

       if ~isempty(find(ismember(x_idx,iROI)))

           hold on;

           plot([clusters_so_far clusters_so_far],[1 nClusters],'k-');

           plot([1 nClusters],[clusters_so_far clusters_so_far],'k-');

       end

    end

    ax = gca;
    %ax.XTick = x_idx_jump;
    %ax.XTickLabel = x_label;
    set(ax,'XTick',x_idx_jump);
    set(ax,'XTickLabel',x_label);

    %ax.YTick = x_idx_jump;
    %ax.YTickLabel = x_label;
    set(ax,'YTick',x_idx_jump);
    set(ax,'YTickLabel',x_label);

    xticklabel_rotate;

    print(f,'-depsc',strcat('FC-All-Clusters-Mean-Runs-',label,'.eps'));

end

end

function plotFCOnAllClustersPerRun

load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nRun = 32;
nTotalClusters = 758;
nTR = 150;
pcriterion = 0.01;
nROI = 90;

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

x_idx(1) = 29;
x_idx(2) = 43;
x_idx(3) = 55;
x_idx(4) = 56;
x_idx(5) = 68;
x_idx(6) = 70;
x_idx(7) = 78;
x_idx(8) = 90;

x_label{1} = 'frontal';
x_label{2} = 'subcortical';
x_label{3} = 'occipital';
x_label{4} = 'temporal';
x_label{5} = 'parietal';
x_label{6} = 'frontal';
x_label{7} = 'subcortical';
x_label{8} = 'temporal';

for iRun=1:nRun

    track_rho = FC_clusters(iRun).track_rho;
    passive_rho = FC_clusters(iRun).passive_rho;
    rest_rho = FC_clusters(iRun).rest_rho;
    
    track_pval = FC_clusters(iRun).track_pval;
    passive_pval = FC_clusters(iRun).passive_pval;
    rest_pval = FC_clusters(iRun).rest_pval;

    track_rho(track_pval > pcriterion) = 0;
    passive_rho(passive_pval > pcriterion) = 0;
    rest_rho(rest_pval > pcriterion) = 0;


for iCondition=1:3
    
    f = figure;
    
    if iCondition == 1; rho = track_rho; label = strcat('Track-',int2str(iRun)); end
    if iCondition == 2; rho = passive_rho; label = strcat('Passive-',int2str(iRun)); end
    if iCondition == 3; rho = rest_rho; label = strcat('Rest-',int2str(iRun)); end

    imagesc(rho);
    hold on

    title(label);

    clrmp = colormap('jet');
    clrmp(33,:) = [1 1 1];
    clrmp(65,:) = clrmp(64,:);

    min_C = -1;
    max_C = 1;

    %colorbar;
    caxis([min_C max_C]);
    %colorbar('Ticks',[min_C max_C/2 max_C]);
    colormap(clrmp);
    
    jump = 0;
    last_clusters_so_far = 0;
    for iROI=1:nROI

       if ~isempty(find(ismember(x_idx,iROI)))

           jump = jump + 1;

           if jump == 1; roi_distance = 0; else roi_distance = round( x_idx(jump) - x_idx(jump-1) )/2; end

           if jump == 1; limit_roi_clusters = round(x_idx(jump)/2); else limit_roi_clusters = x_idx(jump-1); end

           clusters_so_far = 0;
           for iiROI=1:(limit_roi_clusters+roi_distance)

               clusters_so_far = clusters_so_far + ROI(iiROI).nClusters;

           end

           x_idx_jump(jump) = clusters_so_far;

       end

    end

    clusters_so_far = 0;
    nClusters = 758;
    for iROI=1:nROI

       clusters_so_far = ROI(iROI).nClusters + clusters_so_far;

       if ~isempty(find(ismember(x_idx,iROI)))

           hold on;

           plot([clusters_so_far clusters_so_far],[1 nClusters],'k-');

           plot([1 nClusters],[clusters_so_far clusters_so_far],'k-');

       end

    end

    ax = gca;
    %ax.XTick = x_idx_jump;
    %ax.XTickLabel = x_label;
    set(ax,'XTick',x_idx_jump);
    set(ax,'XTickLabel',x_label);

    %ax.YTick = x_idx_jump;
    %ax.YTickLabel = x_label;
    set(ax,'YTick',x_idx_jump);
    set(ax,'YTickLabel',x_label);

    xticklabel_rotate;

    print(f,'-depsc',strcat('FC-All-Clusters-Per-Run','-',label,'.eps'));

end

end

end

function doContrastsFCOnAllClusters

load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nRuns = 32;
nClusters = 758;
pcriterion = 0.01;

nROI = 90;

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

x_idx(1) = 29;
x_idx(2) = 43;
x_idx(3) = 55;
x_idx(4) = 56;
x_idx(5) = 68;
x_idx(6) = 70;
x_idx(7) = 78;
x_idx(8) = 90;

x_label{1} = 'frontal';
x_label{2} = 'subcortical';
x_label{3} = 'occipital';
x_label{4} = 'temporal';
x_label{5} = 'parietal';
x_label{6} = 'frontal';
x_label{7} = 'subcortical';
x_label{8} = 'temporal';

track_corr = zeros(nRuns,nClusters*nClusters);
passive_corr = zeros(nRuns,nClusters*nClusters);
rest_corr = zeros(nRuns,nClusters*nClusters);

for iRun=1:nRuns
    
    track_corr(iRun,:) = squeeze(FC_clusters(iRun).track_rho(:));
    passive_corr(iRun,:) = squeeze(FC_clusters(iRun).passive_rho(:));
    rest_corr(iRun,:) = squeeze(FC_clusters(iRun).rest_rho(:));
    
end

p_attention = zeros(1,nClusters*nClusters);
h_attention = zeros(1,nClusters*nClusters);
p_stimulus = zeros(1,nClusters*nClusters);
h_stimulus = zeros(1,nClusters*nClusters);

%% FISHER TRANSFORM
% track_corr_f = (1/2) * log( (1 + track_corr) ./ (1 - track_corr) );
% passive_corr_f = (1/2) * log( (1 + passive_corr) ./ (1 - passive_corr) );
% rest_corr_f = (1/2) * log( (1 + rest_corr) ./ (1 - rest_corr) );

for iCorr=1:(nClusters*nClusters)
    
   if mod(iCorr,1000) == 0; disp(int2str(iCorr)); end
    
   p_attention(iCorr) = ranksum(squeeze(track_corr(:,iCorr)),squeeze(passive_corr(:,iCorr)));
   %[h, p_attention(iCorr)] = ttest(squeeze(track_corr_f(:,iCorr)),squeeze(passive_corr_f(:,iCorr)));
   
   if median(squeeze(track_corr(:,iCorr))) > median(squeeze(passive_corr(:,iCorr))); h_attention(iCorr) = -1; end
   %if mean(squeeze(track_corr_f(:,iCorr))) > mean(squeeze(passive_corr_f(:,iCorr))); h_attention(iCorr) = 1; end
   if median(squeeze(track_corr(:,iCorr))) < median(squeeze(passive_corr(:,iCorr))); h_attention(iCorr) = 1; end
   %if mean(squeeze(track_corr_f(:,iCorr))) < mean(squeeze(passive_corr_f(:,iCorr))); h_attention(iCorr) = -1; end
    
   p_stimulus(iCorr) = ranksum(squeeze(passive_corr(:,iCorr)),squeeze(rest_corr(:,iCorr)));
   %[h, p_stimulus(iCorr)] = ttest(squeeze(passive_corr_f(:,iCorr)),squeeze(rest_corr_f(:,iCorr)));
   
   if median(squeeze(passive_corr(:,iCorr))) > median(squeeze(rest_corr(:,iCorr))); h_stimulus(iCorr) = -1; end
   %if mean(squeeze(passive_corr_f(:,iCorr))) > mean(squeeze(rest_corr_f(:,iCorr))); h_stimulus(iCorr) = 1; end
  if median(squeeze(passive_corr(:,iCorr))) < median(squeeze(rest_corr(:,iCorr))); h_stimulus(iCorr) = 1; end
   %if mean(squeeze(passive_corr_f(:,iCorr))) < mean(squeeze(rest_corr_f(:,iCorr))); h_stimulus(iCorr) = -1; end

end

contrast_attention = zeros(nClusters);
contrast_stimulus = zeros(nClusters);

contrast_attention_z = zeros(nClusters);
contrast_stimulus_z = zeros(nClusters);

for iCorr=1:(nClusters*nClusters)
   
    [x,y] = ind2sub(size(contrast_attention),iCorr);
    
    contrast_attention(x,y) = p_attention(iCorr);
    contrast_stimulus(x,y) = p_stimulus(iCorr);
    
end

contrast_attention(contrast_attention>pcriterion) = 0;
contrast_stimulus(contrast_stimulus>pcriterion) = 0;

for iCorr=1:(nClusters*nClusters)
   
    [x,y] = ind2sub(size(contrast_attention),iCorr);
    
    if contrast_attention(x,y) ~= 0; contrast_attention_z(x,y) = norminv(contrast_attention(x,y)); end
    if contrast_stimulus(x,y) ~= 0; contrast_stimulus_z(x,y) = norminv(contrast_stimulus(x,y)); end

    contrast_attention_z(x,y) = h_attention(iCorr) * contrast_attention_z(x,y);
    contrast_stimulus_z(x,y) = h_stimulus(iCorr) * contrast_stimulus_z(x,y);
    
end

save('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-Contrast.mat','contrast_stimulus_z','contrast_attention_z');

end

function plotContrastFCOnAllClusters

load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-Contrast.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nRuns = 32;
nClusters = 758;
pcriterion = 0.01;

nROI = 90;

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

x_idx(1) = 29;
x_idx(2) = 43;
x_idx(3) = 55;
x_idx(4) = 56;
x_idx(5) = 68;
x_idx(6) = 70;
x_idx(7) = 78;
x_idx(8) = 90;

x_label{1} = 'frontal';
x_label{2} = 'subcortical';
x_label{3} = 'occipital';
x_label{4} = 'temporal';
x_label{5} = 'parietal';
x_label{6} = 'frontal';
x_label{7} = 'subcortical';
x_label{8} = 'temporal';

for imatrix=1:2
    
    if imatrix == 1; rho = contrast_attention_z; label = 'Attention'; end
    if imatrix == 2; rho = contrast_stimulus_z; label = 'Stimulus'; end

    %rho(rho<0) = 0;

    f = figure;
    rho(isnan(rho)) = 0;
    imagesc(rho);
    hold on

    title(label);

    clrmp = colormap('jet');
    clrmp(33,:) = [1 1 1];
    clrmp(65,:) = clrmp(64,:);

    min_C = -4.5;
    max_C = 4.5;

    caxis([min_C max_C]);
    %colorbar('Ticks',[min_C max_C/2 max_C]);
    colormap(clrmp);
    colorbar;
    
    jump = 0;
    last_clusters_so_far = 0;
    for iROI=1:nROI

       if ~isempty(find(ismember(x_idx,iROI)))

           jump = jump + 1;

           if jump == 1; roi_distance = 0; else roi_distance = round( x_idx(jump) - x_idx(jump-1) )/2; end

           if jump == 1; limit_roi_clusters = round(x_idx(jump)/2); else limit_roi_clusters = x_idx(jump-1); end

           clusters_so_far = 0;
           for iiROI=1:(limit_roi_clusters+roi_distance)

               clusters_so_far = clusters_so_far + ROI(iiROI).nClusters;

           end

           x_idx_jump(jump) = clusters_so_far;

       end

    end

    clusters_so_far = 0;
    nClusters = 758;
    for iROI=1:nROI

       clusters_so_far = ROI(iROI).nClusters + clusters_so_far;

       if ~isempty(find(ismember(x_idx,iROI)))

           hold on;

           plot([clusters_so_far clusters_so_far],[1 nClusters],'k-');

           plot([1 nClusters],[clusters_so_far clusters_so_far],'k-');

       end

    end

    ax = gca;
    %ax.XTick = x_idx_jump;
    %ax.XTickLabel = x_label;
    set(ax,'XTick',x_idx_jump);
    set(ax,'XTickLabel',x_label);

    %ax.YTick = x_idx_jump;
    %ax.YTickLabel = x_label;
    set(ax,'YTick',x_idx_jump);
    set(ax,'YTickLabel',x_label);

    xticklabel_rotate;

    print(f,'-depsc',strcat('FC-All-Clusters-Contrast','-',label,'.eps'));

end


end

function computeGrangerAmongClusters

load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat');

nClusters = 758;
nROI = 90;
nRuns = 32;
nTRs = 150;

X_track = zeros(nClusters,nTRs,nRuns);
X_passive = zeros(nClusters,nTRs,nRuns);
X_rest = zeros(nClusters,nTRs,nRuns);

iiCluster = 0;
for iROI=1:nROI
    
    nnClusters = ROI(iROI).nClusters;
    
    for iCluster=1:nnClusters
        
        iiCluster = iiCluster + 1;
        
        for iRun=1:nRuns
           
            X_track(iiCluster,:,iRun) = squeeze(ROI(iROI).clusters(iCluster).track_mean(iRun,:));
            X_passive(iiCluster,:,iRun) = squeeze(ROI(iROI).clusters(iCluster).passive_mean(iRun,:));
            X_rest(iiCluster,:,iRun) = squeeze(ROI(iROI).clusters(iCluster).rest_mean(iRun,:));
            
        end
        
    end

end

for iRun=1:nRuns

    disp(strcat('Run:',int2str(iRun)));
    
    track_GF = zeros(nClusters,nClusters);
    track_Gpval = zeros(nClusters,nClusters);
    passive_GF = zeros(nClusters,nClusters);
    passive_Gpval = zeros(nClusters,nClusters);
    rest_GF = zeros(nClusters,nClusters);
    rest_Gpval = zeros(nClusters,nClusters);
    
    for iCluster=1:(nClusters-1)
    
        disp(strcat('Cluster:',int2str(iCluster)));
        
        for iiCluster=(iCluster+1):nClusters
        
            %disp('Track');
            
            [GF, Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(X_track(iCluster,:,iRun),X_track(iiCluster,:,iRun));
            track_GF(iCluster,iiCluster) = GF(1,2);
            track_GF(iiCluster,iCluster) = GF(2,1);
            track_Gpval(iCluster,iiCluster) = Gpval(1,2);
            track_Gpval(iiCluster,iCluster) = Gpval(2,1);
            
            %disp('Passive');
            
            [GF, Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(X_passive(iCluster,:,iRun),X_passive(iiCluster,:,iRun));
            passive_GF(iCluster,iiCluster) = GF(1,2);
            passive_GF(iiCluster,iCluster) = GF(2,1);
            passive_Gpval(iCluster,iiCluster) = Gpval(1,2);
            passive_Gpval(iiCluster,iCluster) = Gpval(2,1);
            
            %disp('Rest');
            
            [GF, Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(X_rest(iCluster,:,iRun),X_rest(iiCluster,:,iRun));
            rest_GF(iCluster,iiCluster) = GF(1,2);
            rest_GF(iiCluster,iCluster) = GF(2,1);
            rest_Gpval(iCluster,iiCluster) = Gpval(1,2);
            rest_Gpval(iiCluster,iCluster) = Gpval(2,1);
            
        end
        
    end
    
    save(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Run','-',int2str(iRun),'.mat'),'track_GF','track_Gpval','passive_GF','passive_Gpval','rest_GF','rest_Gpval');
    
end

%[G_track.GF, G_track.Gpval, G_track.GSig, G_track.morder, G_track.A, G_track.G, G_track.info] = getGrangerX(X_track);

end

function computeGrangerAtGroupLevel

nRuns = 32;
nClusters = 758;

track_mean = zeros(nClusters,nClusters);
passive_mean = zeros(nClusters,nClusters);
rest_mean = zeros(nClusters,nClusters);

for iRun=1:nRuns
    
    load(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Run','-',int2str(iRun),'.mat'));
    
    tracking(iRun).pval = track_Gpval;
    passive(iRun).pval = passive_Gpval;
    rest(iRun).pval = rest_Gpval;
    
    tracking(iRun).zscore = -1 * norminv(track_Gpval);
    passive(iRun).zscore = -1 * norminv(passive_Gpval);
    rest(iRun).zscore = -1 * norminv(rest_Gpval);
    
    track_mean = tracking(iRun).zscore + track_mean;
    passive_mean = passive(iRun).zscore + passive_mean;
    rest_mean = rest(iRun).zscore + rest_mean;
    
    clear rest_GF
    clear rest_Gpval
    clear passive_GF
    clear passive_Gpval
    clear track_GF
    clear track_Gpval
    
end

track_mean(isinf(track_mean)) = 0;
track_mean(isnan(track_mean)) = 0;

passive_mean(isinf(passive_mean)) = 0;
passive_mean(isnan(passive_mean)) = 0;

rest_mean(isinf(rest_mean)) = 0;
rest_mean(isnan(rest_mean)) = 0;

track_mean = track_mean ./ nClusters;
passive_mean = passive_mean ./ nClusters;
rest_mean = rest_mean ./ nClusters;

for iCluster=1:nClusters
   
    for iiCluster=1:nClusters
       
        track_samples = zeros(1,nRuns);
        passive_samples = zeros(1,nRuns);
        rest_samples = zeros(1,nRuns);
        
        for iRun=1:nRuns
           
            track_samples(iRun) = tracking(iRun).zscore(iCluster,iiCluster);
            passive_samples(iRun) = passive(iRun).zscore(iCluster,iiCluster);
            rest_samples(iRun) = rest(iRun).zscore(iCluster,iiCluster);
            
        end
        
        [H, P] = ttest(track_samples,passive_samples);
        
        if mean(track_samples) > mean(passive_samples); signal = -1; else signal = 1; end
        
        Attention_Contrast.P(iCluster,iiCluster) = P;
        Attention_Contrast.Z(iCluster,iiCluster) = signal * norminv(P);
        
        [H, P] = ttest(passive_samples,rest_samples);
        
        if mean(passive_samples) > mean(rest_samples); signal = -1; else signal = 1; end
        
        Stimulus_Contrast.P(iCluster,iiCluster) = P;
        Stimulus_Contrast.Z(iCluster,iiCluster) = signal * norminv(P);
        
    end
    
end

save(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast','.mat'),'Attention_Contrast','Stimulus_Contrast','track_mean','passive_mean','rest_mean','tracking','passive','rest');

end

function computeGrangerAtGroupLevelCorbettasMethod

nRuns = 32;
nClusters = 758;

track_mean = zeros(nClusters,nClusters);
passive_mean = zeros(nClusters,nClusters);
rest_mean = zeros(nClusters,nClusters);

for iRun=1:nRuns
    
    load(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Run','-',int2str(iRun),'.mat'));
    
    tracking(iRun).pval = track_Gpval;
    passive(iRun).pval = passive_Gpval;
    rest(iRun).pval = rest_Gpval;
    
    tracking(iRun).zscore = -1 * norminv(track_Gpval);
    passive(iRun).zscore = -1 * norminv(passive_Gpval);
    rest(iRun).zscore = -1 * norminv(rest_Gpval);
    
    track_mean = tracking(iRun).zscore + track_mean;
    passive_mean = passive(iRun).zscore + passive_mean;
    rest_mean = rest(iRun).zscore + rest_mean;
    
    clear rest_GF
    clear rest_Gpval
    clear passive_GF
    clear passive_Gpval
    clear track_GF
    clear track_Gpval
    
end

track_mean(isinf(track_mean)) = 0;
track_mean(isnan(track_mean)) = 0;

passive_mean(isinf(passive_mean)) = 0;
passive_mean(isnan(passive_mean)) = 0;

rest_mean(isinf(rest_mean)) = 0;
rest_mean(isnan(rest_mean)) = 0;

track_mean = track_mean ./ nClusters;
passive_mean = passive_mean ./ nClusters;
rest_mean = rest_mean ./ nClusters;

for iCluster=1:nClusters
   
    for iiCluster=1:nClusters
       
        track_samples = zeros(1,nRuns);
        passive_samples = zeros(1,nRuns);
        rest_samples = zeros(1,nRuns);
        
        for iRun=1:nRuns
           
            track_samples(iRun) = tracking(iRun).zscore(iCluster,iiCluster);
            passive_samples(iRun) = passive(iRun).zscore(iCluster,iiCluster);
            rest_samples(iRun) = rest(iRun).zscore(iCluster,iiCluster);
            
        end
        
        [H, P] = ttest(track_samples,passive_samples);
        
        if mean(track_samples) > mean(passive_samples); signal = -1; else signal = 1; end
        
        Attention_Contrast.P(iCluster,iiCluster) = P;
        Attention_Contrast.Z(iCluster,iiCluster) = signal * norminv(P);
        
        [H, P] = ttest(passive_samples,rest_samples);
        
        if mean(passive_samples) > mean(rest_samples); signal = -1; else signal = 1; end
        
        Stimulus_Contrast.P(iCluster,iiCluster) = P;
        Stimulus_Contrast.Z(iCluster,iiCluster) = signal * norminv(P);
        
    end
    
end

save(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast','.mat'),'Attention_Contrast','Stimulus_Contrast','track_mean','passive_mean','rest_mean','tracking','passive','rest');

end

function plotGrangerAtGroupLevelMean

load(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast','.mat'));
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nRun = 32;
nTotalClusters = 758;
nTR = 150;
pcriterion = 0.01;
nROI = 90;

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

x_idx(1) = 29;
x_idx(2) = 43;
x_idx(3) = 55;
x_idx(4) = 56;
x_idx(5) = 68;
x_idx(6) = 70;
x_idx(7) = 78;
x_idx(8) = 90;

x_label{1} = 'frontal';
x_label{2} = 'subcortical';
x_label{3} = 'occipital';
x_label{4} = 'temporal';
x_label{5} = 'parietal';
x_label{6} = 'frontal';
x_label{7} = 'subcortical';
x_label{8} = 'temporal';

all_rhos = [track_mean(:), passive_mean(:), rest_mean(:)];

for iCondition=1:3
    
    f = figure;
    
    if iCondition == 1; rho = track_mean; label = 'Track'; end
    if iCondition == 2; rho = passive_mean; label = 'Passive'; end
    if iCondition == 3; rho = rest_mean; label = 'Rest'; end
    
    imagesc(rho);
    hold on

    title(label);

    clrmp = colormap('jet');
    clrmp(1,:) = [1 1 1];

    min_C = 0;
    max_C = max(all_rhos(:));

    colorbar;
    caxis([min_C max_C]);
    %colorbar('Ticks',[min_C max_C/2 max_C]);
    colormap(clrmp);
    
    jump = 0;
    last_clusters_so_far = 0;
    for iROI=1:nROI

       if ~isempty(find(ismember(x_idx,iROI)))

           jump = jump + 1;

           if jump == 1; roi_distance = 0; else roi_distance = round( x_idx(jump) - x_idx(jump-1) )/2; end

           if jump == 1; limit_roi_clusters = round(x_idx(jump)/2); else limit_roi_clusters = x_idx(jump-1); end

           clusters_so_far = 0;
           for iiROI=1:(limit_roi_clusters+roi_distance)

               clusters_so_far = clusters_so_far + ROI(iiROI).nClusters;

           end

           x_idx_jump(jump) = clusters_so_far;

       end

    end

    clusters_so_far = 0;
    nClusters = 758;
    for iROI=1:nROI

       clusters_so_far = ROI(iROI).nClusters + clusters_so_far;

       if ~isempty(find(ismember(x_idx,iROI)))

           hold on;

           plot([clusters_so_far clusters_so_far],[1 nClusters],'k-');

           plot([1 nClusters],[clusters_so_far clusters_so_far],'k-');

       end

    end

    ax = gca;
    %ax.XTick = x_idx_jump;
    %ax.XTickLabel = x_label;
    set(ax,'XTick',x_idx_jump);
    set(ax,'XTickLabel',x_label);

    %ax.YTick = x_idx_jump;
    %ax.YTickLabel = x_label;
    set(ax,'YTick',x_idx_jump);
    set(ax,'YTickLabel',x_label);

    xticklabel_rotate;

    print(f,'-depsc',strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast-',label,'.eps'));

end

end

function plotGrangerAtGroupLevelContrast

load(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast','.mat'));
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nRun = 32;
nTotalClusters = 758;
nTR = 150;
pcriterion = 0.01;
nROI = 90;

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

x_idx(1) = 29;
x_idx(2) = 43;
x_idx(3) = 55;
x_idx(4) = 56;
x_idx(5) = 68;
x_idx(6) = 70;
x_idx(7) = 78;
x_idx(8) = 90;

x_label{1} = 'frontal';
x_label{2} = 'subcortical';
x_label{3} = 'occipital';
x_label{4} = 'temporal';
x_label{5} = 'parietal';
x_label{6} = 'frontal';
x_label{7} = 'subcortical';
x_label{8} = 'temporal';

%zcriterion = 1.6;
zcriterion = 2.3;

for iCondition=1:2
    
    f = figure;
    
    if iCondition == 1; rho = Attention_Contrast.Z; label = 'Attention'; end
    if iCondition == 2; rho = Stimulus_Contrast.Z; label = 'Stimulus'; end
    
    rho(find(rho>(-1)*zcriterion & rho<zcriterion)) = 0; 
    
    imagesc(rho);
    hold on
   
    title(label);

    clrmp = colormap('jet');
    clrmp(65,:) = clrmp(64,:);
    clrmp(33,:) = [1 1 1];

    min_C = -3;
    max_C = 3;

    %colorbar;
    caxis([min_C max_C]);
    %colorbar('Ticks',[min_C max_C/2 max_C]);
    colormap(clrmp);
    colorbar;
    
    jump = 0;
    last_clusters_so_far = 0;
    for iROI=1:nROI

       if ~isempty(find(ismember(x_idx,iROI)))

           jump = jump + 1;

           if jump == 1; roi_distance = 0; else roi_distance = round( x_idx(jump) - x_idx(jump-1) )/2; end

           if jump == 1; limit_roi_clusters = round(x_idx(jump)/2); else limit_roi_clusters = x_idx(jump-1); end

           clusters_so_far = 0;
           for iiROI=1:(limit_roi_clusters+roi_distance)

               clusters_so_far = clusters_so_far + ROI(iiROI).nClusters;

           end

           x_idx_jump(jump) = clusters_so_far;

       end

    end

    clusters_so_far = 0;
    nClusters = 758;
    for iROI=1:nROI

       clusters_so_far = ROI(iROI).nClusters + clusters_so_far;

       if ~isempty(find(ismember(x_idx,iROI)))

           hold on;

           plot([clusters_so_far clusters_so_far],[1 nClusters],'k-');

           plot([1 nClusters],[clusters_so_far clusters_so_far],'k-');

       end

    end

    ax = gca;
    %ax.XTick = x_idx_jump;
    %ax.XTickLabel = x_label;
    set(ax,'XTick',x_idx_jump);
    set(ax,'XTickLabel',x_label);

    %ax.YTick = x_idx_jump;
    %ax.YTickLabel = x_label;
    set(ax,'YTick',x_idx_jump);
    set(ax,'YTickLabel',x_label);

    xticklabel_rotate;

    print(f,'-depsc',strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast-',label,'-z-',int2str(zcriterion),'.eps'));

end

end

function plotGrangerAtGroupLevelContrastDTI

load(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast','.mat'));
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nRun = 32;
nTotalClusters = 758;
nTR = 150;
pcriterion = 0.01;
nROI = 90;

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

x_idx(1) = 29;
x_idx(2) = 43;
x_idx(3) = 55;
x_idx(4) = 56;
x_idx(5) = 68;
x_idx(6) = 70;
x_idx(7) = 78;
x_idx(8) = 90;

x_label{1} = 'frontal';
x_label{2} = 'subcortical';
x_label{3} = 'occipital';
x_label{4} = 'temporal';
x_label{5} = 'parietal';
x_label{6} = 'frontal';
x_label{7} = 'subcortical';
x_label{8} = 'temporal';

prefix = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION';
sufix = 'preprocessed\B0-DTI\6.forbedpost.bedpostX\Und_PRE_xx.mat';

DTI(1) = load(strcat(prefix,'\','SUBJECT-1-22-10-2015','\',sufix));
DTI(2) = load(strcat(prefix,'\','SUBJECT-2-26-10-2015','\',sufix));
DTI(3) = load(strcat(prefix,'\','SUBJECT-3-3-11-2015','\',sufix));
DTI(4) = load(strcat(prefix,'\','SUBJECT-4-2-11-2015','\',sufix));
DTI(5) = load(strcat(prefix,'\','SUBJECT-5-2-11-2015','\',sufix));
DTI(6) = load(strcat(prefix,'\','SUBJECT-6-24-11-2015','\',sufix));
DTI(7) = load(strcat(prefix,'\','SUBJECT-7-14-01-2016','\',sufix));
DTI(8) = load(strcat(prefix,'\','SUBJECT-8-14-01-2016','\',sufix));

ave_DTI = zeros(size(DTI(1).C));
for iDTI=1:8
    
    ave_DTI = ave_DTI + DTI(iDTI).C;
    
end
ave_DTI = ave_DTI ./ 8;

common_DTI = ave_DTI;
common_DTI(~(DTI(1).C & DTI(2).C & DTI(3).C & DTI(4).C & DTI(5).C & DTI(6).C & DTI(7).C & DTI(8).C)) = 0;

zcriterion = 1.6;

for iCondition=1:2
    
    if iCondition == 1; rho = Attention_Contrast.Z; label = 'Attention'; end
    if iCondition == 2; rho = Stimulus_Contrast.Z; label = 'Stimulus'; end
    
    rho(isnan(rho)) = 0;
    rho(isinf(rho)) = 0;
    
    for iDirection=1:2
        
        f = figure;
        
        rho(find(rho>(-1)*zcriterion & rho<zcriterion)) = 0; 

        rho_direct = rho;
        rho_direct(~common_DTI) = 0;

        rho_indirect = rho;
        rho_indirect(find(common_DTI)) = 0;
        
        if iDirection == 1; matrix = rho_direct; sub_label = 'direct'; end
        if iDirection == 2; matrix = rho_indirect; sub_label = 'indirect'; end

        imagesc(matrix);
        hold on

        title(strcat(label,'-',sub_label));

        clrmp = colormap('jet');
        clrmp(65,:) = clrmp(64,:);
        clrmp(33,:) = [1 1 1];

        min_C = -3;
        max_C = 3;

        %colorbar;
        caxis([min_C max_C]);
        %colorbar('Ticks',[min_C max_C/2 max_C]);
        colormap(clrmp);
        colorbar;

        jump = 0;
        last_clusters_so_far = 0;
        for iROI=1:nROI

           if ~isempty(find(ismember(x_idx,iROI)))

               jump = jump + 1;

               if jump == 1; roi_distance = 0; else roi_distance = round( x_idx(jump) - x_idx(jump-1) )/2; end

               if jump == 1; limit_roi_clusters = round(x_idx(jump)/2); else limit_roi_clusters = x_idx(jump-1); end

               clusters_so_far = 0;
               for iiROI=1:(limit_roi_clusters+roi_distance)

                   clusters_so_far = clusters_so_far + ROI(iiROI).nClusters;

               end

               x_idx_jump(jump) = clusters_so_far;

           end

        end

        clusters_so_far = 0;
        nClusters = 758;
        for iROI=1:nROI

           clusters_so_far = ROI(iROI).nClusters + clusters_so_far;

           if ~isempty(find(ismember(x_idx,iROI)))

               hold on;

               plot([clusters_so_far clusters_so_far],[1 nClusters],'k-');

               plot([1 nClusters],[clusters_so_far clusters_so_far],'k-');

           end

        end

        ax = gca;
        %ax.XTick = x_idx_jump;
        %ax.XTickLabel = x_label;
        set(ax,'XTick',x_idx_jump);
        set(ax,'XTickLabel',x_label);

        %ax.YTick = x_idx_jump;
        %ax.YTickLabel = x_label;
        set(ax,'YTick',x_idx_jump);
        set(ax,'YTickLabel',x_label);

        xticklabel_rotate;

        print(f,'-depsc',strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast-',strcat(label,'-',sub_label),'.eps'));
        
    end

end

end

function computeFunctionalParcellationForClusters

load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nROI = 90;

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','DMN','AUD'};

nNet = 8;

for iNet=1:length(all_networks_labels)
    
    Net_file = nifti(strcat('LHR-All-Subjects-Functional-Parcels-Pop-Map-',all_networks_labels(iNet),'-seed.nii'));
    Net_file.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually\','LHR-All-Subjects-Functional-Parcels-Pop-Map-',all_networks_labels{iNet},'-seed.nii');
    Net(iNet).img = Net_file.dat(:,:,:);

end

for iROI=1:nROI
    
    nClusters = length(ROI(iROI).clusters);
    
    for iCluster=1:nClusters
    
        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
        
        net_count_label = cell(nNet,2);
        net_count = zeros(nNet,1);
        seed_count = cell.empty;
        
        for iNet=1:nNet
            
            net_count_label(iNet,1) = all_networks_labels(iNet);
            
        end
        
        iSeedIdx = 0;
        seeds_idx = 0;
        
        for iVoxel=1:length(idx_voxels)
            
            [idxx,idxy,idxz] = ind2sub(size(Net(1).img),idx_voxels(iVoxel));
            
            for iNet=1:nNet
              
                if Net(iNet).img(idxx,idxy,idxz) ~= 0
                    
                    iSeedIdx = iSeedIdx + 1;
                    
                    net_count(iNet) = net_count(iNet) + 1;
                    
                    seeds_idx(iSeedIdx) =  squeeze(Net(iNet).img(idxx,idxy,idxz));
                    
                end
                   
            end
            
        end
        
        u_seeds = unique(seeds_idx);
        
        iiSeed = 0;
        
        for iNet=1:nNet
        
            seeds = getFunctionalSeeds_v5(all_networks_labels{iNet});

            nSeeds = length(seeds.ROI);

            for iSeed=1:nSeeds

                if ismember(seeds.ROI(iSeed).idx,u_seeds)

                    iiSeed = iiSeed + 1;
                    seed_count{iiSeed,1} = seeds.ROI(iSeed).label;

                end

            end
            
        end
        
        for iNet=1:nNet
            
            net_count_label{iNet,2} = net_count(iNet);

        end
        
        ROI(iROI).clusters(iCluster).net_count = net_count_label;
        
        ROI(iROI).clusters(iCluster).seed_count = seed_count;
            
    end
    
end

save('FC-Voxels-AAL-ROI-corr-KMeans-Info-Functional-Parcels.mat','ROI');

end

function plotFCAndGrangerAtGroupLevelContrastJB

load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-Contrast.mat');
load(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast','.mat'));
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

do = 'FC';
% do = 'GC';

if strcmp(do,'FC'); attention_z = contrast_attention_z; stimulus_z = contrast_stimulus_z; Analysis_Label = 'FC'; end
if strcmp(do,'GC'); attention_z = Attention_Contrast.Z; stimulus_z = Stimulus_Contrast.Z; Analysis_Label = 'Granger'; end

nRun = 32;
nTotalClusters = 758;
nClusters = 758;
nTR = 150;
pcriterion = 0.01;
nROI = 90;

fs = 5;

% idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
% idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
% idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
% idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
% idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 69 70];
idx_rec_ins_cin = [27 28 29 30 31 32 33 34 35 36];
idx_hc_amyg = [37 38 39 40 41 42];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [71 72 73 74 75 76 77 78];

% x_idx(1) = 29;
% x_idx(2) = 43;
% x_idx(3) = 55;
% x_idx(4) = 56;
% x_idx(5) = 68;
% x_idx(6) = 70;
% x_idx(7) = 78;
% x_idx(8) = 90;
% 
% x_label{1} = 'frontal';
% x_label{2} = 'subcortical';
% x_label{3} = 'occipital';
% x_label{4} = 'temporal';
% x_label{5} = 'parietal';
% x_label{6} = 'frontal';
% x_label{7} = 'subcortical';
% x_label{8} = 'temporal';

cnt = 1;
xtic(cnt) = 2;
xticlb{cnt} = 'precent';
cnt=cnt+1;
xtic(cnt) = 6;
xticlb{cnt} = 'frontsup';
cnt=cnt+1;
xtic(cnt) = 10;
xticlb{cnt} = 'frontmid';
cnt=cnt+1;
xtic(cnt) = 16;
xticlb{cnt} = 'frontinf';
cnt=cnt+1;
xtic(cnt) = 18;
xticlb{5} = 'roland';
cnt=cnt+1;
xtic(cnt) = 20;
xticlb{cnt} = 'suppmot';
cnt=cnt+1;
xtic(cnt) = 22;
xticlb{cnt} = 'olfac';
cnt=cnt+1;
xtic(cnt) = 24;
xticlb{cnt} = 'frontsup';
cnt=cnt+1;
xtic(cnt) = 26;
xticlb{cnt} = 'frontmed';
cnt=cnt+1;
xtic(cnt) = 28;
xticlb{cnt} = 'rectus';
cnt=cnt+1;
xtic(cnt) = 30;
xticlb{cnt} = 'insula';
cnt=cnt+1;
xtic(cnt) = 36;
xticlb{cnt} = 'cingul';
cnt=cnt+1;
xtic(cnt) = 40;
xticlb{13} = 'HC';
cnt=cnt+1;
xtic(cnt) = 42;
xticlb{cnt} = 'amygd';
cnt=cnt+1;
xtic(cnt) = 44;
xticlb{cnt} = 'calcarine';
cnt=cnt+1;
xtic(cnt) = 46;
xticlb{cnt} = 'cuneus';
cnt=cnt+1;
xtic(cnt) = 48;
xticlb{cnt} = 'lingual';
cnt=cnt+1;
xtic(cnt) = 54;
xticlb{cnt} = 'occipit';
cnt=cnt+1;
xtic(cnt) = 56;
xticlb{cnt} = 'fusi';
cnt=cnt+1;
xtic(cnt) = 58;
xticlb{cnt} = 'postcent';
cnt=cnt+1;
xtic(cnt) = 62;
xticlb{cnt} = 'parietal';
cnt=cnt+1;
xtic(cnt) = 64;
xticlb{cnt} = 'suprmarg';
cnt=cnt+1;
xtic(cnt) = 66;
xticlb{cnt} = 'angul';
cnt=cnt+1;
xtic(cnt) = 68;
xticlb{cnt} = 'precun';
cnt=cnt+1;
xtic(cnt) = 70;
xticlb{cnt} = 'paracent';
cnt=cnt+1;
xtic(cnt) = 72;
xticlb{cnt} = 'caudate';
cnt=cnt+1;
xtic(cnt) = 74;
xticlb{cnt} = 'putamen';
cnt=cnt+1;
xtic(cnt) = 76;
xticlb{cnt} = 'pallidum';
cnt=cnt+1;
xtic(cnt) = 78;
xticlb{cnt} = 'thalamus';
cnt=cnt+1;
xtic(cnt) = 84;
xticlb{cnt} = 'tempsup';
cnt=cnt+1;
xtic(cnt) = 88;
xticlb{cnt} = 'tempmid';
cnt=cnt+1;
xtic(cnt) = 90;
xticlb{cnt} = 'tempinf';

nROI = 90;
for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

% zcriterion = 1.6;
zcriterion = 2.3;

for iCondition=1:9
    
    f = figure;
    
    Attention_rho = attention_z;
    Attention_rho(find(Attention_rho>(-1)*zcriterion & Attention_rho<zcriterion)) = 0; 
        
    Stimulus_rho = stimulus_z;
    Stimulus_rho(find(Stimulus_rho>(-1)*zcriterion & Stimulus_rho<zcriterion)) = 0; 
    
    Attention_rho(isnan(Attention_rho)) = 0;
    Stimulus_rho(isnan(Stimulus_rho)) = 0;
        
    if iCondition == 1; rho = attention_z; label = 'Attention'; end
    if iCondition == 2; rho = stimulus_z; label = 'Stimulus'; end
    
    if iCondition == 3; rho = attention_z; label = 'Attention-Only'; end
    if iCondition == 4; rho = stimulus_z; label = 'Stimulus-Only'; end
    
    if iCondition == 5; label = 'Attention-Increase-Stimulus-Increase'; end
    if iCondition == 6; label = 'Attention-Decrease-Stimulus-Decrease'; end
    
    if iCondition == 7; label = 'Attention-Increase-Stimulus-Decrease'; end
    if iCondition == 8; label = 'Attention-Decrease-Stimulus-Increase'; end
    
    if iCondition == 9; label = 'All-Combinations-Together'; end
     
    if iCondition == 3

        rho(find(Stimulus_rho)) = 0; 
    
    end
    
    if iCondition == 4 
        
        rho(find(Attention_rho)) = 0; 
    
    end
%     
%     if iCondition == 5
%         
%         rho = zeros(size(Attention_Contrast.Z));
%         
%         for i=1:size(rho,1)
%             
%             for j=1:size(rho,2)
%                 
%                 if ( Attention_rho(i,j) > 0 ) & ( Stimulus_rho(i,j) > 0 )
%                 
%                     rho(i,j) = max(Attention_rho(i,j),Stimulus_rho(i,j));
%                 
%                 end
%                 
%             end
%             
%         end
%         
%     end
%     
%     if iCondition == 6
%         
%         rho = zeros(size(Attention_Contrast.Z));
%         
%         for i=1:size(rho,1)
%             
%             for j=1:size(rho,2)
%                 
%                 if ( Attention_rho(i,j) < 0 ) & ( Stimulus_rho(i,j) < 0 )
%                     
%                     rho(i,j) = min(Attention_rho(i,j),Stimulus_rho(i,j));
%                 
%                 end
%                 
%             end
%             
%         end
%         
%     end

    if iCondition == 5; rho = Attention_Contrast.Z; rho( ~( ( Attention_rho > 0 ) & ( Stimulus_rho > 0 ) ) ) = 0; end
    
    if iCondition == 6; rho = Attention_Contrast.Z; rho( ~( ( Attention_rho < 0 ) & ( Stimulus_rho < 0 ) ) ) = 0; end
    
    if iCondition == 7; rho = Attention_Contrast.Z; rho( ~( ( Attention_rho > 0 ) & ( Stimulus_rho < 0 ) ) ) = 0; end
    
    if iCondition == 8; rho = Attention_Contrast.Z; rho( ~( ( Attention_rho < 0 ) & ( Stimulus_rho > 0 ) ) ) = 0; end
    
    rho(find(rho>(-1)*zcriterion & rho<zcriterion)) = 0;
    
    if iCondition == 9
        
       rho = zeros(size(Attention_rho));
        
       %%% ALL COMBINATIONS
%        rho(Attention_rho & ~Stimulus_rho) = 3;
%        rho(~Attention_rho & Stimulus_rho) = -3;
%        
%        rho(( Attention_rho > 0 ) & ( Stimulus_rho > 0 )) = 2;
%        rho(( Attention_rho < 0 ) & ( Stimulus_rho < 0 )) = 0.7;
%        
%        rho(( Attention_rho > 0 ) & ( Stimulus_rho < 0 )) = -1;
%        rho(( Attention_rho < 0 ) & ( Stimulus_rho > 0 )) = -0.3;
       
       %%% ONLY ON BOTH CONTRASTS
       rho(Attention_rho & ~Stimulus_rho) = 0;
       rho(~Attention_rho & Stimulus_rho) = 0;
       
       rho(( Attention_rho > 0 ) & ( Stimulus_rho > 0 )) = 3;
       rho(( Attention_rho < 0 ) & ( Stimulus_rho < 0 )) = 3;
       
       rho(( Attention_rho > 0 ) & ( Stimulus_rho < 0 )) = 3;
       rho(( Attention_rho < 0 ) & ( Stimulus_rho > 0 )) = 3;
        
    end
    
    rho(isnan(rho)) = 0;
     
    imagesc(rho);
    hold on
   
    title(label);

    clrmp = colormap('jet');
    clrmp(65,:) = clrmp(64,:);
    clrmp(33,:) = [1 1 1];

    min_C = -3;
    max_C = 3;

    %colorbar;
    caxis([min_C max_C]);
    %colorbar('Ticks',[min_C max_C/2 max_C]);
    colormap(clrmp);
    colorbar;
    
   for i = 1 : length(xtic)      % draw anatomical boundaries
    
        jump = ROI_info{ xtic(i), 5};

        plot(0.5+jump*[1 1],0.5+[0 nClusters],'k-');
        plot([0.5+0 nClusters],0.5+jump*[1 1],'k-');
    
        if i>1
            xtic_jump(i) = 0.5 * ( ROI_info{ xtic(i-1), 5} + ROI_info{ xtic(i), 5} );
        else
            xtic_jump(i) = 0.5 * ( 1 + ROI_info{ xtic(i), 5} );
        end
    end

    hold off;

    ax = gca;
    axis 'square';
    axis([ 0.5 nClusters+0.5 0.5 nClusters+0.5] );
    %set(ax,'XTick', [] );
    set(ax,'XTick',xtic_jump);
    set(ax,'XTickLabel',xticlb);
    %rotateXLabels( ax, 45 );

    set(ax,'YTick',xtic_jump);
    set(ax,'YTickLabel',xticlb);
    %rotateYLabels( ax, 45 );
    
     set(ax,'FontSize',fs);

    xticklabel_rotate([],90,[],'Fontsize',fs);

    print(f,'-depsc',strcat('FC_Voxel_AAL_ROI_kmeans_',Analysis_Label,'_Clusters','-','Mean-Contrast-',label,'.eps'));

end

end

function plotGrangerAtGroupLevelContrastJBSubcortical

load(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast','.mat'));
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nRun = 32;
nTotalClusters = 758;
nClusters = 758;
nTR = 150;
pcriterion = 0.01;
nROI = 90;

% idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
% idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
% idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
% idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
% idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 69 70];
idx_rec_ins_cin = [27 28 29 30 31 32 33 34 35 36];
idx_hc_amyg = [37 38 39 40 41 42];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [71 72 73 74 75 76 77 78];

% x_idx(1) = 29;
% x_idx(2) = 43;
% x_idx(3) = 55;
% x_idx(4) = 56;
% x_idx(5) = 68;
% x_idx(6) = 70;
% x_idx(7) = 78;
% x_idx(8) = 90;
% 
% x_label{1} = 'frontal';
% x_label{2} = 'subcortical';
% x_label{3} = 'occipital';
% x_label{4} = 'temporal';
% x_label{5} = 'parietal';
% x_label{6} = 'frontal';
% x_label{7} = 'subcortical';
% x_label{8} = 'temporal';

cnt = 1;
ytic(cnt) = 2;
yticlb{cnt} = 'precent';
cnt=cnt+1;
ytic(cnt) = 6;
yticlb{cnt} = 'frontsup';
cnt=cnt+1;
ytic(cnt) = 10;
yticlb{cnt} = 'frontmid';
cnt=cnt+1;
ytic(cnt) = 16;
yticlb{cnt} = 'frontinf';
cnt=cnt+1;
ytic(cnt) = 18;
yticlb{5} = 'roland';
cnt=cnt+1;
ytic(cnt) = 20;
yticlb{cnt} = 'suppmot';
cnt=cnt+1;
ytic(cnt) = 22;
yticlb{cnt} = 'olfac';
cnt=cnt+1;
ytic(cnt) = 24;
yticlb{cnt} = 'frontsup';
cnt=cnt+1;
ytic(cnt) = 26;
yticlb{cnt} = 'frontmed';
cnt=cnt+1;
ytic(cnt) = 28;
yticlb{cnt} = 'rectus';
cnt=cnt+1;
ytic(cnt) = 30;
yticlb{cnt} = 'insula';
cnt=cnt+1;
ytic(cnt) = 36;
yticlb{cnt} = 'cingul';
cnt=cnt+1;
ytic(cnt) = 40;
yticlb{13} = 'HC';
cnt=cnt+1;
ytic(cnt) = 42;
yticlb{cnt} = 'amygd';
cnt=cnt+1;
ytic(cnt) = 44;
yticlb{cnt} = 'calcarine';
cnt=cnt+1;
ytic(cnt) = 46;
yticlb{cnt} = 'cuneus';
cnt=cnt+1;
ytic(cnt) = 48;
yticlb{cnt} = 'lingual';
cnt=cnt+1;
ytic(cnt) = 54;
yticlb{cnt} = 'occipit';
cnt=cnt+1;
ytic(cnt) = 56;
yticlb{cnt} = 'fusi';
cnt=cnt+1;
ytic(cnt) = 58;
yticlb{cnt} = 'postcent';
cnt=cnt+1;
ytic(cnt) = 62;
yticlb{cnt} = 'parietal';
cnt=cnt+1;
ytic(cnt) = 64;
yticlb{cnt} = 'suprmarg';
cnt=cnt+1;
ytic(cnt) = 66;
yticlb{cnt} = 'angul';
cnt=cnt+1;
ytic(cnt) = 68;
yticlb{cnt} = 'precun';
cnt=cnt+1;
ytic(cnt) = 70;
yticlb{cnt} = 'paracent';
cnt=cnt+1;
ytic(cnt) = 72;
yticlb{cnt} = 'caudate';
cnt=cnt+1;
ytic(cnt) = 74;
yticlb{cnt} = 'putamen';
cnt=cnt+1;
ytic(cnt) = 76;
yticlb{cnt} = 'pallidum';
cnt=cnt+1;
ytic(cnt) = 78;
yticlb{cnt} = 'thalamus';
cnt=cnt+1;
ytic(cnt) = 84;
yticlb{cnt} = 'tempsup';
cnt=cnt+1;
ytic(cnt) = 88;
yticlb{cnt} = 'tempmid';
cnt=cnt+1;
ytic(cnt) = 90;
yticlb{cnt} = 'tempinf';

xticlb = {'Caudate-L';'Caudate-R';'Putamen-L';'Putamen-R';'Pallidum-L';'Pallidum-R';'Thalamus-L';'Thalamus-R'};

nROI = 90;
for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

zcriterion = 1.6;
%zcriterion = 2.3;

for iCondition=1:4
    
    f = figure;
    
    if mod(iCondition,2) ~= 0; rho = Attention_Contrast.Z; label = 'Attention'; end
    if mod(iCondition,2) == 0; rho = Stimulus_Contrast.Z; label = 'Stimulus'; end
    
    rho(find(rho>(-1)*zcriterion & rho<zcriterion)) = 0; 
    
    if iCondition == 1 || iCondition == 2
        
        imagesc(rho(ROI_info{idx_subcortical(1),4}:ROI_info{idx_subcortical(end),5},:));
        hold on
        
    else
       
        imagesc(rho(:,ROI_info{idx_subcortical(1),4}:ROI_info{idx_subcortical(end),5}));
        hold on
        
    end
   
    if iCondition == 1 || iCondition == 2
    
        direction = 'FROM';
        title(strcat(label,'-',direction));
        
    else
       
        direction = 'TO';
        title(strcat(label,'-',direction));
        
    end

    clrmp = colormap('jet');
    clrmp(65,:) = clrmp(64,:);
    clrmp(33,:) = [1 1 1];

    min_C = -3;
    max_C = 3;

    %colorbar;
    caxis([min_C max_C]);
    %colorbar('Ticks',[min_C max_C/2 max_C]);
    colormap(clrmp);
    colorbar;
    
   for i = 1 : length(ytic)      % draw anatomical boundaries
    
        jump = ROI_info{ ytic(i), 5};

        if iCondition == 1 || iCondition == 2
        
            plot([0.5+0 nClusters],0.5+jump*[1 1],'k-');
            
        else
            
            plot(0.5+jump*[1 1],0.5+[0 nClusters],'k-');
            
        end
    
        if i>1
            ytic_jump(i) = 0.5 * ( ROI_info{ ytic(i-1), 5} + ROI_info{ ytic(i), 5} );
        else
            ytic_jump(i) = 0.5 * ( 1 + ROI_info{ ytic(i), 5} );
        end
   end

    for i = 1 : length(idx_subcortical)      % draw anatomical boundaries
    
        if i>1
            xtic_jump(i) = 0.5 * ( ROI_info{ idx_subcortical(i-1), 5} + ROI_info{ idx_subcortical(i), 5} - 2*ROI_info{ idx_subcortical(1)-1, 5} );
            jump = ROI_info{ idx_subcortical(i), 5} - ROI_info{ idx_subcortical(1)-1, 5};
        else
            xtic_jump(i) = 0.5 * ( 1 + ROI_info{ idx_subcortical(i), 5} - ROI_info{ idx_subcortical(1)-1, 5});
            jump = ROI_info{ idx_subcortical(i), 5} - ROI_info{ idx_subcortical(1)-1, 5};
        end
    
        if iCondition == 1 || iCondition == 2
        
            plot([0.5+0 nClusters],0.5+jump*[1 1],'k-');
            
        else
            
            plot(0.5+jump*[1 1],0.5+[0 nClusters],'k-');
            
        end
        
    end
    hold off;

    if iCondition == 1 || iCondition == 2
        
        ax = gca;
        %axis 'square';
        %axis([ 0.5 nClusters+0.5 0.5 nClusters+0.5] );
        %set(ax,'XTick', [] );
        set(ax,'XTick',ytic_jump);
        set(ax,'XTickLabel',yticlb);
        %rotateXLabels( ax, 45 );

        set(ax,'YTick',xtic_jump);
        set(ax,'YTickLabel',xticlb);
        %rotateYLabels( ax, 45 );
        
        set(ax,'FontSize',5);

        xticklabel_rotate([],90,[],'Fontsize',5);
        
    else
       
        ax = gca;
        %axis 'square';
        %axis([ 0.5 nClusters+0.5 0.5 nClusters+0.5] );
        %set(ax,'XTick', [] );
        set(ax,'XTick',xtic_jump);
        set(ax,'XTickLabel',xticlb);
        %rotateXLabels( ax, 45 );

        set(ax,'YTick',ytic_jump);
        set(ax,'YTickLabel',yticlb);
        %rotateYLabels( ax, 45 );
        
        set(ax,'FontSize',5);

        xticklabel_rotate([],90,[],'Fontsize',5);
        
    end

    print(f,'-depsc',strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast-',label,'-subcortical-',direction,'-z-',int2str(zcriterion),'.eps'));

end

end

function plotFCClusterContrastJBSubcortical

load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-Contrast.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nRun = 32;
nTotalClusters = 758;
nClusters = 758;
nTR = 150;
pcriterion = 0.01;
nROI = 90;

% idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
% idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
% idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
% idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
% idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 69 70];
idx_rec_ins_cin = [27 28 29 30 31 32 33 34 35 36];
idx_hc_amyg = [37 38 39 40 41 42];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [71 72 73 74 75 76 77 78];

% x_idx(1) = 29;
% x_idx(2) = 43;
% x_idx(3) = 55;
% x_idx(4) = 56;
% x_idx(5) = 68;
% x_idx(6) = 70;
% x_idx(7) = 78;
% x_idx(8) = 90;
% 
% x_label{1} = 'frontal';
% x_label{2} = 'subcortical';
% x_label{3} = 'occipital';
% x_label{4} = 'temporal';
% x_label{5} = 'parietal';
% x_label{6} = 'frontal';
% x_label{7} = 'subcortical';
% x_label{8} = 'temporal';

cnt = 1;
ytic(cnt) = 2;
yticlb{cnt} = 'precent';
cnt=cnt+1;
ytic(cnt) = 6;
yticlb{cnt} = 'frontsup';
cnt=cnt+1;
ytic(cnt) = 10;
yticlb{cnt} = 'frontmid';
cnt=cnt+1;
ytic(cnt) = 16;
yticlb{cnt} = 'frontinf';
cnt=cnt+1;
ytic(cnt) = 18;
yticlb{5} = 'roland';
cnt=cnt+1;
ytic(cnt) = 20;
yticlb{cnt} = 'suppmot';
cnt=cnt+1;
ytic(cnt) = 22;
yticlb{cnt} = 'olfac';
cnt=cnt+1;
ytic(cnt) = 24;
yticlb{cnt} = 'frontsup';
cnt=cnt+1;
ytic(cnt) = 26;
yticlb{cnt} = 'frontmed';
cnt=cnt+1;
ytic(cnt) = 28;
yticlb{cnt} = 'rectus';
cnt=cnt+1;
ytic(cnt) = 30;
yticlb{cnt} = 'insula';
cnt=cnt+1;
ytic(cnt) = 36;
yticlb{cnt} = 'cingul';
cnt=cnt+1;
ytic(cnt) = 40;
yticlb{13} = 'HC';
cnt=cnt+1;
ytic(cnt) = 42;
yticlb{cnt} = 'amygd';
cnt=cnt+1;
ytic(cnt) = 44;
yticlb{cnt} = 'calcarine';
cnt=cnt+1;
ytic(cnt) = 46;
yticlb{cnt} = 'cuneus';
cnt=cnt+1;
ytic(cnt) = 48;
yticlb{cnt} = 'lingual';
cnt=cnt+1;
ytic(cnt) = 54;
yticlb{cnt} = 'occipit';
cnt=cnt+1;
ytic(cnt) = 56;
yticlb{cnt} = 'fusi';
cnt=cnt+1;
ytic(cnt) = 58;
yticlb{cnt} = 'postcent';
cnt=cnt+1;
ytic(cnt) = 62;
yticlb{cnt} = 'parietal';
cnt=cnt+1;
ytic(cnt) = 64;
yticlb{cnt} = 'suprmarg';
cnt=cnt+1;
ytic(cnt) = 66;
yticlb{cnt} = 'angul';
cnt=cnt+1;
ytic(cnt) = 68;
yticlb{cnt} = 'precun';
cnt=cnt+1;
ytic(cnt) = 70;
yticlb{cnt} = 'paracent';
cnt=cnt+1;
ytic(cnt) = 72;
yticlb{cnt} = 'caudate';
cnt=cnt+1;
ytic(cnt) = 74;
yticlb{cnt} = 'putamen';
cnt=cnt+1;
ytic(cnt) = 76;
yticlb{cnt} = 'pallidum';
cnt=cnt+1;
ytic(cnt) = 78;
yticlb{cnt} = 'thalamus';
cnt=cnt+1;
ytic(cnt) = 84;
yticlb{cnt} = 'tempsup';
cnt=cnt+1;
ytic(cnt) = 88;
yticlb{cnt} = 'tempmid';
cnt=cnt+1;
ytic(cnt) = 90;
yticlb{cnt} = 'tempinf';

xticlb = {'Caudate-L';'Caudate-R';'Putamen-L';'Putamen-R';'Pallidum-L';'Pallidum-R';'Thalamus-L';'Thalamus-R'};

nROI = 90;
for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

zcriterion = 1.6;
%zcriterion = 2.3;

for iCondition=1:4
    
    f = figure;
    
    if mod(iCondition,2) ~= 0; rho = contrast_attention_z; label = 'Attention'; end
    if mod(iCondition,2) == 0; rho = contrast_stimulus_z; label = 'Stimulus'; end
    
    rho(find(rho>(-1)*zcriterion & rho<zcriterion)) = 0; 
    
    if iCondition == 1 || iCondition == 2
        
        imagesc(rho(ROI_info{idx_subcortical(1),4}:ROI_info{idx_subcortical(end),5},:));
        hold on
        
    else
       
        imagesc(rho(:,ROI_info{idx_subcortical(1),4}:ROI_info{idx_subcortical(end),5}));
        hold on
        
    end
   
    if iCondition == 1 || iCondition == 2
    
        direction = 'FROM';
        title(strcat(label,'-',direction));
        
    else
       
        direction = 'TO';
        title(strcat(label,'-',direction));
        
    end

    clrmp = colormap('jet');
    clrmp(65,:) = clrmp(64,:);
    clrmp(33,:) = [1 1 1];

    min_C = -3;
    max_C = 3;

    %colorbar;
    caxis([min_C max_C]);
    %colorbar('Ticks',[min_C max_C/2 max_C]);
    colormap(clrmp);
    colorbar;
    
   for i = 1 : length(ytic)      % draw anatomical boundaries
    
        jump = ROI_info{ ytic(i), 5};

        if iCondition == 1 || iCondition == 2
        
            plot([0.5+0 nClusters],0.5+jump*[1 1],'k-');
            
        else
            
            plot(0.5+jump*[1 1],0.5+[0 nClusters],'k-');
            
        end
    
        if i>1
            ytic_jump(i) = 0.5 * ( ROI_info{ ytic(i-1), 5} + ROI_info{ ytic(i), 5} );
        else
            ytic_jump(i) = 0.5 * ( 1 + ROI_info{ ytic(i), 5} );
        end
   end

    for i = 1 : length(idx_subcortical)      % draw anatomical boundaries
    
        if i>1
            xtic_jump(i) = 0.5 * ( ROI_info{ idx_subcortical(i-1), 5} + ROI_info{ idx_subcortical(i), 5} - 2*ROI_info{ idx_subcortical(1)-1, 5} );
            jump = ROI_info{ idx_subcortical(i), 5} - ROI_info{ idx_subcortical(1)-1, 5};
        else
            xtic_jump(i) = 0.5 * ( 1 + ROI_info{ idx_subcortical(i), 5} - ROI_info{ idx_subcortical(1)-1, 5});
            jump = ROI_info{ idx_subcortical(i), 5} - ROI_info{ idx_subcortical(1)-1, 5};
        end
    
        if iCondition == 1 || iCondition == 2
        
            plot([0.5+0 nClusters],0.5+jump*[1 1],'k-');
            
        else
            
            plot(0.5+jump*[1 1],0.5+[0 nClusters],'k-');
            
        end
        
    end
    hold off;

    if iCondition == 1 || iCondition == 2
        
        ax = gca;
        %axis 'square';
        %axis([ 0.5 nClusters+0.5 0.5 nClusters+0.5] );
        %set(ax,'XTick', [] );
        set(ax,'XTick',ytic_jump);
        set(ax,'XTickLabel',yticlb);
        %rotateXLabels( ax, 45 );

        set(ax,'YTick',xtic_jump);
        set(ax,'YTickLabel',xticlb);
        %rotateYLabels( ax, 45 );
        
        set(ax,'FontSize',5);

        xticklabel_rotate([],90,[],'Fontsize',5);
        
    else
       
        ax = gca;
        %axis 'square';
        %axis([ 0.5 nClusters+0.5 0.5 nClusters+0.5] );
        %set(ax,'XTick', [] );
        set(ax,'XTick',xtic_jump);
        set(ax,'XTickLabel',xticlb);
        %rotateXLabels( ax, 45 );

        set(ax,'YTick',ytic_jump);
        set(ax,'YTickLabel',yticlb);
        %rotateYLabels( ax, 45 );
        
        set(ax,'FontSize',5);

        xticklabel_rotate([],90,[],'Fontsize',5);
        
    end

    print(f,'-depsc',strcat('FC_Voxel_AAL_ROI_kmeans_Contrast_Clusters','-','Mean-Contrast-',label,'-subcortical-',direction,'.eps'));

end

end

function plotGrangerAtGroupLevelContrastJBAllLobes

load(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast','.mat'));
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRun = 32;
nTotalClusters = 758;
nClusters = 758;
nTR = 150;
pcriterion = 0.01;
nROI = 90;

% idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
% idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
% idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
% idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
% idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

group(1).idx = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 69 70]; % idx_frontal
group(2).idx = [27 28 29 30 31 32 33 34 35 36]; % idx_rec_ins_cin
group(3).idx = [37 38 39 40 41 42]; % idx_hc_amyg
group(4).idx = [43 44 45 46 47 48 49 50 51 52 53 54]; % idx_occipital
group(5).idx = [57 58 59 60 61 62 63 64 65 66 67 68]; % idx_parietal
group(6).idx = [55 56 79 80 81 82 83 84 85 86 87 88 89 90]; % idx_temporal
group(7).idx = [71 72 73 74 75 76 77 78]; % idx_subcortical

group_label{1} = 'Frontal';
group_label{2} = 'Rec_Ins_Cin';
group_label{3} = 'HC_Amyg';
group_label{4} = 'Occipital';
group_label{5} = 'Parietal';
group_label{6} = 'Temporal';
group_label{7} = 'Subcortical';

% x_idx(1) = 29;
% x_idx(2) = 43;
% x_idx(3) = 55;
% x_idx(4) = 56;
% x_idx(5) = 68;
% x_idx(6) = 70;
% x_idx(7) = 78;
% x_idx(8) = 90;
% 
% x_label{1} = 'frontal';
% x_label{2} = 'subcortical';
% x_label{3} = 'occipital';
% x_label{4} = 'temporal';
% x_label{5} = 'parietal';
% x_label{6} = 'frontal';
% x_label{7} = 'subcortical';
% x_label{8} = 'temporal';

cnt = 1;
ytic(cnt) = 2;
yticlb{cnt} = 'precent';
cnt=cnt+1;
ytic(cnt) = 6;
yticlb{cnt} = 'frontsup';
cnt=cnt+1;
ytic(cnt) = 10;
yticlb{cnt} = 'frontmid';
cnt=cnt+1;
ytic(cnt) = 16;
yticlb{cnt} = 'frontinf';
cnt=cnt+1;
ytic(cnt) = 18;
yticlb{5} = 'roland';
cnt=cnt+1;
ytic(cnt) = 20;
yticlb{cnt} = 'suppmot';
cnt=cnt+1;
ytic(cnt) = 22;
yticlb{cnt} = 'olfac';
cnt=cnt+1;
ytic(cnt) = 24;
yticlb{cnt} = 'frontsup';
cnt=cnt+1;
ytic(cnt) = 26;
yticlb{cnt} = 'frontmed';
cnt=cnt+1;
ytic(cnt) = 28;
yticlb{cnt} = 'rectus';
cnt=cnt+1;
ytic(cnt) = 30;
yticlb{cnt} = 'insula';
cnt=cnt+1;
ytic(cnt) = 36;
yticlb{cnt} = 'cingul';
cnt=cnt+1;
ytic(cnt) = 40;
yticlb{13} = 'HC';
cnt=cnt+1;
ytic(cnt) = 42;
yticlb{cnt} = 'amygd';
cnt=cnt+1;
ytic(cnt) = 44;
yticlb{cnt} = 'calcarine';
cnt=cnt+1;
ytic(cnt) = 46;
yticlb{cnt} = 'cuneus';
cnt=cnt+1;
ytic(cnt) = 48;
yticlb{cnt} = 'lingual';
cnt=cnt+1;
ytic(cnt) = 54;
yticlb{cnt} = 'occipit';
cnt=cnt+1;
ytic(cnt) = 56;
yticlb{cnt} = 'fusi';
cnt=cnt+1;
ytic(cnt) = 58;
yticlb{cnt} = 'postcent';
cnt=cnt+1;
ytic(cnt) = 62;
yticlb{cnt} = 'parietal';
cnt=cnt+1;
ytic(cnt) = 64;
yticlb{cnt} = 'suprmarg';
cnt=cnt+1;
ytic(cnt) = 66;
yticlb{cnt} = 'angul';
cnt=cnt+1;
ytic(cnt) = 68;
yticlb{cnt} = 'precun';
cnt=cnt+1;
ytic(cnt) = 70;
yticlb{cnt} = 'paracent';
cnt=cnt+1;
ytic(cnt) = 72;
yticlb{cnt} = 'caudate';
cnt=cnt+1;
ytic(cnt) = 74;
yticlb{cnt} = 'putamen';
cnt=cnt+1;
ytic(cnt) = 76;
yticlb{cnt} = 'pallidum';
cnt=cnt+1;
ytic(cnt) = 78;
yticlb{cnt} = 'thalamus';
cnt=cnt+1;
ytic(cnt) = 84;
yticlb{cnt} = 'tempsup';
cnt=cnt+1;
ytic(cnt) = 88;
yticlb{cnt} = 'tempmid';
cnt=cnt+1;
ytic(cnt) = 90;
yticlb{cnt} = 'tempinf';

% xticlb = {'Caudate-L';'Caudate-R';'Putamen-L';'Putamen-R';'Pallidum-L';'Pallidum-R';'Thalamus-L';'Thalamus-R'};

nROI = 90;
for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

%zcriterion = 1.6;
zcriterion = 2.3;

nGroup = 7;

for iGroup=1:nGroup
    
    indexes = [];
    
    for idx=1:length(group(iGroup).idx)
        
        indexes = [indexes, ROI_info{group(iGroup).idx(idx),4}:ROI_info{group(iGroup).idx(idx),5}];
        
        if idx > 1
            
            ROI_info_v2{idx,4} = ROI_info_v2{idx-1,4} + 1;
            ROI_info_v2{idx,5} = ROI_info_v2{idx-1,5} + length(ROI_info{group(iGroup).idx(idx),4}:ROI_info{group(iGroup).idx(idx),5});
        
        else
        
            ROI_info_v2{idx,4} = 1;
            ROI_info_v2{idx,5} = length(ROI_info{group(iGroup).idx(idx),4}:ROI_info{group(iGroup).idx(idx),5});
            
        end
        
    end

    for idx=1:length(group(iGroup).idx)
       
        xticlb{idx} = AAL_ROI(group(iGroup).idx(idx)).Nom_L;
        
    end
    
    for iCondition=1:4

        f = figure;

        if mod(iCondition,2) ~= 0; rho = Attention_Contrast.Z; label = [group_label{iGroup} ' - ' 'Attention']; end
        if mod(iCondition,2) == 0; rho = Stimulus_Contrast.Z; label = [group_label{iGroup} ' - ' 'Stimulus']; end

        rho(find(rho>(-1)*zcriterion & rho<zcriterion)) = 0; 

        if iCondition == 1 || iCondition == 2

            imagesc(rho(indexes,:));
            hold on

        else

            imagesc(rho(:,indexes));
            hold on

        end

        if iCondition == 1 || iCondition == 2

            direction = 'FROM';
            title(strcat(label,'-',direction));

        else

            direction = 'TO';
            title(strcat(label,'-',direction));

        end

        clrmp = colormap('jet');
        clrmp(65,:) = clrmp(64,:);
        clrmp(33,:) = [1 1 1];

        min_C = -3;
        max_C = 3;

        %colorbar;
        caxis([min_C max_C]);
        %colorbar('Ticks',[min_C max_C/2 max_C]);
        colormap(clrmp);
        colorbar;

       for i = 1 : length(ytic)      % draw anatomical boundaries

            jump = ROI_info{ ytic(i), 5};

            if iCondition == 1 || iCondition == 2

                plot([0.5+0 nClusters],0.5+jump*[1 1],'k-');

            else

                plot(0.5+jump*[1 1],0.5+[0 nClusters],'k-');

            end

            if i>1
                ytic_jump(i) = 0.5 * ( ROI_info{ ytic(i-1), 5} + ROI_info{ ytic(i), 5} );
            else
                ytic_jump(i) = 0.5 * ( 1 + ROI_info{ ytic(i), 5} );
            end
       end

        for i = 1 : length(group(iGroup).idx)      % draw anatomical boundaries

            if i>1
                xtic_jump(i) = 0.5 * ( ROI_info_v2{ i-1, 5} + ROI_info_v2{ i, 5} );
                jump = ROI_info_v2{ i, 5 };
            else
                
                if iGroup > 1
                    
                    xtic_jump(i) = 0.5 * ( 1 + ROI_info_v2{ i, 5});
                    jump = ROI_info_v2{ i, 5};
               
                else
                
                    xtic_jump(i) = 0.5 * ( 1 + ROI_info_v2{ i, 5});
                    jump = ROI_info_v2{ i, 5 };
                    
                end
                
            end

            if iCondition == 1 || iCondition == 2

                plot([0.5+0 nClusters],0.5+jump*[1 1],'k-');

            else

                plot(0.5+jump*[1 1],0.5+[0 nClusters],'k-');

            end

        end
        hold off;

        if iCondition == 1 || iCondition == 2

            ax = gca;
            %axis 'square';
            %axis([ 0.5 nClusters+0.5 0.5 nClusters+0.5] );
            %set(ax,'XTick', [] );
            set(ax,'XTick',ytic_jump);
            set(ax,'XTickLabel',yticlb);
            %rotateXLabels( ax, 45 );

            set(ax,'YTick',xtic_jump);
            set(ax,'YTickLabel',xticlb);
            %rotateYLabels( ax, 45 );

            set(ax,'FontSize',5);

            xticklabel_rotate([],90,[],'Fontsize',5);

        else

            ax = gca;
            %axis 'square';
            %axis([ 0.5 nClusters+0.5 0.5 nClusters+0.5] );
            %set(ax,'XTick', [] );
            set(ax,'XTick',xtic_jump);
            set(ax,'XTickLabel',xticlb);
            %rotateXLabels( ax, 45 );

            set(ax,'YTick',ytic_jump);
            set(ax,'YTickLabel',yticlb);
            %rotateYLabels( ax, 45 );

            set(ax,'FontSize',5);

            xticklabel_rotate([],90,[],'Fontsize',5);

        end

        print(f,'-depsc',strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast-',label,'-',group_label{iGroup},'-',direction,'-z-',num2str(zcriterion),'.eps'));

    end

end

end

function plotFCClusterContrastJBAllLobes

load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-Contrast.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRun = 32;
nTotalClusters = 758;
nClusters = 758;
nTR = 150;
pcriterion = 0.01;
nROI = 90;

% idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
% idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
% idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
% idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
% idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 69 70];
idx_rec_ins_cin = [27 28 29 30 31 32 33 34 35 36];
idx_hc_amyg = [37 38 39 40 41 42];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [71 72 73 74 75 76 77 78];

group(1).idx = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 69 70]; % idx_frontal
group(2).idx = [27 28 29 30 31 32 33 34 35 36]; % idx_rec_ins_cin
group(3).idx = [37 38 39 40 41 42]; % idx_hc_amyg
group(4).idx = [43 44 45 46 47 48 49 50 51 52 53 54]; % idx_occipital
group(5).idx = [57 58 59 60 61 62 63 64 65 66 67 68]; % idx_parietal
group(6).idx = [55 56 79 80 81 82 83 84 85 86 87 88 89 90]; % idx_temporal
group(7).idx = [71 72 73 74 75 76 77 78]; % idx_subcortical

group_label{1} = 'Frontal';
group_label{2} = 'Rec_Ins_Cin';
group_label{3} = 'HC_Amyg';
group_label{4} = 'Occipital';
group_label{5} = 'Parietal';
group_label{6} = 'Temporal';
group_label{7} = 'Subcortical';

% x_idx(1) = 29;
% x_idx(2) = 43;
% x_idx(3) = 55;
% x_idx(4) = 56;
% x_idx(5) = 68;
% x_idx(6) = 70;
% x_idx(7) = 78;
% x_idx(8) = 90;
% 
% x_label{1} = 'frontal';
% x_label{2} = 'subcortical';
% x_label{3} = 'occipital';
% x_label{4} = 'temporal';
% x_label{5} = 'parietal';
% x_label{6} = 'frontal';
% x_label{7} = 'subcortical';
% x_label{8} = 'temporal';

cnt = 1;
ytic(cnt) = 2;
yticlb{cnt} = 'precent';
cnt=cnt+1;
ytic(cnt) = 6;
yticlb{cnt} = 'frontsup';
cnt=cnt+1;
ytic(cnt) = 10;
yticlb{cnt} = 'frontmid';
cnt=cnt+1;
ytic(cnt) = 16;
yticlb{cnt} = 'frontinf';
cnt=cnt+1;
ytic(cnt) = 18;
yticlb{5} = 'roland';
cnt=cnt+1;
ytic(cnt) = 20;
yticlb{cnt} = 'suppmot';
cnt=cnt+1;
ytic(cnt) = 22;
yticlb{cnt} = 'olfac';
cnt=cnt+1;
ytic(cnt) = 24;
yticlb{cnt} = 'frontsup';
cnt=cnt+1;
ytic(cnt) = 26;
yticlb{cnt} = 'frontmed';
cnt=cnt+1;
ytic(cnt) = 28;
yticlb{cnt} = 'rectus';
cnt=cnt+1;
ytic(cnt) = 30;
yticlb{cnt} = 'insula';
cnt=cnt+1;
ytic(cnt) = 36;
yticlb{cnt} = 'cingul';
cnt=cnt+1;
ytic(cnt) = 40;
yticlb{13} = 'HC';
cnt=cnt+1;
ytic(cnt) = 42;
yticlb{cnt} = 'amygd';
cnt=cnt+1;
ytic(cnt) = 44;
yticlb{cnt} = 'calcarine';
cnt=cnt+1;
ytic(cnt) = 46;
yticlb{cnt} = 'cuneus';
cnt=cnt+1;
ytic(cnt) = 48;
yticlb{cnt} = 'lingual';
cnt=cnt+1;
ytic(cnt) = 54;
yticlb{cnt} = 'occipit';
cnt=cnt+1;
ytic(cnt) = 56;
yticlb{cnt} = 'fusi';
cnt=cnt+1;
ytic(cnt) = 58;
yticlb{cnt} = 'postcent';
cnt=cnt+1;
ytic(cnt) = 62;
yticlb{cnt} = 'parietal';
cnt=cnt+1;
ytic(cnt) = 64;
yticlb{cnt} = 'suprmarg';
cnt=cnt+1;
ytic(cnt) = 66;
yticlb{cnt} = 'angul';
cnt=cnt+1;
ytic(cnt) = 68;
yticlb{cnt} = 'precun';
cnt=cnt+1;
ytic(cnt) = 70;
yticlb{cnt} = 'paracent';
cnt=cnt+1;
ytic(cnt) = 72;
yticlb{cnt} = 'caudate';
cnt=cnt+1;
ytic(cnt) = 74;
yticlb{cnt} = 'putamen';
cnt=cnt+1;
ytic(cnt) = 76;
yticlb{cnt} = 'pallidum';
cnt=cnt+1;
ytic(cnt) = 78;
yticlb{cnt} = 'thalamus';
cnt=cnt+1;
ytic(cnt) = 84;
yticlb{cnt} = 'tempsup';
cnt=cnt+1;
ytic(cnt) = 88;
yticlb{cnt} = 'tempmid';
cnt=cnt+1;
ytic(cnt) = 90;
yticlb{cnt} = 'tempinf';

% xticlb = {'Caudate-L';'Caudate-R';'Putamen-L';'Putamen-R';'Pallidum-L';'Pallidum-R';'Thalamus-L';'Thalamus-R'};

nROI = 90;
for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

%zcriterion = 1.6;
zcriterion = 2.3;

nGroup = 7;

for iGroup=1:nGroup
    
    indexes = [];
    
    for idx=1:length(group(iGroup).idx)
        
        indexes = [indexes, ROI_info{group(iGroup).idx(idx),4}:ROI_info{group(iGroup).idx(idx),5}];
        
    end
    
    
    if idx > 1

        ROI_info_v2{idx,4} = ROI_info_v2{idx-1,4} + 1;
        ROI_info_v2{idx,5} = ROI_info_v2{idx-1,5} + length(ROI_info{group(iGroup).idx(idx),4}:ROI_info{group(iGroup).idx(idx),5});

    else

        ROI_info_v2{idx,4} = 1;
        ROI_info_v2{idx,5} = length(ROI_info{group(iGroup).idx(idx),4}:ROI_info{group(iGroup).idx(idx),5});

    end
    
    
    for idx=1:length(group(iGroup).idx)
       
        xticlb{idx} = AAL_ROI(group(iGroup).idx(idx)).Nom_L;
        
    end

    for iCondition=1:4

        f = figure;

        if mod(iCondition,2) ~= 0; rho = contrast_attention_z; label = [group_label{iGroup} ' - ' 'Attention']; end
        if mod(iCondition,2) == 0; rho = contrast_stimulus_z; label = [group_label{iGroup} ' - ' 'Stimulus']; end

        rho(find(rho>(-1)*zcriterion & rho<zcriterion)) = 0; 

        if iCondition == 1 || iCondition == 2

            imagesc(rho(indexes,:));
            hold on

        else

            imagesc(rho(:,indexes));
            hold on

        end

        if iCondition == 1 || iCondition == 2

            direction = 'FROM';
            title(strcat(label,'-',direction));

        else

            direction = 'TO';
            title(strcat(label,'-',direction));

        end

        clrmp = colormap('jet');
        clrmp(65,:) = clrmp(64,:);
        clrmp(33,:) = [1 1 1];

        min_C = -3;
        max_C = 3;

        %colorbar;
        caxis([min_C max_C]);
        %colorbar('Ticks',[min_C max_C/2 max_C]);
        colormap(clrmp);
        colorbar;

       for i = 1 : length(ytic)      % draw anatomical boundaries

            jump = ROI_info{ ytic(i), 5};

            if iCondition == 1 || iCondition == 2

                plot([0.5+0 nClusters],0.5+jump*[1 1],'k-');

            else

                plot(0.5+jump*[1 1],0.5+[0 nClusters],'k-');

            end

            if i>1
                ytic_jump(i) = 0.5 * ( ROI_info{ ytic(i-1), 5} + ROI_info{ ytic(i), 5} );
            else
                ytic_jump(i) = 0.5 * ( 1 + ROI_info{ ytic(i), 5} );
            end
       end

        for i = 1 : length(group(iGroup).idx)      % draw anatomical boundaries

            if i>1
                xtic_jump(i) = 0.5 * ( ROI_info_v2{ i, 5} + ROI_info_v2{ i, 5} );
                jump = ROI_info_v2{ i, 5};
            else
                
                if iGroup > 1
                    
                    xtic_jump(i) = 0.5 * ( 1 + ROI_info_v2{ i, 5});
                    jump = ROI_info_v2{ i, 5};
               
                else
                
                    xtic_jump(i) = 0.5 * ( 1 + ROI_info_v2{ i, 5});
                    jump = ROI_info_v2{ i, 5 };
                    
                end
                
            end

            if iCondition == 1 || iCondition == 2

                plot([0.5+0 nClusters],0.5+jump*[1 1],'k-');

            else

                plot(0.5+jump*[1 1],0.5+[0 nClusters],'k-');

            end

        end
        hold off;

        if iCondition == 1 || iCondition == 2

            ax = gca;
            %axis 'square';
            %axis([ 0.5 nClusters+0.5 0.5 nClusters+0.5] );
            %set(ax,'XTick', [] );
            set(ax,'XTick',ytic_jump);
            set(ax,'XTickLabel',yticlb);
            %rotateXLabels( ax, 45 );

            set(ax,'YTick',xtic_jump);
            set(ax,'YTickLabel',xticlb);
            %rotateYLabels( ax, 45 );

            set(ax,'FontSize',5);

            xticklabel_rotate([],90,[],'Fontsize',5);

        else

            ax = gca;
            %axis 'square';
            %axis([ 0.5 nClusters+0.5 0.5 nClusters+0.5] );
            %set(ax,'XTick', [] );
            set(ax,'XTick',xtic_jump);
            set(ax,'XTickLabel',xticlb);
            %rotateXLabels( ax, 45 );

            set(ax,'YTick',ytic_jump);
            set(ax,'YTickLabel',yticlb);
            %rotateYLabels( ax, 45 );

            set(ax,'FontSize',5);

            xticklabel_rotate([],90,[],'Fontsize',5);

        end

        print(f,'-depsc',strcat('FC_Voxel_AAL_ROI_kmeans_Contrast_Clusters','-','Mean-Contrast-',label,'-',group_label{iGroup},'-',direction,'.eps'));

    end

end

end

function getMeaningfulClustersWithChangeBasedOnFC

clear all

load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-Contrast.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Functional-Parcels.mat');
% load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');
load('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters-Mean-Contrast.mat');

% load DTI
prefix = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION';
sufix = 'preprocessed\B0-DTI\6.forbedpost.bedpostX\Und_PRE_xx.mat';

DTI(1) = load(strcat(prefix,'\','SUBJECT-1-22-10-2015','\',sufix));
DTI(2) = load(strcat(prefix,'\','SUBJECT-2-26-10-2015','\',sufix));
DTI(3) = load(strcat(prefix,'\','SUBJECT-3-3-11-2015','\',sufix));
DTI(4) = load(strcat(prefix,'\','SUBJECT-4-2-11-2015','\',sufix));
DTI(5) = load(strcat(prefix,'\','SUBJECT-5-2-11-2015','\',sufix));
DTI(6) = load(strcat(prefix,'\','SUBJECT-6-24-11-2015','\',sufix));
DTI(7) = load(strcat(prefix,'\','SUBJECT-7-14-01-2016','\',sufix));
DTI(8) = load(strcat(prefix,'\','SUBJECT-8-14-01-2016','\',sufix));

% compute Common DTI
ave_DTI = zeros(size(DTI(1).C));
for iDTI=1:8
    
    ave_DTI = ave_DTI + DTI(iDTI).C;
    
end
ave_DTI = ave_DTI ./ 8;

common_DTI = ave_DTI;
common_DTI(~(DTI(1).C & DTI(2).C & DTI(3).C & DTI(4).C & DTI(5).C & DTI(6).C & DTI(7).C & DTI(8).C)) = 0;

nArea = 7;

area(1).idx = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 69 70]; % frontal
area(2).idx = [27 28 29 30 31 32 33 34 35 36]; % rec_ins_cing
area(3).idx = [37 38 39 40 41 42]; % hc_amyg
area(4).idx = [43 44 45 46 47 48 49 50 51 52 53 54]; % occipital
area(5).idx = [57 58 59 60 61 62 63 64 65 66 67 68]; % parietal
area(6).idx = [55 56 79 80 81 82 83 84 85 86 87 88 89 90]; % temporal
area(7).idx = [71 72 73 74 75 76 77 78]; % subcortical

area(1).label = 'frontal';
area(2).label = 'rec_ins_cing';
area(3).label = 'hc_amyg';
area(4).label = 'occipital';
area(5).label = 'parietal';
area(6).label = 'temporal';
area(7).label = 'subcortical';

area(1).color = 'r';
area(2).color = 'b';
area(3).color = 'g';
area(4).color = 'y';
area(5).color = 'w';
area(6).color = 'm';
area(7).color = 'k';

nClusters = 758;

FC_contrast = contrast_attention_z;
% contrast = contrast_stimulus_z;
GC_contrast = Attention_Contrast.Z;

zcriterion = 2.3;

FC_contrast(isnan(FC_contrast)) = 0;
GC_contrast(isnan(GC_contrast)) = 0;

FC_contrast(find(FC_contrast>(-1)*zcriterion & FC_contrast<zcriterion)) = 0; 
GC_contrast(find(GC_contrast>(-1)*zcriterion & GC_contrast<zcriterion)) = 0; 

FC_contrast_direct = FC_contrast;
FC_contrast_indirect = FC_contrast;

GC_contrast_direct = GC_contrast;
GC_contrast_indirect = GC_contrast;

FC_contrast_direct(~(common_DTI)) = 0;
GC_contrast_direct(~(common_DTI)) = 0;

FC_contrast_indirect(find(common_DTI)) = 0;
GC_contrast_indirect(find(common_DTI)) = 0;

% FC_contrast(~(common_DTI & FC_contrast & GC_contrast)) = 0;
% GC_contrast(~(common_DTI & FC_contrast & GC_contrast)) = 0;

nROI = 90;
for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

all_idx_pos = [];

for iCluster=1:nClusters
    
    FC_cluster = squeeze(FC_contrast(iCluster,:));
    idx_pos = find(FC_cluster>0);
    
    all_idx_pos = [all_idx_pos, idx_pos];

end

unique_idx_pos = sort(unique(all_idx_pos));

MeaningFul{1,1} = 'Cluster';
MeaningFul{1,2} = 'ROI';

for iUnique=1:length(unique_idx_pos)
    
   for iROI=1:nROI
       
      if unique_idx_pos(iUnique) >= ROI_info{iROI,4} && unique_idx_pos(iUnique) <= ROI_info{iROI,5}
          
          MeaningFul{iUnique+1,1} = unique_idx_pos(iUnique);
          MeaningFul{iUnique+1,2} = ROI_info{iROI,2};
          
      end
       
   end
    
end

MeaningFul{1,3} = 'Nets';
MeaningFul{1,4} = 'Nets-seeds';

nNet = 8;

for iUnique=1:length(unique_idx_pos)
    
   FC_parcels = cell.empty;
   
   for iROI=1:nROI
       
       if unique_idx_pos(iUnique) >= ROI_info{iROI,4} && unique_idx_pos(iUnique) <= ROI_info{iROI,5}
           
           for iNet=1:nNet
           
                if ROI(iROI).clusters(unique_idx_pos(iUnique)-ROI_info{iROI,4}+1).net_count{iNet,2} > 0
                    
                    if isempty(FC_parcels)
                        
                        FC_parcels = ROI(iROI).clusters(unique_idx_pos(iUnique)-ROI_info{iROI,4}+1).net_count(iNet,1);
                        
                    else
                        
                        FC_parcels = strcat(FC_parcels(1),',',ROI(iROI).clusters(unique_idx_pos(iUnique)-ROI_info{iROI,4}+1).net_count(iNet,1));
                        
                    end
                    
                end
           
           end
           
           seed_count = ROI(iROI).clusters(unique_idx_pos(iUnique)-ROI_info{iROI,4}+1).seed_count;
           
       end
       
   end
   
   MeaningFul{iUnique+1,3} = FC_parcels(:);
   MeaningFul{iUnique+1,4} = seed_count(:);
    
end

MeaningFul{1,5} = 'FC:+d';
MeaningFul{1,6} = 'FC:+i';
MeaningFul{1,7} = 'FC:+d-clusters';
MeaningFul{1,8} = 'FC:+i-clusters';

for iUnique=1:length(unique_idx_pos)
    
    FC_cluster = squeeze(FC_contrast_direct(unique_idx_pos(iUnique),:));
    idx_pos = find(FC_cluster>0);

    FC_areas = cell.empty;
    for idx=1:length(idx_pos)
        
        for iROI=1:nROI
       
            if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}
          
                FC_areas = [FC_areas(:)' ROI_info{iROI,2}];
                         
            end
       
        end
        
    end

    if ~isempty(FC_areas)
        
        [u_FC_areas, c_FC_areas] = count_unique(FC_areas);
        
    else
        
        u_FC_areas = cell.empty;
        
    end
    
    MeaningFul{iUnique+1,5} = u_FC_areas(:);
    MeaningFul{iUnique+1,7} = idx_pos;
    
    FC_cluster = squeeze(FC_contrast_indirect(unique_idx_pos(iUnique),:));
    idx_pos = find(FC_cluster>0);

    FC_areas = cell.empty;
    for idx=1:length(idx_pos)
        
        for iROI=1:nROI
       
            if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}
          
                FC_areas = [FC_areas(:)' ROI_info{iROI,2}];
                         
            end
       
        end
        
    end

    if ~isempty(FC_areas)
        
        [u_FC_areas, c_FC_areas] = count_unique(FC_areas);
        
    else
        
        u_FC_areas = cell.empty;
        
    end
    
    MeaningFul{iUnique+1,6} = u_FC_areas(:);
    MeaningFul{iUnique+1,8} = idx_pos;
    
end

MeaningFul{1,9} = 'GC:+d';
MeaningFul{1,10} = 'GC:+i';
MeaningFul{1,11} = 'GC:+d-clusters';
MeaningFul{1,12} = 'GC:+i-clusters';

for iUnique=1:length(unique_idx_pos)
    
    GC_cluster = squeeze(GC_contrast_direct(unique_idx_pos(iUnique),:));
    idx_pos = find(GC_cluster>0);
    idx_neg = find(GC_cluster<0);

    GC_areas_pos = cell.empty;
    for idx=1:length(idx_pos)
        
        for iROI=1:nROI
       
            if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}
          
                GC_areas_pos = [GC_areas_pos(:)' ROI_info{iROI,2}];
                         
            end
       
        end
        
    end

    if ~isempty(GC_areas_pos)
        
        [u_GC_areas, c_GC_areas] = count_unique(GC_areas_pos);
        
    else
        
        u_GC_areas = cell.empty;
        
    end
    
    MeaningFul{iUnique+1,9} = u_GC_areas(:);
    MeaningFul{iUnique+1,11} = idx_pos;

    GC_cluster = squeeze(GC_contrast_indirect(unique_idx_pos(iUnique),:));
    idx_pos = find(GC_cluster>0);
    idx_neg = find(GC_cluster<0);

    GC_areas_pos = cell.empty;
    for idx=1:length(idx_pos)
        
        for iROI=1:nROI
       
            if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}
          
                GC_areas_pos = [GC_areas_pos(:)' ROI_info{iROI,2}];
                         
            end
       
        end
        
    end

    if ~isempty(GC_areas_pos)
        
        [u_GC_areas, c_GC_areas] = count_unique(GC_areas_pos);
    
    else
        
        u_GC_areas = cell.empty;
        
    end
    
    MeaningFul{iUnique+1,10} = u_GC_areas(:);
    MeaningFul{iUnique+1,12} = idx_pos;
    
end

MeaningFul{1,13} = 'Seeds';

for iUnique=1:length(unique_idx_pos)
    
    idx_cluster = 0;
    
    for iROI=1:nROI
        
        for iCluster=1:length(ROI(iROI).clusters)
            
            idx_cluster = idx_cluster + 1;
            
            if idx_cluster == unique_idx_pos(iUnique)
    
                MeaningFul{iUnique+1,13} = ROI(iROI).clusters(iCluster).seed_count;
                
            end
            
        end
        
    end
    
end

save('FC-Voxels-AAL-ROI-corr-KMeans-MeaningFul-FC.mat','MeaningFul');

% all_meaningful_AAL = MeaningFul(2:end,2);
% [u_all_meaningful_AAL c_all_meaningful_AAL] = count_unique(all_meaningful_AAL);
% 
% for iMean=1:length(u_all_meaningful_AAL)
%     
%     for iROI=1:nROI
%    
%         if strcmp(ROI_info{iROI,2},u_all_meaningful_AAL{iMean,1})
%             
%             u_all_meaningful_AAL{iMean,2} = ROI_info{iROI,3};
%             
%             for iArea=1:nArea
%                 
%                if ismember(iROI,area(iArea).idx)
%                    
%                    u_all_meaningful_AAL{iMean,3} = area(iArea).label;
%                    
%                end
%                
%             end
%     
%         end
%     
%     end
%     
% end
% 
% icomb_in = 0;
% icomb_di = 0;
% for iMean=1:(size(MeaningFul,1)-1)
%     
%    conn_in = MeaningFul{iMean+1,6};
%     
%    for iCon=1:length(conn_in)
%        
%        icomb_in = icomb_in + 1;
%        
%        all_connections_indirect(icomb_in,1) = MeaningFul{iMean+1,1};
%        all_connections_indirect(icomb_in,2) = conn_in(iCon);
%         
%    end
%    
%    conn_di = MeaningFul{iMean+1,5};
%     
%    for iCon=1:length(conn_di)
%    
%        icomb_di = icomb_di + 1;
%        
%        all_connections_direct(icomb_di,1) = MeaningFul{iMean+1,1};
%        all_connections_direct(icomb_di,2) = conn_di(iCon);
%      
%    end
%    
% end
% 
% all_AAL_mean = [];
% all_clusters_mean = [];
% for iMean=[10 20 30 86]
%     
%     all_AAL_mean = [all_AAL_mean; MeaningFul{iMean,3}];
%     all_AAL_mean = [all_AAL_mean; MeaningFul{iMean,4}];
%     all_AAL_mean = [all_AAL_mean; MeaningFul{iMean,7}];
%     all_AAL_mean = [all_AAL_mean; MeaningFul{iMean,8}];
%     
%     all_clusters_mean = [all_clusters_mean, MeaningFul{iMean,5}];
%     all_clusters_mean = [all_clusters_mean, MeaningFul{iMean,6}];
%     
% end
% 
% [u_all_AAL_mean, c_all_AAL_mean] = count_unique(all_AAL_mean);
% [u_all_clusters_mean, c_all_clusters_mean] = count_unique(all_clusters_mean);
% 
% for iMean=1:length(u_all_AAL_mean)
%     
%     for iROI=1:nROI
%    
%         if strcmp(ROI_info{iROI,2},u_all_AAL_mean{iMean,1})
%             
%             u_all_AAL_mean{iMean,2} = ROI_info{iROI,3};
%             
%             for iArea=1:nArea
%                 
%                if ismember(iROI,area(iArea).idx)
%                    
%                    u_all_AAL_mean{iMean,2} = area(iArea).label;
%                    
%                end
%                
%             end
%     
%         end
%     
%     end
%     
% end

end

function getMeaningfulClustersWithChangeBasedOnGranger

clear all

load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-Contrast.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Functional-Parcels.mat');
% load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');
load('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters-Mean-Contrast.mat');

% load DTI
prefix = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION';
sufix = 'preprocessed\B0-DTI\6.forbedpost.bedpostX\Und_PRE_xx.mat';

DTI(1) = load(strcat(prefix,'\','SUBJECT-1-22-10-2015','\',sufix));
DTI(2) = load(strcat(prefix,'\','SUBJECT-2-26-10-2015','\',sufix));
DTI(3) = load(strcat(prefix,'\','SUBJECT-3-3-11-2015','\',sufix));
DTI(4) = load(strcat(prefix,'\','SUBJECT-4-2-11-2015','\',sufix));
DTI(5) = load(strcat(prefix,'\','SUBJECT-5-2-11-2015','\',sufix));
DTI(6) = load(strcat(prefix,'\','SUBJECT-6-24-11-2015','\',sufix));
DTI(7) = load(strcat(prefix,'\','SUBJECT-7-14-01-2016','\',sufix));
DTI(8) = load(strcat(prefix,'\','SUBJECT-8-14-01-2016','\',sufix));

% compute Common DTI
ave_DTI = zeros(size(DTI(1).C));
for iDTI=1:8
    
    ave_DTI = ave_DTI + DTI(iDTI).C;
    
end
ave_DTI = ave_DTI ./ 8;

common_DTI = ave_DTI;
common_DTI(~(DTI(1).C & DTI(2).C & DTI(3).C & DTI(4).C & DTI(5).C & DTI(6).C & DTI(7).C & DTI(8).C)) = 0;

nArea = 7;

area(1).idx = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 69 70]; % frontal
area(2).idx = [27 28 29 30 31 32 33 34 35 36]; % rec_ins_cing
area(3).idx = [37 38 39 40 41 42]; % hc_amyg
area(4).idx = [43 44 45 46 47 48 49 50 51 52 53 54]; % occipital
area(5).idx = [57 58 59 60 61 62 63 64 65 66 67 68]; % parietal
area(6).idx = [55 56 79 80 81 82 83 84 85 86 87 88 89 90]; % temporal
area(7).idx = [71 72 73 74 75 76 77 78]; % subcortical

area(1).label = 'frontal';
area(2).label = 'rec_ins_cing';
area(3).label = 'hc_amyg';
area(4).label = 'occipital';
area(5).label = 'parietal';
area(6).label = 'temporal';
area(7).label = 'subcortical';

area(1).color = 'r';
area(2).color = 'b';
area(3).color = 'g';
area(4).color = 'y';
area(5).color = 'w';
area(6).color = 'm';
area(7).color = 'k';

nClusters = 758;

FC_contrast = contrast_attention_z;
% contrast = contrast_stimulus_z;
GC_contrast = Attention_Contrast.Z;

zcriterion = 2.3;

FC_contrast(isnan(FC_contrast)) = 0;
GC_contrast(isnan(GC_contrast)) = 0;

FC_contrast(find(FC_contrast>(-1)*zcriterion & FC_contrast<zcriterion)) = 0; 
GC_contrast(find(GC_contrast>(-1)*zcriterion & GC_contrast<zcriterion)) = 0; 

FC_contrast_direct = FC_contrast;
FC_contrast_indirect = FC_contrast;

GC_contrast_direct = GC_contrast;
GC_contrast_indirect = GC_contrast;

FC_contrast_direct(~(common_DTI)) = 0;
GC_contrast_direct(~(common_DTI)) = 0;

FC_contrast_indirect(find(common_DTI)) = 0;
GC_contrast_indirect(find(common_DTI)) = 0;

% FC_contrast(~(common_DTI & FC_contrast & GC_contrast)) = 0;
% GC_contrast(~(common_DTI & FC_contrast & GC_contrast)) = 0;

nROI = 90;
for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

all_idx_pos = [];

for iCluster=1:nClusters
    
    GC_cluster = squeeze(GC_contrast(iCluster,:));
    idx_pos = find(GC_cluster>0);
    
    all_idx_pos = [all_idx_pos, idx_pos];

end

unique_idx_pos = sort(unique(all_idx_pos));

MeaningFul{1,1} = 'Cluster';
MeaningFul{1,2} = 'ROI';

for iUnique=1:length(unique_idx_pos)
    
   for iROI=1:nROI
       
      if unique_idx_pos(iUnique) >= ROI_info{iROI,4} && unique_idx_pos(iUnique) <= ROI_info{iROI,5}
          
          MeaningFul{iUnique+1,1} = unique_idx_pos(iUnique);
          MeaningFul{iUnique+1,2} = ROI_info{iROI,2};
          
      end
       
   end
    
end

MeaningFul{1,3} = 'Nets';
MeaningFul{1,4} = 'Nets-seeds';

nNet = 8;

for iUnique=1:length(unique_idx_pos)
    
   FC_parcels = cell.empty;
   
   for iROI=1:nROI
       
       if unique_idx_pos(iUnique) >= ROI_info{iROI,4} && unique_idx_pos(iUnique) <= ROI_info{iROI,5}
           
           for iNet=1:nNet
           
                if ROI(iROI).clusters(unique_idx_pos(iUnique)-ROI_info{iROI,4}+1).net_count{iNet,2} > 0
                    
                    if isempty(FC_parcels)
                        
                        FC_parcels = ROI(iROI).clusters(unique_idx_pos(iUnique)-ROI_info{iROI,4}+1).net_count(iNet,1);
                        
                    else
                        
                        FC_parcels = strcat(FC_parcels(1),',',ROI(iROI).clusters(unique_idx_pos(iUnique)-ROI_info{iROI,4}+1).net_count(iNet,1));
                        
                    end
                    
                end
           
           end
           
           seed_count = ROI(iROI).clusters(unique_idx_pos(iUnique)-ROI_info{iROI,4}+1).seed_count;
           
       end
       
   end
   
   MeaningFul{iUnique+1,3} = FC_parcels(:);
   MeaningFul{iUnique+1,4} = seed_count(:);
    
end

MeaningFul{1,5} = 'FC:+d';
MeaningFul{1,6} = 'FC:+i';
MeaningFul{1,7} = 'FC:+d-clusters';
MeaningFul{1,8} = 'FC:+i-clusters';

for iUnique=1:length(unique_idx_pos)
    
    FC_cluster = squeeze(FC_contrast_direct(unique_idx_pos(iUnique),:));
    idx_pos = find(FC_cluster>0);

    FC_areas = cell.empty;
    for idx=1:length(idx_pos)
        
        for iROI=1:nROI
       
            if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}
          
                FC_areas = [FC_areas(:)' ROI_info{iROI,2}];
                         
            end
       
        end
        
    end

    if ~isempty(FC_areas)
        
        [u_FC_areas, c_FC_areas] = count_unique(FC_areas);
        
    else
        
        u_FC_areas = cell.empty;
        
    end
    
    MeaningFul{iUnique+1,5} = u_FC_areas(:);
    MeaningFul{iUnique+1,7} = idx_pos;
    
    FC_cluster = squeeze(FC_contrast_indirect(unique_idx_pos(iUnique),:));
    idx_pos = find(FC_cluster>0);

    FC_areas = cell.empty;
    for idx=1:length(idx_pos)
        
        for iROI=1:nROI
       
            if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}
          
                FC_areas = [FC_areas(:)' ROI_info{iROI,2}];
                         
            end
       
        end
        
    end

    if ~isempty(FC_areas)
        
        [u_FC_areas, c_FC_areas] = count_unique(FC_areas);
        
    else
        
        u_FC_areas = cell.empty;
        
    end
    
    MeaningFul{iUnique+1,6} = u_FC_areas(:);
    MeaningFul{iUnique+1,8} = idx_pos;
    
end

MeaningFul{1,9} = 'GC:+d';
MeaningFul{1,10} = 'GC:+i';
MeaningFul{1,11} = 'GC:+d-clusters';
MeaningFul{1,12} = 'GC:+i-clusters';

for iUnique=1:length(unique_idx_pos)
    
    GC_cluster = squeeze(GC_contrast_direct(unique_idx_pos(iUnique),:));
    idx_pos = find(GC_cluster>0);
    idx_neg = find(GC_cluster<0);

    GC_areas_pos = cell.empty;
    for idx=1:length(idx_pos)
        
        for iROI=1:nROI
       
            if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}
          
                GC_areas_pos = [GC_areas_pos(:)' ROI_info{iROI,2}];
                         
            end
       
        end
        
    end

    if ~isempty(GC_areas_pos)
        
        [u_GC_areas, c_GC_areas] = count_unique(GC_areas_pos);
        
    else
        
        u_GC_areas = cell.empty;
        
    end
    
    MeaningFul{iUnique+1,9} = u_GC_areas(:);
    MeaningFul{iUnique+1,11} = idx_pos;

    GC_cluster = squeeze(GC_contrast_indirect(unique_idx_pos(iUnique),:));
    idx_pos = find(GC_cluster>0);
    idx_neg = find(GC_cluster<0);

    GC_areas_pos = cell.empty;
    for idx=1:length(idx_pos)
        
        for iROI=1:nROI
       
            if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}
          
                GC_areas_pos = [GC_areas_pos(:)' ROI_info{iROI,2}];
                         
            end
       
        end
        
    end

    if ~isempty(GC_areas_pos)
        
        [u_GC_areas, c_GC_areas] = count_unique(GC_areas_pos);
    
    else
        
        u_GC_areas = cell.empty;
        
    end
    
    MeaningFul{iUnique+1,10} = u_GC_areas(:);
    MeaningFul{iUnique+1,12} = idx_pos;
    
end

MeaningFul{1,13} = 'Seeds';

for iUnique=1:length(unique_idx_pos)
    
    idx_cluster = 0;
    
    for iROI=1:nROI
        
        for iCluster=1:length(ROI(iROI).clusters)
            
            idx_cluster = idx_cluster + 1;
            
            if idx_cluster == unique_idx_pos(iUnique)
    
                MeaningFul{iUnique+1,13} = ROI(iROI).clusters(iCluster).seed_count;
                
            end
            
        end
        
    end
    
end

save('FC-Voxels-AAL-ROI-corr-KMeans-MeaningFul-Granger.mat','MeaningFul');

% all_meaningful_AAL = MeaningFul(2:end,2);
% [u_all_meaningful_AAL c_all_meaningful_AAL] = count_unique(all_meaningful_AAL);
% 
% for iMean=1:length(u_all_meaningful_AAL)
%     
%     for iROI=1:nROI
%    
%         if strcmp(ROI_info{iROI,2},u_all_meaningful_AAL{iMean,1})
%             
%             u_all_meaningful_AAL{iMean,2} = ROI_info{iROI,3};
%             
%             for iArea=1:nArea
%                 
%                if ismember(iROI,area(iArea).idx)
%                    
%                    u_all_meaningful_AAL{iMean,3} = area(iArea).label;
%                    
%                end
%                
%             end
%     
%         end
%     
%     end
%     
% end
% 
% icomb_in = 0;
% icomb_di = 0;
% for iMean=1:(size(MeaningFul,1)-1)
%     
%    conn_in = MeaningFul{iMean+1,6};
%     
%    for iCon=1:length(conn_in)
%        
%        icomb_in = icomb_in + 1;
%        
%        all_connections_indirect(icomb_in,1) = MeaningFul{iMean+1,1};
%        all_connections_indirect(icomb_in,2) = conn_in(iCon);
%         
%    end
%    
%    conn_di = MeaningFul{iMean+1,5};
%     
%    for iCon=1:length(conn_di)
%    
%        icomb_di = icomb_di + 1;
%        
%        all_connections_direct(icomb_di,1) = MeaningFul{iMean+1,1};
%        all_connections_direct(icomb_di,2) = conn_di(iCon);
%      
%    end
%    
% end
% 
% all_AAL_mean = [];
% all_clusters_mean = [];
% for iMean=[10 20 30 86]
%     
%     all_AAL_mean = [all_AAL_mean; MeaningFul{iMean,3}];
%     all_AAL_mean = [all_AAL_mean; MeaningFul{iMean,4}];
%     all_AAL_mean = [all_AAL_mean; MeaningFul{iMean,7}];
%     all_AAL_mean = [all_AAL_mean; MeaningFul{iMean,8}];
%     
%     all_clusters_mean = [all_clusters_mean, MeaningFul{iMean,5}];
%     all_clusters_mean = [all_clusters_mean, MeaningFul{iMean,6}];
%     
% end
% 
% [u_all_AAL_mean, c_all_AAL_mean] = count_unique(all_AAL_mean);
% [u_all_clusters_mean, c_all_clusters_mean] = count_unique(all_clusters_mean);
% 
% for iMean=1:length(u_all_AAL_mean)
%     
%     for iROI=1:nROI
%    
%         if strcmp(ROI_info{iROI,2},u_all_AAL_mean{iMean,1})
%             
%             u_all_AAL_mean{iMean,2} = ROI_info{iROI,3};
%             
%             for iArea=1:nArea
%                 
%                if ismember(iROI,area(iArea).idx)
%                    
%                    u_all_AAL_mean{iMean,2} = area(iArea).label;
%                    
%                end
%                
%             end
%     
%         end
%     
%     end
%     
% end

end

function getMeaningfulClustersWithChangeBasedOnFCGCAttentionStimulusOnly

clear all

% load Cluster information
load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-Contrast.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Functional-Parcels.mat');
load('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters-Mean-Contrast.mat');

% load DTI
prefix = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION';
sufix = 'preprocessed\B0-DTI\6.forbedpost.bedpostX\Und_PRE_xx.mat';

DTI(1) = load(strcat(prefix,'\','SUBJECT-1-22-10-2015','\',sufix));
DTI(2) = load(strcat(prefix,'\','SUBJECT-2-26-10-2015','\',sufix));
DTI(3) = load(strcat(prefix,'\','SUBJECT-3-3-11-2015','\',sufix));
DTI(4) = load(strcat(prefix,'\','SUBJECT-4-2-11-2015','\',sufix));
DTI(5) = load(strcat(prefix,'\','SUBJECT-5-2-11-2015','\',sufix));
DTI(6) = load(strcat(prefix,'\','SUBJECT-6-24-11-2015','\',sufix));
DTI(7) = load(strcat(prefix,'\','SUBJECT-7-14-01-2016','\',sufix));
DTI(8) = load(strcat(prefix,'\','SUBJECT-8-14-01-2016','\',sufix));

% compute Common DTI
ave_DTI = zeros(size(DTI(1).C));
for iDTI=1:8
    
    ave_DTI = ave_DTI + DTI(iDTI).C;
    
end
ave_DTI = ave_DTI ./ 8;

common_DTI = ave_DTI;
common_DTI(~(DTI(1).C & DTI(2).C & DTI(3).C & DTI(4).C & DTI(5).C & DTI(6).C & DTI(7).C & DTI(8).C)) = 0;

nClusters = 758;

% load Contrast - Attention

FC_contrast_Attention = contrast_attention_z;
GC_contrast_Attention = Attention_Contrast.Z;

zcriterion = 2.3;

FC_contrast_Attention(isnan(FC_contrast_Attention)) = 0;
GC_contrast_Attention(isnan(GC_contrast_Attention)) = 0;

FC_contrast_Attention(find(FC_contrast_Attention>(-1)*zcriterion & FC_contrast_Attention<zcriterion)) = 0; 
GC_contrast_Attention(find(GC_contrast_Attention>(-1)*zcriterion & GC_contrast_Attention<zcriterion)) = 0; 

% load Contrast - Stimulus

FC_contrast_Stimulus = contrast_stimulus_z;
GC_contrast_Stimulus = Stimulus_Contrast.Z;

zcriterion = 2.3;

FC_contrast_Stimulus(isnan(FC_contrast_Stimulus)) = 0;
GC_contrast_Stimulus(isnan(GC_contrast_Stimulus)) = 0;

FC_contrast_Stimulus(find(FC_contrast_Stimulus>(-1)*zcriterion & FC_contrast_Stimulus<zcriterion)) = 0; 
GC_contrast_Stimulus(find(GC_contrast_Stimulus>(-1)*zcriterion & GC_contrast_Stimulus<zcriterion)) = 0; 

% Attention and Stimulus - ONLY

FC_contrast_Attention_tmp = FC_contrast_Attention;
FC_contrast_Stimulus_tmp = FC_contrast_Stimulus;

FC_contrast_Attention(find(FC_contrast_Stimulus_tmp)) = 0;
FC_contrast_Stimulus(find(FC_contrast_Attention_tmp)) = 0;

% Direct and Indirect - Attention

FC_contrast_Attention_direct = FC_contrast_Attention;
FC_contrast_Attention_indirect = FC_contrast_Attention;

GC_contrast_Attention_direct = GC_contrast_Attention;
GC_contrast_Attention_indirect = GC_contrast_Attention;

FC_contrast_Attention_direct(~(common_DTI)) = 0;
GC_contrast_Attention_direct(~(common_DTI)) = 0;

FC_contrast_Attention_indirect(find(common_DTI)) = 0;
GC_contrast_Attention_indirect(find(common_DTI)) = 0;

% Direct and Indirect - Stimulus

FC_contrast_Stimulus_direct = FC_contrast_Stimulus;
FC_contrast_Stimulus_indirect = FC_contrast_Stimulus;

GC_contrast_Stimulus_direct = GC_contrast_Stimulus;
GC_contrast_Stimulus_indirect = GC_contrast_Stimulus;

FC_contrast_Stimulus_direct(~(common_DTI)) = 0;
GC_contrast_Stimulus_direct(~(common_DTI)) = 0;

FC_contrast_Stimulus_indirect(find(common_DTI)) = 0;
GC_contrast_Stimulus_indirect(find(common_DTI)) = 0;

% ROI and Clusters information

nROI = 90;
for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

% MeaningFul Clusters - ATTENTION

all_idx_pos = [];
all_idx_neg = [];

for iCluster=1:nClusters
    
    FC_cluster = squeeze(FC_contrast_Attention(iCluster,:));
    idx_pos = find(FC_cluster>0);
    idx_neg = find(FC_cluster<0);
    
    all_idx_pos = [all_idx_pos, idx_pos];
    all_idx_neg = [all_idx_neg, idx_neg];
    
%     GC_cluster = squeeze(GC_contrast_Attention(iCluster,:));
%     idx_pos = find(GC_cluster>0);
%     
%     all_idx_pos = [all_idx_pos, idx_pos];

end

for iChange=1:2
    
    if iChange == 1; unique_idx = sort(unique(all_idx_pos)); change_label = 'Increase'; end
    if iChange == 2; unique_idx = sort(unique(all_idx_neg)); change_label = 'Decrease'; end

    MeaningFul{1,1} = 'Cluster';
    MeaningFul{1,2} = 'ROI';

    for iUnique=1:length(unique_idx)

       for iROI=1:nROI

          if unique_idx(iUnique) >= ROI_info{iROI,4} && unique_idx(iUnique) <= ROI_info{iROI,5}

              MeaningFul{iUnique+1,1} = unique_idx(iUnique);
              MeaningFul{iUnique+1,2} = ROI_info{iROI,2};

          end

       end

    end

    MeaningFul{1,3} = 'Nets';
    MeaningFul{1,4} = 'Nets-seeds';

    nNet = 8;

    for iUnique=1:length(unique_idx)

       FC_parcels = cell.empty;

       for iROI=1:nROI

           if unique_idx(iUnique) >= ROI_info{iROI,4} && unique_idx(iUnique) <= ROI_info{iROI,5}

               for iNet=1:nNet

                    if ROI(iROI).clusters(unique_idx(iUnique)-ROI_info{iROI,4}+1).net_count{iNet,2} > 0

                        if isempty(FC_parcels)

                            FC_parcels = ROI(iROI).clusters(unique_idx(iUnique)-ROI_info{iROI,4}+1).net_count(iNet,1);

                        else

                            FC_parcels = strcat(FC_parcels(1),',',ROI(iROI).clusters(unique_idx(iUnique)-ROI_info{iROI,4}+1).net_count(iNet,1));

                        end

                    end

               end

               seed_count = ROI(iROI).clusters(unique_idx(iUnique)-ROI_info{iROI,4}+1).seed_count;

           end

       end

       MeaningFul{iUnique+1,3} = FC_parcels(:);
       MeaningFul{iUnique+1,4} = seed_count(:);

    end

    MeaningFul{1,5} = 'FC:+d';
    MeaningFul{1,6} = 'FC:+i';
    MeaningFul{1,7} = 'FC:+d-clusters';
    MeaningFul{1,8} = 'FC:+i-clusters';

    for iUnique=1:length(unique_idx)

        FC_cluster = squeeze(FC_contrast_Attention_direct(unique_idx(iUnique),:));
        if iChange == 1; idx_change = find(FC_cluster>0); end
        if iChange == 2; idx_change = find(FC_cluster<0); end

        FC_areas = cell.empty;
        for idx=1:length(idx_change)

            for iROI=1:nROI

                if idx_change(idx) >= ROI_info{iROI,4} && idx_change(idx) <= ROI_info{iROI,5}

                    FC_areas = [FC_areas(:)' ROI_info{iROI,2}];

                end

            end

        end

        if ~isempty(FC_areas)

            [u_FC_areas, c_FC_areas] = count_unique(FC_areas);

        else

            u_FC_areas = cell.empty;

        end

        MeaningFul{iUnique+1,5} = u_FC_areas(:);
        MeaningFul{iUnique+1,7} = idx_change;

        FC_cluster = squeeze(FC_contrast_Attention_indirect(unique_idx(iUnique),:));
        if iChange == 1; idx_change = find(FC_cluster>0); end
        if iChange == 2; idx_change = find(FC_cluster<0); end

        FC_areas = cell.empty;
        for idx=1:length(idx_change)

            for iROI=1:nROI

                if idx_change(idx) >= ROI_info{iROI,4} && idx_change(idx) <= ROI_info{iROI,5}

                    FC_areas = [FC_areas(:)' ROI_info{iROI,2}];

                end

            end

        end

        if ~isempty(FC_areas)

            [u_FC_areas, c_FC_areas] = count_unique(FC_areas);

        else

            u_FC_areas = cell.empty;

        end

        MeaningFul{iUnique+1,6} = u_FC_areas(:);
        MeaningFul{iUnique+1,8} = idx_change;

    end

    MeaningFul{1,9} = 'FROM:GC:+d';
    MeaningFul{1,10} = 'FROM:GC:+i';
    MeaningFul{1,11} = 'FROM:GC:+d-clusters';
    MeaningFul{1,12} = 'FROM:GC:+i-clusters';

    for iUnique=1:length(unique_idx)

        GC_cluster = squeeze(GC_contrast_Attention_direct(unique_idx(iUnique),:));
        idx_pos = find(GC_cluster>0);
        idx_neg = find(GC_cluster<0);

        GC_areas_pos = cell.empty;
        for idx=1:length(idx_pos)

            for iROI=1:nROI

                if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}

                    GC_areas_pos = [GC_areas_pos(:)' ROI_info{iROI,2}];

                end

            end

        end

        if ~isempty(GC_areas_pos)

            [u_GC_areas, c_GC_areas] = count_unique(GC_areas_pos);

        else

            u_GC_areas = cell.empty;

        end

        MeaningFul{iUnique+1,9} = u_GC_areas(:);
        MeaningFul{iUnique+1,11} = idx_pos;

        GC_cluster = squeeze(GC_contrast_Attention_indirect(unique_idx(iUnique),:));
        idx_pos = find(GC_cluster>0);
        idx_neg = find(GC_cluster<0);

        GC_areas_pos = cell.empty;
        for idx=1:length(idx_pos)

            for iROI=1:nROI

                if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}

                    GC_areas_pos = [GC_areas_pos(:)' ROI_info{iROI,2}];

                end

            end

        end

        if ~isempty(GC_areas_pos)

            [u_GC_areas, c_GC_areas] = count_unique(GC_areas_pos);

        else

            u_GC_areas = cell.empty;

        end

        MeaningFul{iUnique+1,10} = u_GC_areas(:);
        MeaningFul{iUnique+1,12} = idx_pos;

    end

    MeaningFul{1,13} = 'TO:GC:+d';
    MeaningFul{1,14} = 'TO:GC:+i';
    MeaningFul{1,15} = 'TO:GC:+d-clusters';
    MeaningFul{1,16} = 'TO:GC:+i-clusters';

    for iUnique=1:length(unique_idx)

        GC_cluster = squeeze(GC_contrast_Attention_direct(:,unique_idx(iUnique)));
        idx_pos = find(GC_cluster>0);
        idx_neg = find(GC_cluster<0);

        GC_areas_pos = cell.empty;
        for idx=1:length(idx_pos)

            for iROI=1:nROI

                if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}

                    GC_areas_pos = [GC_areas_pos(:)' ROI_info{iROI,2}];

                end

            end

        end

        if ~isempty(GC_areas_pos)

            [u_GC_areas, c_GC_areas] = count_unique(GC_areas_pos);

        else

            u_GC_areas = cell.empty;

        end

        MeaningFul{iUnique+1,13} = u_GC_areas(:);
        MeaningFul{iUnique+1,15} = idx_pos;

        GC_cluster = squeeze(GC_contrast_Attention_indirect(:,unique_idx(iUnique)));
        idx_pos = find(GC_cluster>0);
        idx_neg = find(GC_cluster<0);

        GC_areas_pos = cell.empty;
        for idx=1:length(idx_pos)

            for iROI=1:nROI

                if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}

                    GC_areas_pos = [GC_areas_pos(:)' ROI_info{iROI,2}];

                end

            end

        end

        if ~isempty(GC_areas_pos)

            [u_GC_areas, c_GC_areas] = count_unique(GC_areas_pos);

        else

            u_GC_areas = cell.empty;

        end

        MeaningFul{iUnique+1,14} = u_GC_areas(:);
        MeaningFul{iUnique+1,16} = idx_pos;

    end

    MeaningFul{1,17} = 'Seeds';

    for iUnique=1:length(unique_idx)

        idx_cluster = 0;

        seeds = [];

        for iROI=1:nROI

            for iCluster=1:length(ROI(iROI).clusters)

                idx_cluster = idx_cluster + 1;

                if idx_cluster == unique_idx(iUnique)

                    %MeaningFul{iUnique+1,17} = ROI(iROI).clusters(iCluster).seed_count;

                    idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

                    Voxels = zeros(length(idx_voxels),3);

                    for iVoxel=1:length(idx_voxels)

                        [Voxels(iVoxel,1) Voxels(iVoxel,2) Voxels(iVoxel,3)] = ind2sub([91 109 91],idx_voxels(iVoxel));

                    end

                    seeds = seeIfAVoxelIsInsideASeed(Voxels);

                end

            end

        end

        MeaningFul{iUnique+1,17} = seeds;

    end

    MeaningFul_Attention = MeaningFul;
    clear MeaningFul

    % MeaningFul Clusters - STIMULUS

    all_idx = [];

    for iCluster=1:nClusters

        FC_cluster = squeeze(FC_contrast_Stimulus(iCluster,:));
        if iChange == 1; idx_change = find(FC_cluster>0); end
        if iChange == 2; idx_change = find(FC_cluster<0); end

        all_idx = [all_idx, idx_change];

%         GC_cluster = squeeze(GC_contrast_Stimulus(iCluster,:));
%         idx_pos = find(GC_cluster>0);
% 
%         all_idx_pos = [all_idx_pos, idx_pos];

    end

    unique_idx = sort(unique(all_idx));

    MeaningFul{1,1} = 'Cluster';
    MeaningFul{1,2} = 'ROI';

    for iUnique=1:length(unique_idx)

       for iROI=1:nROI

          if unique_idx(iUnique) >= ROI_info{iROI,4} && unique_idx(iUnique) <= ROI_info{iROI,5}

              MeaningFul{iUnique+1,1} = unique_idx(iUnique);
              MeaningFul{iUnique+1,2} = ROI_info{iROI,2};

          end

       end

    end

    MeaningFul{1,3} = 'Nets';
    MeaningFul{1,4} = 'Nets-seeds';

    nNet = 8;

    for iUnique=1:length(unique_idx)

       FC_parcels = cell.empty;

       for iROI=1:nROI

           if unique_idx(iUnique) >= ROI_info{iROI,4} && unique_idx(iUnique) <= ROI_info{iROI,5}

               for iNet=1:nNet

                    if ROI(iROI).clusters(unique_idx(iUnique)-ROI_info{iROI,4}+1).net_count{iNet,2} > 0

                        if isempty(FC_parcels)

                            FC_parcels = ROI(iROI).clusters(unique_idx(iUnique)-ROI_info{iROI,4}+1).net_count(iNet,1);

                        else

                            FC_parcels = strcat(FC_parcels(1),',',ROI(iROI).clusters(unique_idx(iUnique)-ROI_info{iROI,4}+1).net_count(iNet,1));

                        end

                    end

               end

               seed_count = ROI(iROI).clusters(unique_idx(iUnique)-ROI_info{iROI,4}+1).seed_count;

           end

       end

       MeaningFul{iUnique+1,3} = FC_parcels(:);
       MeaningFul{iUnique+1,4} = seed_count(:);

    end

    MeaningFul{1,5} = 'FC:+d';
    MeaningFul{1,6} = 'FC:+i';
    MeaningFul{1,7} = 'FC:+d-clusters';
    MeaningFul{1,8} = 'FC:+i-clusters';

    for iUnique=1:length(unique_idx)

        FC_cluster = squeeze(FC_contrast_Stimulus_direct(unique_idx(iUnique),:));
        if iChange == 1; idx_change = find(FC_cluster>0); end
        if iChange == 2; idx_change = find(FC_cluster<0); end

        FC_areas = cell.empty;
        for idx=1:length(idx_change)

            for iROI=1:nROI

                if idx_change(idx) >= ROI_info{iROI,4} && idx_change(idx) <= ROI_info{iROI,5}

                    FC_areas = [FC_areas(:)' ROI_info{iROI,2}];

                end

            end

        end

        if ~isempty(FC_areas)

            [u_FC_areas, c_FC_areas] = count_unique(FC_areas);

        else

            u_FC_areas = cell.empty;

        end

        MeaningFul{iUnique+1,5} = u_FC_areas(:);
        MeaningFul{iUnique+1,7} = idx_change;

        FC_cluster = squeeze(FC_contrast_Stimulus_indirect(unique_idx(iUnique),:));
        if iChange == 1; idx_change = find(FC_cluster>0); end
        if iChange == 2; idx_change = find(FC_cluster<0); end

        FC_areas = cell.empty;
        for idx=1:length(idx_change)

            for iROI=1:nROI

                if idx_change(idx) >= ROI_info{iROI,4} && idx_change(idx) <= ROI_info{iROI,5}

                    FC_areas = [FC_areas(:)' ROI_info{iROI,2}];

                end

            end

        end

        if ~isempty(FC_areas)

            [u_FC_areas, c_FC_areas] = count_unique(FC_areas);

        else

            u_FC_areas = cell.empty;

        end

        MeaningFul{iUnique+1,6} = u_FC_areas(:);
        MeaningFul{iUnique+1,8} = idx_change;

    end

    MeaningFul{1,9} = 'FROM:GC:+d';
    MeaningFul{1,10} = 'FROM:GC:+i';
    MeaningFul{1,11} = 'FROM:GC:+d-clusters';
    MeaningFul{1,12} = 'FROM:GC:+i-clusters';

    for iUnique=1:length(unique_idx)

        GC_cluster = squeeze(GC_contrast_Stimulus_direct(unique_idx(iUnique),:));
        idx_pos = find(GC_cluster>0);
        idx_neg = find(GC_cluster<0);

        GC_areas_pos = cell.empty;
        for idx=1:length(idx_pos)

            for iROI=1:nROI

                if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}

                    GC_areas_pos = [GC_areas_pos(:)' ROI_info{iROI,2}];

                end

            end

        end

        if ~isempty(GC_areas_pos)

            [u_GC_areas, c_GC_areas] = count_unique(GC_areas_pos);

        else

            u_GC_areas = cell.empty;

        end

        MeaningFul{iUnique+1,9} = u_GC_areas(:);
        MeaningFul{iUnique+1,11} = idx_pos;

        GC_cluster = squeeze(GC_contrast_Stimulus_indirect(unique_idx(iUnique),:));
        idx_pos = find(GC_cluster>0);
        idx_neg = find(GC_cluster<0);

        GC_areas_pos = cell.empty;
        for idx=1:length(idx_pos)

            for iROI=1:nROI

                if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}

                    GC_areas_pos = [GC_areas_pos(:)' ROI_info{iROI,2}];

                end

            end

        end

        if ~isempty(GC_areas_pos)

            [u_GC_areas, c_GC_areas] = count_unique(GC_areas_pos);

        else

            u_GC_areas = cell.empty;

        end

        MeaningFul{iUnique+1,10} = u_GC_areas(:);
        MeaningFul{iUnique+1,12} = idx_pos;

    end

    MeaningFul{1,13} = 'TO:GC:+d';
    MeaningFul{1,14} = 'TO:GC:+i';
    MeaningFul{1,15} = 'TO:GC:+d-clusters';
    MeaningFul{1,16} = 'TO:GC:+i-clusters';

    for iUnique=1:length(unique_idx)

        GC_cluster = squeeze(GC_contrast_Stimulus_direct(:,unique_idx(iUnique)));
        idx_pos = find(GC_cluster>0);
        idx_neg = find(GC_cluster<0);

        GC_areas_pos = cell.empty;
        for idx=1:length(idx_pos)

            for iROI=1:nROI

                if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}

                    GC_areas_pos = [GC_areas_pos(:)' ROI_info{iROI,2}];

                end

            end

        end

        if ~isempty(GC_areas_pos)

            [u_GC_areas, c_GC_areas] = count_unique(GC_areas_pos);

        else

            u_GC_areas = cell.empty;

        end

        MeaningFul{iUnique+1,13} = u_GC_areas(:);
        MeaningFul{iUnique+1,15} = idx_pos;

        GC_cluster = squeeze(GC_contrast_Stimulus_indirect(:,unique_idx(iUnique)));
        idx_pos = find(GC_cluster>0);
        idx_neg = find(GC_cluster<0);

        GC_areas_pos = cell.empty;
        for idx=1:length(idx_pos)

            for iROI=1:nROI

                if idx_pos(idx) >= ROI_info{iROI,4} && idx_pos(idx) <= ROI_info{iROI,5}

                    GC_areas_pos = [GC_areas_pos(:)' ROI_info{iROI,2}];

                end

            end

        end

        if ~isempty(GC_areas_pos)

            [u_GC_areas, c_GC_areas] = count_unique(GC_areas_pos);

        else

            u_GC_areas = cell.empty;

        end

        MeaningFul{iUnique+1,14} = u_GC_areas(:);
        MeaningFul{iUnique+1,16} = idx_pos;

    end

    MeaningFul{1,17} = 'Seeds';

    for iUnique=1:length(unique_idx)

        idx_cluster = 0;

        seeds = [];

        for iROI=1:nROI

            for iCluster=1:length(ROI(iROI).clusters)

                idx_cluster = idx_cluster + 1;

                if idx_cluster == unique_idx(iUnique)

                    %MeaningFul{iUnique+1,17} = ROI(iROI).clusters(iCluster).seed_count;

                    idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

                    Voxels = zeros(length(idx_voxels),3);

                    for iVoxel=1:length(idx_voxels)

                        [Voxels(iVoxel,1) Voxels(iVoxel,2) Voxels(iVoxel,3)] = ind2sub([91 109 91],idx_voxels(iVoxel));

                    end

                    seeds = seeIfAVoxelIsInsideASeed(Voxels);

                end

            end

        end

        MeaningFul{iUnique+1,17} = seeds;

    end

    MeaningFul_Stimulus = MeaningFul;
    clear MeaningFul

    save(strcat('FC-Voxels-AAL-ROI-corr-KMeans-MeaningFul-FC-Att-Stim-Only-',change_label,'.mat'),'MeaningFul_Attention','MeaningFul_Stimulus');

end

end

function plotMeaningfulClustersWithChangeBasedOnFCGCAttentionOnly

load('FC-Voxels-AAL-ROI-corr-KMeans-MeaningFul-FC-GC-Att-Stim-Only.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Functional-Parcels.mat');

% ROI and Clusters information

nROI = 90;
for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

nClusters = 758;
nROI = 90;

angle = linspace(0,2*pi,nClusters+1);
x = cos(angle);
y = sin(angle);

plot(x,y,'k');
hold on

for iROI=1:nROI
   
    position = ROI_info{iROI,4} + round((ROI_info{iROI,5} - ROI_info{iROI,4})/2);
    
    t = text(x(position),y(position),ROI_info{iROI,2});
    set(t,'Rotation',angle(position)*(360/2*pi));
    
end

end

function plotFCAndGrangerAtGroupLevelContrastOnlySeeds

load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-Contrast.mat');
load(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast','.mat'));
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-MeaningFul-FC-Att-Stim-Only.mat');

do = 'FC';
% do = 'GC';

if strcmp(do,'FC'); attention_z = contrast_attention_z; stimulus_z = contrast_stimulus_z; Analysis_Label = 'FC'; end
if strcmp(do,'GC'); attention_z = Attention_Contrast.Z; stimulus_z = Stimulus_Contrast.Z; Analysis_Label = 'Granger'; end

nRun = 32;
nTotalClusters = 758;
nClusters = 758;
nTR = 150;
pcriterion = 0.01;
nROI = 90;

fs = 5;

nROI = 90;
for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

cluster_has_seeds = zeros(1,nClusters);
idx_cluster = 0;
seeds_labels = cell.empty;
iiCluster = 0;
for iROI=1:nROI
        
    for iCluster=1:length(ROI(iROI).clusters)

        idx_cluster = idx_cluster + 1;

        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

        Voxels = zeros(length(idx_voxels),3);

        for iVoxel=1:length(idx_voxels)

            [Voxels(iVoxel,1) Voxels(iVoxel,2) Voxels(iVoxel,3)] = ind2sub([91 109 91],idx_voxels(iVoxel));

        end

        seeds = seeIfAVoxelIsInsideASeed(Voxels);
        
        ROI_info{iROI,6} = seeds;
        
        if ~isempty(seeds)
            
            cluster_has_seeds(idx_cluster) = 1; 
        
            if iscell(seeds)
               
                this_seed = seeds{1};
                for i=2:length(seeds)
                    this_seed = strcat(this_seed,' / ',seeds{i});
                end
            else
                this_seed = seeds;
            end
            
            iiCluster = iiCluster + 1;
            
            %check_seeds{iiCluster,1} = seeds;
            
            seeds_labels{iiCluster} = this_seed;
            
        end

    end
        
end

% zcriterion = 1.6;
zcriterion = 2.3;

for iCondition=1:4
    
    attention_z(find(attention_z>(-1)*zcriterion & attention_z<zcriterion)) = 0; 
    stimulus_z(find(stimulus_z>(-1)*zcriterion & stimulus_z<zcriterion)) = 0; 
    
    if iCondition == 1; rho = attention_z; label = 'Attention'; end
    if iCondition == 2; rho = stimulus_z; label = 'Stimulus'; end
    
    if iCondition == 3; rho = attention_z; label = 'Attention-Only'; end
    if iCondition == 4; rho = stimulus_z; label = 'Stimulus-Only'; end
    
    rho(isnan(rho)) = 0;
    
    if iCondition == 3

        rho(find(stimulus_z)) = 0; 
    
    end
    
    if iCondition == 4 
        
        rho(find(attention_z)) = 0; 
    
    end
    
    rho(find(~(cluster_has_seeds)),:) = [];
    rho(:,find(~(cluster_has_seeds))) = [];
    
    f = figure;
    
    imagesc(rho);
    hold on
   
    title(label);

    clrmp = colormap('jet');
    clrmp(65,:) = clrmp(64,:);
    clrmp(33,:) = [1 1 1];

    min_C = -3;
    max_C = 3;

    %colorbar;
    caxis([min_C max_C]);
    %colorbar('Ticks',[min_C max_C/2 max_C]);
    colormap(clrmp);
    colorbar;
    
    hold off;

    ax = gca;
    set(ax,'XTick',1:length(seeds_labels));
    set(ax,'YTick',1:length(seeds_labels));
    set(ax,'XTickLabel',seeds_labels);
    set(ax,'YTickLabel',seeds_labels);
    set(ax,'FontSize',fs);

    xticklabel_rotate([],90,[],'Fontsize',fs);

    print(f,'-depsc',strcat('FC_Voxel_AAL_ROI_kmeans_',Analysis_Label,'_Clusters','-','Mean-Contrast-',label,'-Seeds.eps'));

end

end

function plotFCAndGrangerAtGroupLevelContrastOnlySeedsSortedPerNetwork

load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-Contrast.mat');
load(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast','.mat'));
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-MeaningFul-FC-Att-Stim-Only.mat');

% do = 'FC';
do = 'GC';

if strcmp(do,'FC'); attention_z = contrast_attention_z; stimulus_z = contrast_stimulus_z; Analysis_Label = 'FC'; end
if strcmp(do,'GC'); attention_z = Attention_Contrast.Z; stimulus_z = Stimulus_Contrast.Z; Analysis_Label = 'Granger'; end

nRun = 32;
nTotalClusters = 758;
nClusters = 758;
nTR = 150;
pcriterion = 0.01;
nROI = 90;

fs = 5;

nROI = 90;
for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

cluster_has_seeds = zeros(1,nClusters);
idx_cluster = 0;
seeds_labels = cell.empty;
iiCluster = 0;
iiseed = 0;
for iROI=1:nROI
        
    for iCluster=1:length(ROI(iROI).clusters)

        idx_cluster = idx_cluster + 1;

        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

        Voxels = zeros(length(idx_voxels),3);

        for iVoxel=1:length(idx_voxels)

            [Voxels(iVoxel,1) Voxels(iVoxel,2) Voxels(iVoxel,3)] = ind2sub([91 109 91],idx_voxels(iVoxel));

        end

        seeds = seeIfAVoxelIsInsideASeed(Voxels);
        
        ROI_info{iROI,6} = seeds;
        
        if ~isempty(seeds)
            
            cluster_has_seeds(idx_cluster) = 1; 
        
            if iscell(seeds)
               
                iiseed = iiseed + 1;
                
                this_seed = seeds{1};
                all_seeds{iiseed} = this_seed;
                which_cluster_each_seed(iiseed) = idx_cluster;
                for i=2:length(seeds)
                    
                    iiseed = iiseed + 1;
                    this_seed = strcat(this_seed,' / ',seeds{i});
                    all_seeds{iiseed} = seeds{i};
                    which_cluster_each_seed(iiseed) = idx_cluster;
                    
                end
                
            else
                
                iiseed = iiseed + 1;
                this_seed = seeds;
                all_seeds{iiseed} = this_seed;
                which_cluster_each_seed(iiseed) = idx_cluster;
                
            end
            
            iiCluster = iiCluster + 1;
            
            %check_seeds{iiCluster,1} = seeds;
            
            seeds_labels{iiCluster} = this_seed;
            
        end

    end
        
end

[labels,idx] = sort(seeds_labels);
[new_labels,new_idx] = sort(all_seeds);
which_cluster_each_seed = which_cluster_each_seed(new_idx);

corbetta_networks_labels = {'DAN' 'VAN' 'VIS'};

corbetta_seeds = {'FEF' 'PIPS' 'SPL' 'r-P' 'PRECU' 'TPJ' 'V3' 'V4' 'V7' 'V8'};

idx_out = [];
for iSeed=1:length(which_cluster_each_seed)
   
    seed_name = new_labels{iSeed};
    
    this_seed_is_out = 1;
    
    for iNet=1:length(corbetta_networks_labels)
       
        if strcmp(seed_name(1:3),corbetta_networks_labels{iNet})
            
            for iMatch=1:length(corbetta_seeds)
            
                if strfind(lower(seed_name),lower(corbetta_seeds{iMatch}))
                    
                    this_seed_is_out = 0;
            
                end
                
            end
            
        end
        
    end
    
    if this_seed_is_out
        
        idx_out = [idx_out, iSeed];
        
    end
    
end

which_cluster_each_seed(idx_out) = [];
new_labels(idx_out) = [];

% zcriterion = 1.6;
zcriterion = 2.3;

for iCondition=1:4
    
%     attention_z(find(attention_z>(-1)*zcriterion & attention_z<zcriterion)) = 0; 
%     stimulus_z(find(stimulus_z>(-1)*zcriterion & stimulus_z<zcriterion)) = 0; 

    attention_sig = attention_z;
    stimulus_sig = stimulus_z;
    attention_sig(find(attention_z>(-1)*zcriterion & attention_z<zcriterion)) = 0; 
    stimulus_sig(find(stimulus_z>(-1)*zcriterion & stimulus_z<zcriterion)) = 0; 
    
    if iCondition == 1; rho = attention_z; rho_sig = attention_sig; label = 'Attention'; end
    if iCondition == 2; rho = stimulus_z; rho_sig = stimulus_sig; label = 'Stimulus'; end
    
    if iCondition == 3; rho = attention_z; rho_sig = attention_sig; label = 'Attention-Only'; end
    if iCondition == 4; rho = stimulus_z; rho_sig = stimulus_sig; label = 'Stimulus-Only'; end
    
    rho(isnan(rho)) = 0;
    rho_sig(isnan(rho_sig)) = 0;
    
    attention_sig(find(attention_sig)) = 1;
    stimulus_sig(find(stimulus_sig)) = 1;
    
    if iCondition == 3

        rho(find(stimulus_z)) = 0; 
        rho_sig(find(stimulus_z)) = 0;
    
    end
    
    if iCondition == 4 
        
        rho(find(attention_z)) = 0; 
        rho_sig(find(attention_z)) = 0;
    
    end
    
%     rho(find(~(cluster_has_seeds)),:) = [];
%     rho(:,find(~(cluster_has_seeds))) = [];

%     rho = rho(idx,:);
%     rho = rho(:,idx);
    
    rho = rho(which_cluster_each_seed,:);
    rho = rho(:,which_cluster_each_seed);

    rho_sig = rho_sig(which_cluster_each_seed,:);
    rho_sig = rho_sig(:,which_cluster_each_seed);
    
    f = figure;
    
    imagesc(rho);
    hold on
    
    [l,c] = size(rho_sig);
    for il=1:l
        for ic=1:c
            if rho_sig(il,ic)
                plot(ic,il,'k*','MarkerSize',3);
            end
        end
    end
    
    title(label);

    clrmp = colormap('jet');
    clrmp(65,:) = clrmp(64,:);
    clrmp(33,:) = [1 1 1];

    min_C = -3;
    max_C = 3;

    %colorbar;
    caxis([min_C max_C]);
    %colorbar('Ticks',[min_C max_C/2 max_C]);
    colormap(clrmp);
    colorbar;
    
    hold off;

    ax = gca;
%    set(ax,'XTick',1:length(labels));
%    set(ax,'YTick',1:length(labels));
    set(ax,'XTick',1:length(new_labels));
    set(ax,'YTick',1:length(new_labels));
%    set(ax,'XTickLabel',labels);
%    set(ax,'YTickLabel',labels);
    set(ax,'XTickLabel',new_labels);
    set(ax,'YTickLabel',new_labels);
    set(ax,'FontSize',fs);

    xticklabel_rotate([],90,[],'Fontsize',fs);

    print(f,'-depsc',strcat('FC_Voxel_AAL_ROI_kmeans_',Analysis_Label,'_Clusters','-','Mean-Contrast-',label,'-Seeds-Sorted-Per-Net-Repeating-Corbetta.eps'));

end

end

function plotFCAndGrangerAtGroupLevelMeanOnlySeedsSortedPerNetwork

nRun = 32;
nTotalClusters = 758;
nClusters = 758;
nTR = 150;
pcriterion = 0.01;
nROI = 90;
fs = 5;

load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

% do = 'FC';
do = 'GC';

if strcmp(do,'FC') 

    load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters.mat');

    mean_track = zeros(nTotalClusters,nTotalClusters);
    mean_passive = zeros(nTotalClusters,nTotalClusters);
    mean_rest = zeros(nTotalClusters,nTotalClusters);

    for iRun=1:nRun

        mean_track = mean_track + FC_clusters(iRun).track_rho;
        mean_passive = mean_passive + FC_clusters(iRun).passive_rho;
        mean_rest = mean_rest + FC_clusters(iRun).rest_rho;

    end

    mean_track = mean_track ./ nRun;
    mean_passive = mean_passive ./ nRun;
    mean_rest = mean_rest ./ nRun;

    mean_track_pval = erfc((abs(mean_track).*sqrt(nTR))./sqrt(2));
    mean_passive_pval = erfc((abs(mean_passive).*sqrt(nTR))./sqrt(2));
    mean_rest_pval = erfc((abs(mean_rest).*sqrt(nTR))./sqrt(2));

    mean_track(mean_track_pval > pcriterion) = 0;
    mean_passive(mean_passive_pval > pcriterion) = 0;
    mean_rest(mean_rest_pval > pcriterion) = 0;

    track_mean = mean_track; 
    passive_mean = mean_passive; 
    rest_mean = mean_rest;
    
    Analysis_Label = 'FC'; 

end

if strcmp(do,'GC')

    load(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast','.mat'));

    track_mean = track_mean; 
    passive_mean = passive_mean;
    rest_mean = rest_mean;
    
    Analysis_Label = 'Granger'; 

end

for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

cluster_has_seeds = zeros(1,nClusters);
idx_cluster = 0;
seeds_labels = cell.empty;
iiCluster = 0;
iiseed = 0;
for iROI=1:nROI
        
    for iCluster=1:length(ROI(iROI).clusters)

        idx_cluster = idx_cluster + 1;

        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;

        Voxels = zeros(length(idx_voxels),3);

        for iVoxel=1:length(idx_voxels)

            [Voxels(iVoxel,1) Voxels(iVoxel,2) Voxels(iVoxel,3)] = ind2sub([91 109 91],idx_voxels(iVoxel));

        end

        seeds = seeIfAVoxelIsInsideASeed(Voxels);
        
        ROI_info{iROI,6} = seeds;
        
        if ~isempty(seeds)
            
            cluster_has_seeds(idx_cluster) = 1; 
        
            if iscell(seeds)
               
                iiseed = iiseed + 1;
                
                this_seed = seeds{1};
                all_seeds{iiseed} = this_seed;
                which_cluster_each_seed(iiseed) = idx_cluster;
                for i=2:length(seeds)
                    
                    iiseed = iiseed + 1;
                    this_seed = strcat(this_seed,' / ',seeds{i});
                    all_seeds{iiseed} = seeds{i};
                    which_cluster_each_seed(iiseed) = idx_cluster;
                    
                end
                
            else
                
                iiseed = iiseed + 1;
                this_seed = seeds;
                all_seeds{iiseed} = this_seed;
                which_cluster_each_seed(iiseed) = idx_cluster;
                
            end
            
            iiCluster = iiCluster + 1;
            
            %check_seeds{iiCluster,1} = seeds;
            
            seeds_labels{iiCluster} = this_seed;
            
        end

    end
        
end

[labels,idx] = sort(seeds_labels);
[new_labels,new_idx] = sort(all_seeds);
which_cluster_each_seed = which_cluster_each_seed(new_idx);

all_networks_labels = {'DAN' 'VAN' 'SMN' 'VIS' 'FPC' 'LAN' 'DMN' 'AUD'};

corbetta_networks_labels = {'DAN' 'VAN' 'VIS'};

corbetta_seeds = {'FEF' 'PIPS' 'SPL' 'r-P' 'PRECU' 'TPJ' 'V3' 'V4' 'V7' 'V8'};

idx_out = [];
for iSeed=1:length(which_cluster_each_seed)
   
    seed_name = new_labels{iSeed};
    
    this_seed_is_out = 1;
    
    for iNet=1:length(corbetta_networks_labels)
       
        if strcmp(seed_name(1:3),corbetta_networks_labels{iNet})
            
            for iMatch=1:length(corbetta_seeds)
            
                if strfind(lower(seed_name),lower(corbetta_seeds{iMatch}))
                    
                    this_seed_is_out = 0;
            
                end
                
            end
            
        end
        
    end
    
    if this_seed_is_out
        
        idx_out = [idx_out, iSeed];
        
    end
    
end

which_cluster_each_seed(idx_out) = [];
new_labels(idx_out) = [];

% zcriterion = 1.6;
zcriterion = 2.3;

for iCondition=1:3
    
    if iCondition == 1; rho = track_mean; label = 'Track'; end
    if iCondition == 2; rho = passive_mean; label = 'Passive'; end
    if iCondition == 3; rho = rest_mean; label = 'Rest'; end
       
    rho = rho(which_cluster_each_seed,:);
    rho = rho(:,which_cluster_each_seed);

    f = figure;
    
    rho(isnan(rho)) = 0;
    
    imagesc(rho);
    hold on
   
    title(label);

    clrmp = colormap('jet');
%     clrmp(65,:) = clrmp(64,:);
%     clrmp(33,:) = [1 1 1];

    min_C = 0;
    max_C = max(rho(:));

    %colorbar;
    caxis([min_C max_C]);
    %colorbar('Ticks',[min_C max_C/2 max_C]);
    colormap(clrmp);
    colorbar;
    
    hold off;

    ax = gca;
%    set(ax,'XTick',1:length(labels));
%    set(ax,'YTick',1:length(labels));
    set(ax,'XTick',1:length(new_labels));
    set(ax,'YTick',1:length(new_labels));
%    set(ax,'XTickLabel',labels);
%    set(ax,'YTickLabel',labels);
    set(ax,'XTickLabel',new_labels);
    set(ax,'YTickLabel',new_labels);
    set(ax,'FontSize',fs);

    xticklabel_rotate([],90,[],'Fontsize',fs);

    print(f,'-depsc',strcat('FC_Voxel_AAL_ROI_kmeans_',Analysis_Label,'_Clusters','-','Mean-',label,'-Seeds-Sorted-Per-Net-Repeating-Corbetta.eps'));

end

end

function whoIsTheStrongestConnection

which_cluster = [264 266 511 573 358];

load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-Contrast.mat');
load(strcat('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters','-','Mean-Contrast','.mat'));
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

% load DTI
prefix = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION';
sufix = 'preprocessed\B0-DTI\6.forbedpost.bedpostX\Und_PRE_xx.mat';

DTI(1) = load(strcat(prefix,'\','SUBJECT-1-22-10-2015','\',sufix));
DTI(2) = load(strcat(prefix,'\','SUBJECT-2-26-10-2015','\',sufix));
DTI(3) = load(strcat(prefix,'\','SUBJECT-3-3-11-2015','\',sufix));
DTI(4) = load(strcat(prefix,'\','SUBJECT-4-2-11-2015','\',sufix));
DTI(5) = load(strcat(prefix,'\','SUBJECT-5-2-11-2015','\',sufix));
DTI(6) = load(strcat(prefix,'\','SUBJECT-6-24-11-2015','\',sufix));
DTI(7) = load(strcat(prefix,'\','SUBJECT-7-14-01-2016','\',sufix));
DTI(8) = load(strcat(prefix,'\','SUBJECT-8-14-01-2016','\',sufix));

% compute Common DTI
ave_DTI = zeros(size(DTI(1).C));
for iDTI=1:8
    
    ave_DTI = ave_DTI + DTI(iDTI).C;
    
end
ave_DTI = ave_DTI ./ 8;

common_DTI = ave_DTI;
common_DTI(~(DTI(1).C & DTI(2).C & DTI(3).C & DTI(4).C & DTI(5).C & DTI(6).C & DTI(7).C & DTI(8).C)) = 0;

nClusters = 758;

% load Contrast - Attention

FC_contrast_Attention = contrast_attention_z;
GC_contrast_Attention = Attention_Contrast.Z;

zcriterion = 2.3;

FC_contrast_Attention(isnan(FC_contrast_Attention)) = 0;
GC_contrast_Attention(isnan(GC_contrast_Attention)) = 0;

FC_contrast_Attention(find(FC_contrast_Attention>(-1)*zcriterion & FC_contrast_Attention<zcriterion)) = 0; 
GC_contrast_Attention(find(GC_contrast_Attention>(-1)*zcriterion & GC_contrast_Attention<zcriterion)) = 0; 

% load Contrast - Stimulus

FC_contrast_Stimulus = contrast_stimulus_z;
GC_contrast_Stimulus = Stimulus_Contrast.Z;

zcriterion = 2.3;

FC_contrast_Stimulus(isnan(FC_contrast_Stimulus)) = 0;
GC_contrast_Stimulus(isnan(GC_contrast_Stimulus)) = 0;

FC_contrast_Stimulus(find(FC_contrast_Stimulus>(-1)*zcriterion & FC_contrast_Stimulus<zcriterion)) = 0; 
GC_contrast_Stimulus(find(GC_contrast_Stimulus>(-1)*zcriterion & GC_contrast_Stimulus<zcriterion)) = 0; 

% Attention and Stimulus - ONLY

FC_contrast_Attention_tmp = FC_contrast_Attention;
FC_contrast_Stimulus_tmp = FC_contrast_Stimulus;

FC_contrast_Attention(find(FC_contrast_Stimulus_tmp)) = 0;
FC_contrast_Stimulus(find(FC_contrast_Attention_tmp)) = 0;

% Direct and Indirect - Attention

FC_contrast_Attention_direct = FC_contrast_Attention;
FC_contrast_Attention_indirect = FC_contrast_Attention;

GC_contrast_Attention_direct = GC_contrast_Attention;
GC_contrast_Attention_indirect = GC_contrast_Attention;

FC_contrast_Attention_direct(~(common_DTI)) = 0;
GC_contrast_Attention_direct(~(common_DTI)) = 0;

FC_contrast_Attention_indirect(find(common_DTI)) = 0;
GC_contrast_Attention_indirect(find(common_DTI)) = 0;

% Direct and Indirect - Stimulus

FC_contrast_Stimulus_direct = FC_contrast_Stimulus;
FC_contrast_Stimulus_indirect = FC_contrast_Stimulus;

GC_contrast_Stimulus_direct = GC_contrast_Stimulus;
GC_contrast_Stimulus_indirect = GC_contrast_Stimulus;

FC_contrast_Stimulus_direct(~(common_DTI)) = 0;
GC_contrast_Stimulus_direct(~(common_DTI)) = 0;

FC_contrast_Stimulus_indirect(find(common_DTI)) = 0;
GC_contrast_Stimulus_indirect(find(common_DTI)) = 0;


for iCluster=which_cluster

    conn(1,:) = FC_contrast_Attention_direct(iCluster,:);
    conn(2,:) = GC_contrast_Attention_direct(iCluster,:);
    conn(3,:) = GC_contrast_Attention_direct(:,iCluster);
    
    tmp = max(conn,[],2);
    [tmp, i] = max(tmp);
    tmp = conn(i,:);
    [tmp_sorted,connected_clusters] = sort(tmp,'descend');
    
    if i == 1; label = 'FC'; end
    if i == 2; label = 'GC-FROM'; end
    if i == 3; label = 'GC-TO'; end
    
    disp(strcat(int2str(iCluster),'-',label,'-',int2str(connected_clusters(1))));
    
end

end

function plotSignOfChangeAndCommonOrigin

%%% LOAD CONTRASTS

% load('FC-Voxels-AAL-ROI-corr-KMeans-MeaningFul-FC-GC-Att-Stim-Only.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-Contrast.mat');
load('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters-Mean-Contrast.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Functional-Parcels.mat');

%%% DEFINE PARAMETERS

nClusters = 758;
nROI = 90;
zcriterion = 2.3;

%%% PROCESS MATRICES

FC_Contrast_Attention = contrast_attention_z;
GC_Contrast_Attention = Attention_Contrast.Z;
FC_Contrast_Stimulus = contrast_stimulus_z;
GC_Contrast_Stimulus = Stimulus_Contrast.Z;

FC_Contrast_Attention(isnan(FC_Contrast_Attention)) = 0;
GC_Contrast_Attention(isnan(GC_Contrast_Attention)) = 0;
FC_Contrast_Attention(find(FC_Contrast_Attention>(-1)*zcriterion & FC_Contrast_Attention<zcriterion)) = 0; 
GC_Contrast_Attention(find(GC_Contrast_Attention>(-1)*zcriterion & GC_Contrast_Attention<zcriterion)) = 0; 

FC_Contrast_Stimulus(isnan(FC_Contrast_Stimulus)) = 0;
GC_Contrast_Stimulus(isnan(GC_Contrast_Stimulus)) = 0;
FC_Contrast_Stimulus(find(FC_Contrast_Stimulus>(-1)*zcriterion & FC_Contrast_Stimulus<zcriterion)) = 0; 
GC_Contrast_Stimulus(find(GC_Contrast_Stimulus>(-1)*zcriterion & GC_Contrast_Stimulus<zcriterion)) = 0; 

FC_Attention_Only = FC_Contrast_Attention;
FC_Attention_Only(find(FC_Contrast_Stimulus)) = 0;
FC_Stimulus_Only = FC_Contrast_Stimulus;
FC_Stimulus_Only(find(FC_Contrast_Attention)) = 0;
FC_Common = FC_Contrast_Attention & FC_Contrast_Stimulus;

GC_Attention_Only = GC_Contrast_Attention;
GC_Attention_Only(find(GC_Contrast_Stimulus)) = 0;
GC_Stimulus_Only = GC_Contrast_Stimulus;
GC_Stimulus_Only(find(GC_Contrast_Attention)) = 0;
GC_Common = GC_Contrast_Attention & GC_Contrast_Stimulus;

%%% CALCULATE DENSITIES

for iCluster=1:nClusters
    
   vector_FC_A = FC_Contrast_Attention(iCluster,:);
   vector_FC_S = FC_Contrast_Stimulus(iCluster,:);
   vector_FC_A_Only = FC_Attention_Only(iCluster,:);
   vector_FC_S_Only = FC_Stimulus_Only(iCluster,:);
   vector_FC_C = FC_Common(iCluster,:);
   
   n_FC_A_p = length(find(vector_FC_A > 0));
   n_FC_A_n = length(find(vector_FC_A < 0));
   n_FC_S_p = length(find(vector_FC_S > 0));
   n_FC_S_n = length(find(vector_FC_S < 0));
   n_FC_A_o = length(find(vector_FC_A_Only));
   n_FC_S_o = length(find(vector_FC_S_Only));
   n_FC_C_p = length(find(vector_FC_C > 0));
   
   cluster_FC(1,iCluster) = n_FC_A_p;
   cluster_FC(2,iCluster) = n_FC_A_n;
   cluster_FC(3,iCluster) = n_FC_S_p;
   cluster_FC(4,iCluster) = n_FC_S_n;
   cluster_FC(5,iCluster) = n_FC_A_o;
   cluster_FC(6,iCluster) = n_FC_S_o;
   cluster_FC(7,iCluster) = n_FC_C_p;
   
end

for iCluster=1:nClusters
    
   vector_GC_A_FROM = GC_Contrast_Attention(iCluster,:);
   vector_GC_S_FROM = GC_Contrast_Stimulus(iCluster,:);
   vector_GC_A_Only_FROM = GC_Attention_Only(iCluster,:);
   vector_GC_S_Only_FROM = GC_Stimulus_Only(iCluster,:);
   vector_GC_C_FROM = GC_Common(iCluster,:);
   
   vector_GC_A_TO = GC_Contrast_Attention(:,iCluster);
   vector_GC_S_TO = GC_Contrast_Stimulus(:,iCluster);
   vector_GC_A_Only_TO = GC_Attention_Only(:,iCluster);
   vector_GC_S_Only_TO = GC_Stimulus_Only(:,iCluster);
   vector_GC_C_TO = GC_Common(:,iCluster);
   
   n_GC_A_p = length(find(vector_GC_A_FROM > 0)) + length(find(vector_GC_A_TO > 0));
   n_GC_A_n = length(find(vector_GC_A_FROM < 0)) + length(find(vector_GC_A_TO < 0));
   n_GC_S_p = length(find(vector_GC_S_FROM > 0)) + length(find(vector_GC_S_TO > 0));
   n_GC_S_n = length(find(vector_GC_S_FROM < 0)) + length(find(vector_GC_S_TO < 0));
   n_GC_A_o = length(find(vector_GC_A_Only_FROM)) + length(find(vector_GC_A_Only_TO));
   n_GC_S_o = length(find(vector_GC_S_Only_FROM)) + length(find(vector_GC_S_Only_TO));
   n_GC_C_p = length(find(vector_GC_C_FROM > 0)) + length(find(vector_GC_C_TO > 0));
   
   cluster_GC(1,iCluster) = n_GC_A_p;
   cluster_GC(2,iCluster) = n_GC_A_n;
   cluster_GC(3,iCluster) = n_GC_S_p;
   cluster_GC(4,iCluster) = n_GC_S_n;
   cluster_GC(5,iCluster) = n_GC_A_o;
   cluster_GC(6,iCluster) = n_GC_S_o;
   cluster_GC(7,iCluster) = n_GC_C_p;
   
end

cluster_C_label{1} = 'Increase with Attention';
cluster_C_label{2} = 'Decrease with Attention';
cluster_C_label{3} = 'Increase with Stimulus';
cluster_C_label{4} = 'Decrease with Stimulus';
cluster_C_label{5} = 'Unique to Attention';
cluster_C_label{6} = 'Unique to Stimulus';
cluster_C_label{7} = 'Common to Both';

max_val_FC = max(cluster_FC(:));
min_val_FC = min(cluster_FC(:));

max_val_GC = max(cluster_GC(:));
min_val_GC = min(cluster_GC(:));

plotDistributionDensity(min_val_FC,max_val_FC,cluster_FC(1,:),cluster_FC(2,:),cluster_C_label{1},cluster_C_label{2},'Attention Changes (FC)');

plotDistributionDensity(min_val_FC,max_val_FC,cluster_FC(3,:),cluster_FC(4,:),cluster_C_label{3},cluster_C_label{4},'Stimulus Changes (FC)');

plotDistributionDensity(min_val_FC,max_val_FC,cluster_FC(5,:),cluster_FC(7,:),cluster_C_label{5},cluster_C_label{7},'Attention Commonality (FC)');

plotDistributionDensity(min_val_FC,max_val_FC,cluster_FC(6,:),cluster_FC(7,:),cluster_C_label{6},cluster_C_label{7},'Stimulus Commonality (FC)');

plotDistributionDensity(min_val_GC,max_val_GC,cluster_GC(1,:),cluster_GC(2,:),cluster_C_label{1},cluster_C_label{2},'Attention Changes (GC)');

plotDistributionDensity(min_val_GC,max_val_GC,cluster_GC(3,:),cluster_GC(4,:),cluster_C_label{3},cluster_C_label{4},'Stimulus Changes (GC)');

plotDistributionDensity(min_val_GC,max_val_GC,cluster_GC(5,:),cluster_GC(7,:),cluster_C_label{5},cluster_C_label{7},'Attention Commonality (GC)');

plotDistributionDensity(min_val_GC,max_val_GC,cluster_GC(6,:),cluster_GC(7,:),cluster_C_label{6},cluster_C_label{7},'Stimulus Commonality (GC)');


%%% PLOT FUNCTIONAL CONNECTIVITY

% map = brewermap(2,'Set1'); 
% f = figure;
% vector_one = cluster_FC(1,:);
% vector_two = cluster_FC(2,:);
% vector_one(vector_one==0) = [];
% vector_two(vector_two==0) = [];
% max_plot = max([max(vector_one(:)) max(vector_two(:))]);
% min_plot = max([min(vector_one(:)) min(vector_two(:))]);
% histf(vector_one(:),min_plot:.01:max_plot,'facecolor',map(1,:),'facealpha',.5,'edgecolor','none')
% hold on
% histf(vector_two(:),min_plot:.01:max_plot,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none')
% box off
% axis tight
% legalpha(cluster_C_label{1},cluster_C_label{2},'location','northwest')
% legend boxoff

% map = brewermap(2,'Set1'); 
% f = figure;
% max_plot = max([max(cluster_FC(3,:)) max(cluster_FC(4,:))]);
% min_plot = max([min(cluster_FC(3,:)) min(cluster_FC(4,:))]);
% histf(cluster_FC(3,:),min_plot:.01:max_plot,'facecolor',map(1,:),'facealpha',.5,'edgecolor','none')
% hold on
% histf(cluster_FC(4,:),min_plot:.01:max_plot,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none')
% box off
% axis tight
% legalpha(cluster_C_label{3},cluster_C_label{4},'location','northwest')
% legend boxoff
% 
% map = brewermap(2,'Set1'); 
% f = figure;
% max_plot = max([max(cluster_FC(5,:)) max(cluster_FC(7,:))]);
% min_plot = max([min(cluster_FC(5,:)) min(cluster_FC(7,:))]);
% histf(cluster_FC(5,:),min_plot:.01:max_plot,'facecolor',map(1,:),'facealpha',.5,'edgecolor','none')
% hold on
% histf(cluster_FC(7,:),min_plot:.01:max_plot,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none')
% box off
% axis tight
% legalpha(cluster_C_label{5},cluster_C_label{5},'location','northwest')
% legend boxoff
% 
% map = brewermap(2,'Set1'); 
% f = figure;
% max_plot = max([max(cluster_FC(6,:)) max(cluster_FC(7,:))]);
% min_plot = max([min(cluster_FC(6,:)) min(cluster_FC(7,:))]);
% histf(cluster_FC(6,:),min_plot:.01:max_plot,'facecolor',map(1,:),'facealpha',.5,'edgecolor','none')
% hold on
% histf(cluster_FC(7,:),min_plot:.01:max_plot,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none')
% box off
% axis tight
% legalpha(cluster_C_label{6},cluster_C_label{7},'location','northwest')
% legend boxoff
% 
% %%% PLOT GRANGER CAUSALITY
% 
% map = brewermap(2,'Set1'); 
% f = figure;
% max_plot = max([max(cluster_GC(1,:)) max(cluster_GC(2,:))]);
% min_plot = max([min(cluster_GC(1,:)) min(cluster_GC(2,:))]);
% histf(cluster_FC(1,:),min_plot:.01:max_plot,'facecolor',map(1,:),'facealpha',.5,'edgecolor','none')
% hold on
% histf(cluster_FC(2,:),min_plot:.01:max_plot,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none')
% box off
% axis tight
% legalpha(cluster_C_label{1},cluster_C_label{2},'location','northwest')
% legend boxoff
% 
% map = brewermap(2,'Set1'); 
% f = figure;
% max_plot = max([max(cluster_GC(3,:)) max(cluster_GC(4,:))]);
% min_plot = max([min(cluster_GC(3,:)) min(cluster_GC(4,:))]);
% histf(cluster_GC(3,:),min_plot:.01:max_plot,'facecolor',map(1,:),'facealpha',.5,'edgecolor','none')
% hold on
% histf(cluster_GC(4,:),min_plot:.01:max_plot,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none')
% box off
% axis tight
% legalpha(cluster_C_label{3},cluster_C_label{4},'location','northwest')
% legend boxoff
% 
% map = brewermap(2,'Set1'); 
% f = figure;
% max_plot = max([max(cluster_GC(5,:)) max(cluster_GC(7,:))]);
% min_plot = max([min(cluster_GC(5,:)) min(cluster_GC(7,:))]);
% histf(cluster_GC(5,:),min_plot:.01:max_plot,'facecolor',map(1,:),'facealpha',.5,'edgecolor','none')
% hold on
% histf(cluster_GC(7,:),min_plot:.01:max_plot,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none')
% box off
% axis tight
% legalpha(cluster_C_label{5},cluster_C_label{7},'location','northwest')
% legend boxoff
% 
% map = brewermap(2,'Set1'); 
% f = figure;
% max_plot = max([max(cluster_GC(6,:)) max(cluster_GC(7,:))]);
% min_plot = max([min(cluster_GC(6,:)) min(cluster_GC(7,:))]);
% histf(cluster_FC(6,:),min_plot:.01:max_plot,'facecolor',map(1,:),'facealpha',.5,'edgecolor','none')
% hold on
% histf(cluster_FC(7,:),min_plot:.01:max_plot,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none')
% box off
% axis tight
% legalpha(cluster_C_label{6},cluster_C_label{7},'location','northwest')
% legend boxoff

end

function plotDistributionDensity(min_val,max_val,vector_one,vector_two,label_one,label_two,title_label)

f = figure;

nPoints = 1000;

vector_one(vector_one==0) = [];
vector_two(vector_two==0) = [];
XI= linspace(min_val,max_val,nPoints);

kvone = ksdensity(vector_one,XI);
kvtwo = ksdensity(vector_two,XI);
kvone = kvone .* (sum(vector_one(:)));
kvtwo = kvtwo .* (sum(vector_two(:)));

idx = find(kvone==0);
kvone(idx) = [];
XIone = XI;
XIone(idx) = [];
plot(XIone,kvone,'r');

hold on

idx = find(kvtwo==0);
kvtwo(idx) = [];
XItwo = XI;
XItwo(idx) = [];
plot(XItwo,kvtwo,'b');

xlabel('# of connections');
ylabel('# of clusters');

min_y = min([kvone(:); kvtwo(:)]) - 50;
max_y = max([kvone(:); kvtwo(:)]) + 50;

box off
axis tight
legalpha(strcat(label_one,'-','(Total # of Connections:',int2str(sum(vector_one(:))),')'),strcat(label_two,'-','(Total # of Connections:',int2str(sum(vector_two(:))),')'),'location','northeast')
legend boxoff

title(title_label);

xlim([(min_val - 5) (max_val + 5)]);
ylim([min_y max_y]);

print(f,'-depsc',strcat('Cluster-Distribution-',strrep(title_label,' ','-'),'.eps'));

end


