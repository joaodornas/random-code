function real_FC_voxel_AAL_ROI_plot

settings_subj1_2210;
%settings_subj2_2610;
%settings_subj3_0311;
%settings_subj4_0211;
%settings_subj5_0211;
%settings_subj6_2411;

%%% #59-Parietal_Sup_L #25-Frontal_Med_Orb_L #5-Frontal_Sup_Orb_L
%%% #9-Frontal_Mid_Orb_L #15-Frontal_Inf_Orb_L #57-Postcentral_L
%%% #89-Temporal_Inf_L #47-Lingual_L #53-Occipital_Inf_L
%%% #65-Angular_L
%%% #19-SMA_L #73-Putamen_L #82-Temporal_Sup_R 
%%% #91-Cerebellum #92-Cerebellum #98-Cerebellum #99-Cerebellum
%%% idx_ROI = [59 25 5 9 15 57 89 47 53 65 19 73 82];
idx_ROI = [25 5 9 15 57 89 47 53 65 19 73 82];
plotConditionsRunsClusteredPerSession(settings,idx_ROI,1);

end

function plotConditionsRunsClusteredPerSession(settings,idx_ROI,shouldIPlotChanges)
     
%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRuns = 4;

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    TrackRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'.mat'));
    PassiveRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'.mat'));
    RestingStateRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'.mat'));
        
    %for irun=1:nRuns
        for irun=1:1
       
        disp(strcat('irun=',int2str(irun)));
        
        track_rho = TrackRHOs.FC_Voxels.run(irun).rho_Track;
        passive_rho = PassiveRHOs.FC_Voxels.run(irun).rho_Passive;
        rest_rho = RestingStateRHOs.FC_Voxels.run(irun).rho_RestingState;
        
        track_pval = TrackRHOs.FC_Voxels.run(irun).pval_Track;
        passive_pval = PassiveRHOs.FC_Voxels.run(irun).pval_Passive;
        rest_pval = RestingStateRHOs.FC_Voxels.run(irun).pval_RestingState;
        
        disp('cluster Resting State');
        
        [IdxClusters, Tidx, Ncluster] = ClusterWithKmeans( rest_rho, rest_pval );
        
        Ncomponent = size(rest_rho,1);
        cluster_assignment = IdxClusters;
    
%         plot_label = strcat(area_label,'-','Track','-','Run','-',int2str(irun),'-FC-Clustered');
%         plotFCClustered(settings, track_rho, cluster_assignment, Ncluster, Ncomponent, plot_label);
%         
%         plot_label = strcat(area_label,'-','Passive','-','Run','-',int2str(irun),'-FC-Clustered');
%         plotFCClustered(settings, passive_rho, cluster_assignment, Ncluster, Ncomponent, plot_label);
%         
%         plot_label = strcat(area_label,'-','Rest','-','Run','-',int2str(irun),'-FC-Clustered');
%         plotFCClustered(settings, rest_rho, cluster_assignment, Ncluster, Ncomponent, plot_label);
        
        if shouldIPlotChanges
            
           shouldIDoCluster = 0; 
           [a, b, c, d, e, f, track_increase, passive_increase, track_decrease, passive_decrease, g, h, i, j, k, l, m, n, o, p, q, r] = getFCDifferencePositiveNegativeFraction(track_rho,track_pval,passive_rho,passive_pval,shouldIDoCluster); 
        
           plot_label = strcat(area_label,'-','T-P','-','Track-Increase','-','Run','-',int2str(irun),'-FC-Clustered');
           plotFCClustered(settings, track_increase, cluster_assignment, Ncluster, Ncomponent, plot_label);
        
           plot_label = strcat(area_label,'-','T-P','-','Track-Decrease','-','Run','-',int2str(irun),'-FC-Clustered');
           plotFCClustered(settings, track_decrease, cluster_assignment, Ncluster, Ncomponent, plot_label);
        
           plot_label = strcat(area_label,'-','T-P','-','Passive-Increase','-','Run','-',int2str(irun),'-FC-Clustered');
           plotFCClustered(settings, passive_increase, cluster_assignment, Ncluster, Ncomponent, plot_label);
        
           plot_label = strcat(area_label,'-','T-P','-','Passive-Decrease','-','Run','-',int2str(irun),'-FC-Clustered');
           plotFCClustered(settings, passive_decrease, cluster_assignment, Ncluster, Ncomponent, plot_label);
        
           shouldIDoCluster = 0; 
           [a, b, c, d, e, f, passive_increase, rest_increase, passive_decrease, rest_decrease, g, h, i, j, k, l, m, n, o, p, q, r] = getFCDifferencePositiveNegativeFraction(passive_rho,passive_pval,rest_rho,rest_pval,shouldIDoCluster); 
        
           plot_label = strcat(area_label,'-','P-R','-','Passive-Increase','-','Run','-',int2str(irun),'-FC-Clustered');
           plotFCClustered(settings, passive_increase, cluster_assignment, Ncluster, Ncomponent, plot_label);
        
           plot_label = strcat(area_label,'-','P-R','-','Passive-Decrease','-','Run','-',int2str(irun),'-FC-Clustered');
           plotFCClustered(settings, passive_decrease, cluster_assignment, Ncluster, Ncomponent, plot_label);
        
           plot_label = strcat(area_label,'-','P-R','-','Rest-Increase','-','Run','-',int2str(irun),'-FC-Clustered');
           plotFCClustered(settings, rest_increase, cluster_assignment, Ncluster, Ncomponent, plot_label);
        
           plot_label = strcat(area_label,'-','P-R','-','Rest-Decrease','-','Run','-',int2str(irun),'-FC-Clustered');
           plotFCClustered(settings, rest_decrease, cluster_assignment, Ncluster, Ncomponent, plot_label);
        
        end
    
    end
    
end

end

function plotFCClustered(settings,rho_mat,cluster_assignment, Ncluster, Ncomponent, plot_label)

    source_order = 1:Ncomponent;
    [target_order, cluster_idx_sorted] = TargetOrder( cluster_assignment, Ncluster );

    [rho_sorted, ~] = SortMatrix( rho_mat, source_order, target_order, Ncomponent );
    
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
    title(plot_label);
    set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

    h = colorbar;
    set(h,'YLim',rlimit*[(-1+1/32) 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1]);

    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-',plot_label,'.jpg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-',plot_label,'.eps'));
    print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-',plot_label,'.pdf'));
    
    close all
    
end

function [IdxClusters, Tidx, Ncluster] = ClusterWithKmeans( rho_input, pval_input )

tic

nVoxels = size(rho_input,1);

%% keep only significant correlations
    
    pcriterion = 0.01;

    rho_input( pval_input > pcriterion ) = 0;
    
%% find NaNs

    kk = find(~any(rho_input,2));
    
    rho_input(kk,:) = [];
    rho_input(:,kk) = [];

%% compute the clusters

    num_Clust = 10;

    [Tidx,C,sumd,D] = kmeans(rho_input, num_Clust , 'distance', 'correlation', 'display', 'final','replicate',20);

    Ncluster = max(Tidx);

%% repositions voxels clusters with NaNs

IdxClusters = zeros(nVoxels,1);
ik = 0;

for iVoxel=1:nVoxels
    
    if find(kk==iVoxel)
        
        IdxClusters(iVoxel) = NaN;
        
    else
        
        ik = ik + 1;
        
        IdxClusters(iVoxel) = Tidx(ik);

    end
    
end

toc

end

%% funtions for sorting clusters of rows/columns by size
function [target_order, cluster_idx_sorted] = TargetOrder( cluster_assignment, ncluster )

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

%% funtions for sorting correlation matrices
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

%% remove nan rows/columns 
function [rho_output,pval_output] = removeNaNs(rho_input,pval_input)

[Ncomponent, ~] = size( rho_input );

source_order = 1 : Ncomponent;

kk = find( isnan( rho_input(1,:) ) );

['remove ' num2str(length(kk)) ' NaN rows/columns']

for i = length(kk) : -1 : 1
    
    idelete = kk(i);
    
    [temp_rr, ~]  = MSwap( rho_input, source_order, idelete, Ncomponent );
    [temp_pr, ~]  = MSwap( pval_input, source_order, idelete, Ncomponent );
    
    clear rho_input pval_input;
    
    Ncomponent = Ncomponent - 1;
    source_order = 1 : Ncomponent;
    
    rho_input  = temp_rr( 1:Ncomponent, 1: Ncomponent);
    pval_input = temp_pr( 1:Ncomponent, 1: Ncomponent);
    
end

rho_ouput = rho_input;
pval_output = pval_input;

end

% swap locations i and k, plus update order accordingly
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


