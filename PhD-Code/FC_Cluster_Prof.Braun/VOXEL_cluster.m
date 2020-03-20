function VOXEL_cluster

clear all;
close all;

global fs;

fs = 14;


%% load and clearn correlation matrix

% load('LOW-HIGH-ATTENTION-Jan-08-05-2015-FC-Voxels-per-ROI-AAL-Parietal-Sup.mat');

load('LOW-HIGH-ATTENTION-Jan-08-05-2015-FC-Voxels-per-ROI-AAL.mat');

iROI = 5;

who

ROI_label = char( label_ROI( iROI ) )

filenameH = [ROI_label '_hist'];

filenameFC = [ROI_label '_FC'];

filenameFCcond = [ROI_label '_FCcond'];


%% take a first look

rho_att = FC_Voxels.ROI(iROI).run.rho_MOT4;
pval_att = FC_Voxels.ROI(iROI).run.pval_MOT4;

rho_rest = FC_Voxels.ROI(iROI).run.rho_RestingState;
pval_rest = FC_Voxels.ROI(iROI).run.pval_RestingState;

[Ncomponent, ~] = size( rho_rest )

clear FC_Voxels;

View2Histograms( rho_att, rho_rest, pval_att, pval_rest );

suptitle( [ROI_label ' ' num2str(Ncomponent) ' voxels'] );

print( filenameH, '-depsc2') ;

%% remove nan rows/columns

source_order = 1 : Ncomponent;

kr = find( isnan( rho_rest(1,:) ) );

ka = find( isnan( rho_att(1,:) ) );

kk = unique([ kr ka ] );

['remove ' num2str(length(kk)) ' NaN rows/columns']

for i = length(kk) : -1 : 1
    
    idelete = kk(i);
    
    [temp_rr, ~]  = MSwap( rho_rest, source_order, idelete, Ncomponent );
    [temp_pr, ~]  = MSwap( pval_rest, source_order, idelete, Ncomponent );
    [temp_ra,  ~]  = MSwap( rho_att, source_order, idelete, Ncomponent );
    [temp_pa,  ~]  = MSwap( pval_att, source_order, idelete, Ncomponent );
    
    clear rho_rest pval_rest rho_att pval_att;
    
    Ncomponent = Ncomponent - 1;
    source_order = 1 : Ncomponent;
    
    rho_rest  = temp_rr( 1:Ncomponent, 1: Ncomponent);
    pval_rest = temp_pr( 1:Ncomponent, 1: Ncomponent);
    rho_att   = temp_ra( 1:Ncomponent, 1: Ncomponent);
    pval_att  = temp_pa( 1:Ncomponent, 1: Ncomponent);
    
end

%% keep only significant correlations
    
pcriterion = 0.01;

rho_att( pval_att > pcriterion ) = 0;

rho_rest( pval_rest > pcriterion ) = 0;

Avalue = rho_att;   % attention 

Rvalue = rho_rest;  % resting


%% Kmeans clustering of correlation matrix

[cluster_assignment, Ncluster] = ClusterWithKmeans( Rvalue, Ncomponent );  % kmeans cluster original correlation matrix (for resting)

source_order = 1:Ncomponent;
[target_order, cluster_idx_sorted] = TargetOrder( cluster_assignment, Ncluster );

[Rsorted, ~] = SortMatrix( Rvalue, source_order, target_order, Ncomponent );

[Asorted, ~] = SortMatrix( Avalue, source_order, target_order, Ncomponent );

View4Matrices( Avalue, Asorted, Rvalue, Rsorted, Ncomponent );

suptitle([ROI_label ' ' num2str(Ncomponent) ' voxels']);

print( filenameFC, '-depsc2') ;


%% condense FC matrix

%  cluster_idx_sorted lists the cluster numbers (chosen arbitrarily by kmeans) in descending order of size
%  
%  cluster_assignment assigns the nc rows/columns of the ORIGINAL matrix to
%  different clusters
%

[Rcondensed, Acondensed, Warea, cluster_assignment_sorted] = CondenseFC( Rvalue, Avalue, Ncomponent, cluster_assignment, cluster_idx_sorted, Ncluster );

% Rcondensed: Rvalue shrunk to Ncluster x Ncluster

% Acondensed: Avalue shrunk to Ncluster x Ncluster

% Warea:      matrix of weights (product of cluster sizes)

% cluster_assignment_sorted assigns the nc rows/columns of the ORIGINAL
% matrix to the sorted clusters (in the same order as in the condensed
% matrices

% cluster_size_sorted lists the number of voxels per cluster

for ic = 1 : Ncluster
     cluster_size_sorted( ic ) = length( find( cluster_assignment_sorted == ic ) ); 
end

cluster_size_sorted

View2Matrices( Rcondensed, Acondensed, Ncluster );

suptitle([ROI_label ' ' num2str(Ncluster) ' clusters']);

print( filenameFCcond, '-depsc2') ;

return;

%% condense voxel correlations to cluster correlations
%
%  for each cluster pair, compute mean of SIGNIFICANT voxel correlations
%

function [Rcondensed, Acondensed, Warea, cluster_assignment_sorted] = CondenseFC( Rvalue, Avalue, Ncomponent, cluster_assignment, cluster_idx_sorted, Ncluster )

cluster_assignment_sorted = nan(1,Ncomponent);

for ic = 1 : Ncluster    % loop over clusters in order of descending size
    kk = find( cluster_assignment == cluster_idx_sorted(ic) );   % find locations of rows/columns in original matrix
     
    cluster_size(ic) = length( kk );   % for debugging only
    
    cluster_assignment_sorted(kk) = ic;   % assign sorted cluster numbers to locationns of rows/columns
end

Rcondensed = nan( Ncluster, Ncluster );
Acondensed = nan( Ncluster, Ncluster );
Warea      = nan( Ncluster, Ncluster );

for ic = 1 : Ncluster   % fill in diagonal
    
    i_idx  = find( cluster_idx_sorted( ic ) == cluster_assignment );
    
    i_size = length( i_idx );
    
    R = Rvalue( i_idx, i_idx );
    kkR = find(R ~= 0 );
    
    A = Avalue( i_idx, i_idx );
    kkA = find(A ~= 0 );
    
    Rcondensed( ic, ic ) = mean( R( kkR ) );
    
    Acondensed( ic, ic ) = mean( A( kkA ) );
    
    Warea( ic, ic )      = i_size * i_size;
    
end

for ic = 1 : Ncluster-1   % fill in off-diagonal
    
    i_idx  = find( cluster_idx_sorted( ic ) == cluster_assignment );
    
    i_size = length( i_idx );
    
    for jc = ic+1 : Ncluster
        
        j_idx  = find( cluster_idx_sorted( jc ) == cluster_assignment );
        
        j_size = length( j_idx );
        
        R = Rvalue( i_idx, j_idx );
        kkR = find(R ~= 0 );
        
        A = Avalue( i_idx, j_idx );
        kkA = find(A ~= 0 );
        
        Rcondensed( ic, jc ) = mean( R( kkR ) );
        
        Acondensed( ic, jc ) = mean( A( kkA ) );
        
        Warea( ic, jc )      = i_size * j_size;
        
        Rcondensed( jc, ic ) = Rcondensed( ic, jc );
        
        Acondensed( jc, ic ) = Acondensed( ic, jc );
        
        Warea( jc, ic )      = Warea( ic, jc );
        
    end
    
end


return;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cluster with Kmeans
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Tidx, Ncluster] = ClusterWithKmeans( Rvalue, Ncomponent )

num_Clust = 10;

[Tidx,C,sumd,D] = kmeans(Rvalue, num_Clust , 'distance', 'correlation', 'display', 'final','replicate',20);

Ncluster = max(Tidx);

return;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% funtions for sorting correlation matrices
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

return;

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

return;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% funtions for sorting clusters of rows/columns by size
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [target_order, cluster_idx_sorted] = TargetOrder( cluster_assignment, ncluster )

% cluster_assignment assigns nc rows/columns to ncluster clusters

% target_order assigns nc matrix locations to nc rows/columns

% the first locations are filled with rows/columns from the largest cluster

% the next locations are filled with rows/columns from the next smaller
% cluster

% and so on

% cluster_idx_sorted lists the clusters in descending order of size ...

nc = length( cluster_assignment );


for kc = 1 : ncluster   % loop over clusters
    
    cluster_size(kc) = length( find( cluster_assignment == kc ) );
    
end

[cluster_size_sorted, cluster_idx_sorted] = sort( cluster_size, 'descend' );   % c_size_sorted = c_size( c_idx_sorted )

cluster_size_sorted

target_order  = [];  

for kc = 1 : ncluster % loop over sorted clusters
    
    target_order = [target_order; find( cluster_idx_sorted(kc) == cluster_assignment )];  % append indices of matrix rows/columns in current cluster
       
end

target_order = reshape( target_order, 1, nc );

return;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% funtions for visualization
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function View2Histograms( rho_att, rho_rest, pval_att, pval_rest )

xa = linspace(-1,1.5,101);
xb = linspace(0,1,101);
XLima = [-1 2];
XLimb = [0 1];
YLim = [0 3e5];

figure;
subplot(1,2,1);
hold on;
na = hist( rho_att(:), xa );
plot( xa, (na), 'r', 'LineWidth', 2 );
na = hist( rho_rest(:), xa );
plot( xa, (na), 'b', 'LineWidth', 2 );
hold off;
h = legend('attention', 'resting' );
set( gca, 'XLim', XLima, 'YLim', YLim );
axis 'square';


subplot(1,2,2);
hold on;
nb = hist( pval_att(:), xb );
plot( xb, (nb), 'r', 'LineWidth', 2 );
nb = hist( pval_rest(:), xb );
plot( xb, (nb), 'b', 'LineWidth', 2 );
hold off;
set( gca, 'XLim', XLimb, 'YLim', YLim );
axis 'square';

return;

function View2Matrices( Aval, Bval, Ncomponent )

rlimit = 0.5;

rmin = -rlimit + rlimit/32;
rmax =  rlimit;

Tick = [1 Ncomponent];
TickLabel = { '1' ; num2str(Ncomponent) };

A = zeros( Ncomponent+1, Ncomponent+1 );
A(1:Ncomponent, 1:Ncomponent) = Aval;

B = zeros( Ncomponent+1, Ncomponent+1 );
B(1:Ncomponent, 1:Ncomponent) = Bval;

figure;

subplot(1,2,1);
clrmp = colormap('jet');
clrmp(32,:) = [1 1 1];

hold on;
caxis([rmin rmax]);
h = pcolor(0.5:Ncomponent+0.5, 0.5:Ncomponent+0.5, A );
set( h, 'EdgeColor', 'none');
colormap( clrmp);
hold off;
axis 'square'; 
axis([0 Ncomponent+1 0 Ncomponent+1]);
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

h = colorbar;
set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1]);

subplot(1,2,2);
clrmp = colormap('jet');
clrmp(32,:) = [1 1 1];

hold on;
caxis([rmin rmax]);
h = pcolor( 0.5:Ncomponent+0.5, 0.5:Ncomponent+0.5, B );
set( h, 'EdgeColor', 'none');
colormap(clrmp);
hold off;
axis 'square'; 
axis([0 Ncomponent+1 0 Ncomponent+1]);
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

h = colorbar;
set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1]);
    
return;

%% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function View4Matrices( Aval, Bval, Cval, Dval, Ncomponent )

rlimit = 0.5;

rmin = -rlimit + rlimit/32;
rmax =  rlimit;

Tick = [1 Ncomponent];
TickLabel = { '1' ; num2str(Ncomponent) };

figure;

subplot(2,2,1);
clrmp = colormap('jet');
clrmp(32,:) = [1 1 1];

hold on;
caxis([rmin rmax]);
h = pcolor( 0.5:Ncomponent-0.5, 0.5:Ncomponent-0.5, Aval );
set( h, 'EdgeColor', 'none');
colormap( clrmp);
hold off;
axis 'square'; 
axis([0 Ncomponent 0 Ncomponent]);
xlabel('attention');
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

h = colorbar;
set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1]);

subplot(2,2,2);
clrmp = colormap('jet');
clrmp(32,:) = [1 1 1];

hold on;
caxis([rmin rmax]);
h = pcolor( 0.5:Ncomponent-0.5, 0.5:Ncomponent-0.5, Bval );
set( h, 'EdgeColor', 'none');
colormap(clrmp);
hold off;
axis 'square'; 
axis([0 Ncomponent 0 Ncomponent]);
xlabel('attention');
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

h = colorbar;
set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1]);

subplot(2,2,3);
clrmp = colormap('jet');
clrmp(32,:) = [1 1 1];

hold on;
caxis([rmin rmax]);
h = pcolor( 0.5:Ncomponent-0.5, 0.5:Ncomponent-0.5, Cval );
set( h, 'EdgeColor', 'none');
colormap( clrmp);
hold off;
axis 'square'; 
axis([0 Ncomponent 0 Ncomponent]);
xlabel('resting');
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );
h = colorbar;
set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1]);

subplot(2,2,4);
clrmp = colormap('jet');
clrmp(32,:) = [1 1 1];

hold on;
caxis([rmin rmax]);
h = pcolor( 0.5:Ncomponent-0.5, 0.5:Ncomponent-0.5, Dval );
set( h, 'EdgeColor', 'none');
colormap(clrmp);
hold off;
axis 'square'; 
axis([0 Ncomponent 0 Ncomponent]);
xlabel('resting');
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

h = colorbar;
set(h,'YLim',rlimit*[-1 1], 'YTick',rlimit*[-1 0 1], 'PlotBoxAspectRatio', [1 20 1]);
    
return;
