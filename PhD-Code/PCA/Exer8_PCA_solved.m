function Exer8_PCA_solved

clear all;
close all;

%% load Christoph's files

load('MUA_b_t_g');

% Nb number of bursts
% Nt number of time points
% ti negative latencies before peak activity [ms]
% Ng number of neuron groups (20 neurons each)

who


%% compute and plot mean timecourse, averaging over bursts

mu_MUA_g_t    = nan( Ng, Nt );

for ig = 1 : Ng
    mu_MUA_g_t( ig, : ) = mean( squeeze( MUA_b_t_g( :, :, ig ) ), 1 ); % mean over bursts
end

figure;
hold on;
for ig = 1 : Ng
    
    cig = getcolor( ig );
    plot( ti, mu_MUA_g_t( ig, :), cig, 'LineWidth', 2 );
    
end;
hold off;
axis([min(ti) max(ti) -50 max(mu_MUA_g_t(:))]);

title( 'average time course of group activities' );

%% compute and plot covariance matrix

X = mu_MUA_g_t;                     % Ng variables (rows), Nt observations (columns)

[m,n] = size( X )                  % m variables (rows), n observations (columns)
mn    = mean( X, 2);                % compute row mean
X0    = X - repmat(mn,1,n);         % subtract row mean

CX = X0 * X0' / (n-1);

figure;
hold on;
pcolor( CX );
caxis([min(CX(:)) max(CX(:))]);
colormap(hot);
hold off;
axis([1 m 1 m]);
axis 'square';
xlabel( 'group' );
ylabel( 'group' );
title(' covariance of group activity' );


%% identify principal components in mean timecourses





[U, S, V] = svd( X0' );             % skinny SVD

size(U)
size(S)
size(V)

figure;
hold on;
bar( log( diag(S) ) );              % plot variance of principal components
hold off;
xlabel( 'PCs');
ylabel( 'log variance' );

title( 'log variance explained by PCs' );

P   = V';                           % principal components as row vectors (for projection)


Y0 = P*X0;                          % project zero-mean timecourses onto principal components
size( Y0 )

figure;                             % illustrate first 11
subplot(2,2,1);
hold on;
plot3( Y0(1,:), Y0(2,:), Y0(3,:), 'r', 'LineWidth', 2 );    % plot time-development in PC space
plot3( Y0(1,1), Y0(2,1), Y0(3,1), 'ro', 'LineWidth', 2 );
plot_3_axes( gca );
hold off;
axis 'equal';
xlabel('PC 1');
ylabel('PC 2');
zlabel('PC 3');
view( -37.5, 30 );

subplot(2,2,2);
hold on;
plot3( Y0(4,:), Y0(5,:), Y0(6,:), 'b', 'LineWidth', 2 );    % plot time-development in PC space
plot3( Y0(4,1), Y0(5,1), Y0(6,1), 'bo', 'LineWidth', 2 );
plot_3_axes( gca );
hold off;
axis 'equal';
view( -37.5, 30 );
xlabel('PC 4');
ylabel('PC 5');
zlabel('PC 6');

subplot(2,2,3);
hold on;
plot3( Y0(7,:), Y0(8,:), Y0(9,:), 'k', 'LineWidth', 2 );    % plot time-development in PC space
plot3( Y0(7,1), Y0(8,1), Y0(9,1), 'ko', 'LineWidth', 2 );
plot_3_axes( gca );
hold off;

axis 'equal';
view( -37.5, 30 );
xlabel('PC 7');
ylabel('PC 8');
zlabel('PC 9');

subplot(2,2,4);
hold on;
plot( Y0(10,:), Y0(11,:), 'g', 'LineWidth', 2 );    % plot time-development in PC space
plot( Y0(10,1), Y0(11,1), 'go', 'LineWidth', 2 );
plot_2_axes( gca );
hold off;
axis 'equal';
xlabel('PC 10');
ylabel('PC 11');

title( 'average time course, in space of first 11 PCs' );

%% compute and visualize covariance of transformed activity

CY = Y0 * Y0' / (n-1);

figure;
hold on;
pcolor( CY );
caxis([min(CY(:)) max(CY(:))]);
colormap(hot);
hold off;
axis([1 m 1 m]);
axis 'square';
xlabel( 'group' );
ylabel( 'group' );
title(' covariance of group activity' );


%% de-noise mean timecourse by backprojecting first three principal components

Y0dn = Y0;
Y0dn( 4:Ng, 1:Nt ) = zeros(1:Ng-3, 1:Nt);                  % zero later PCs

X0dn = P' * Y0dn;                                          % project back to groups
Xdn  = X0dn + repmat(mn,1,n);                              % add back row mean

figure;
hold on;
for ig = 1 : Ng
    
    cig = getcolor( ig );
    plot( ti, Xdn( ig, :), cig, 'LineWidth', 2 );
    
end;
hold off;
axis([min(ti) max(ti) -50 max(mu_MUA_g_t(:))]);


title( 'average time course of groups, denoised with first 3 PCs' );


%% visualize timecourse of individual bursts in terms of principal components

MUA_g_t    = nan( Ng, Nt * Nb );

for ig = 1 : Ng
    MUA_g_t( ig, : ) = reshape( squeeze( MUA_b_t_g( :, :, ig ) )', 1, Nt * Nb ); % concatenate bursts
end


X = MUA_g_t;                                             % Ng groups (rows), Nt * Nb observations (columns)
[m,n] = size( X );                                       % m variables (rows), n observations (columns)
mn    = mean( X, 2);                                     % compute row mean
X0    = X - repmat(mn,1,n);                              % subtract row mean

Y0 = P*X0;                                                % project all observations projected on PCs
size( Y0 )

figure;

hold on;
for ib = 1:Nb
    
    ci = getcolor( ib );
    tidx = (ib-1)*Nt + [1:Nt];
    
    plot3( Y0(1,tidx), Y0(2,tidx), Y0(3,tidx), ci, 'LineWidth', 2 );    % plot time-development in PC space
    plot3( Y0(1,tidx(1)), Y0(2,tidx(1)), Y0(3,tidx(1)), [ci 'o'], 'LineWidth', 2, 'MarkerSize', 10 );
end
plot_3_axes( gca );
hold off;
axis 'equal';
xlabel('PC 1');
ylabel('PC 2');
zlabel('PC 3');
view( -37.5, 30 );



title( 'individual time course of bursts, in space of first 3 PCs' );


return;

%% some plotting functions

function plot_3_axes( gca )

XLim = get(gca, 'XLim'); YLim = get(gca, 'YLim');  ZLim = get(gca, 'ZLim');

plot3( XLim, [0 0], [0 0], 'k' );
plot3( [0 0], YLim, [0 0], 'k' );
plot3( [0 0], [0 0], [ZLim], 'k' );
axis( [XLim YLim ZLim ] );

return;

function plot_2_axes( gca )

XLim = get(gca, 'XLim'); YLim = get(gca, 'YLim'); 

plot3( XLim, [0 0], [0 0], 'k' );
plot3( [0 0], YLim, [0 0], 'k' );
axis( [XLim YLim ] );

return;

function ci = getcolor( idx )

clr = 'rgbcmyk';

cidx = mod( idx, length(clr) ) + 1;

ci = clr( cidx );

return;
    
%% save data for exercise

function SaveForExercise( Nburst, Nt, Ngroup, ti, MUA_b_t_g )

filename = 'MUA_b_t_g.mat';

save( filename, 'Nburst', 'Nt', 'Ngroup', 'ti', 'MUA_b_t_g' );

return;
    