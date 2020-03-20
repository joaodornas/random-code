function PCA_of_MRI_voxels

clear all;
close all;

%% load Joao's files

% fs  = { 'trk4'; 'trk2'; 'rest' };
% rs  = {'1'; '2' );

lobe = 'occipital';

fs = {'rest'};
rs = {'1'};

cube = LoadCube( ['cube_' lobe '.mat'], fs, rs );

[Nx, Ny, Nz, Nt] = size( cube );


%% plot mean timecourses

Nv = Nx * Ny * Nz;

A_v_t    = reshape( cube, Nv, Nt );

[lobe ' ' char(fs) ' ' char(rs)]

size( A_v_t )

figure;
hold on;
for iv = 1 : Nv
    civ = getcolor( iv );
    plot( A_v_t( iv, :), civ, 'LineWidth', 2 );   
end;
hold off;
%axis([min(ti) max(ti) -50 max(mu_MUA_g_t(:))]);

title( 'voxel time courses' );




%% compute and plot covariance matrix

X = A_v_t;                     % Nv variables (rows), Nt observations (columns)

[m,n] = size( X );                  % m variables (rows), n observations (columns)
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
title(' covariance of voxel' );


%% identify principal components 


[U, S, V] = svd( X0' );             % skinny SVD

figure;
hold on;
bar( log( diag(S) ) );              % plot variance of principal components
hold off;
xlabel( 'PCs');
ylabel( 'log variance' );

title( 'log variance explained by PCs' );




P   = V';                           % principal components as row vectors (for projection)

Y0 = P*X0;                          % project zero-mean timecourses onto principal components

                                    
                                    
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
title(' covariance of transformed activity' );


%% plot mean timecourses of transformed activity

Nv = Nx * Ny * Nz;

Ncutoff = 3;

figure;
subplot(2,1,1);
hold on;
for iv = 1 : Ncutoff
    civ = getcolor( iv );
    plot( Y0( iv, :), civ, 'LineWidth', 2 );   
end;

aveX0 = mean( X0, 1 ); 
    plot( aveX0, 'k--', 'LineWidth', 2 );   
hold off;

subplot(2,1,2);
hold on;
for iv = Ncutoff+1 : Nv
    civ = getcolor( iv );
    plot( Y0( iv, :), civ, 'LineWidth', 2 );   
end;
hold off;
%axis([min(ti) max(ti) -50 max(mu_MUA_g_t(:))]);

suptitle( 'transformed time courses' );




%% de-noise mean timecourse by backprojecting selected principal components

Y0dn = Y0;
Y0dn( 1:Ncutoff, 1:Nt ) = zeros(Ncutoff, Nt);                  % zero later PCs

X0dn = P' * Y0dn;                                          % project back to groups
Xdn  = X0dn + repmat(mn,1,n);                              % add back row mean

figure;
hold on;
for iv = 1 : Nv
    
    civ = getcolor( iv );
    plot( Xdn( iv, :), civ, 'LineWidth', 2 );
    
end;
hold off;
%axis([min(ti) max(ti) -50 max(mu_MUA_g_t(:))]);


title( 'voxel time courses, denoised' );


%% interpret average of zero-mean timecourses

% Y0 (Nv x Nt ) time course in PC basis

aveX0 = mean( X0, 1 );   % aveX0 (1 x Nt) time course of average
aveX0dn = mean( X0dn, 1 );   % dnaveX0 (1 x Nt) time course of denoised average

pave_on_v = nan(1,Nv);
pavedn_on_v = nan(1,Nv);

for iv=1:Nv
    pave_on_v(iv) = Y0(iv,:) * aveX0'; 
    pavedn_on_v(iv) = Y0(iv,:) * aveX0dn'; 
end

figure;
hold on;
bar( abs( pave_on_v ) );              % plot variance of principal components
bar( abs( pavedn_on_v ), 'r' );
hold off;
xlabel( 'PCs');
ylabel( 'projection' );

title( 'timecourses: average projected on PCs' );

figure(6);
figure(5);
figure(4);
figure(3);
figure(2);
figure(1);
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

%% load clubes

function cube = LoadCube( fn, fs, rs )
% fs  = { 'trk4'; 'trk2'; 'rest' };
% rs  = {'1'; '2' );

  load( fn );

% MOT2Run1_cube          RestingStateRun1_cube  subject                
% MOT4Run1_cube          RestingStateRun2_cube  
% MOT4Run2_cube          label_AAL   

cs = 'cube = ';

if     strcmp( fs, 'trk4')
    cs = [cs 'MOT4'];
elseif strcmp( fs, 'trk2')
    cs = [cs 'MOT2'];
elseif strcmp( fs, 'rest')
    cs = [cs 'RestingState'];
else
    'oops 1'
    return;
end

if     strcmp( rs, '1')
    cs = [cs 'Run1_cube;'];
elseif strcmp( rs, '2')
    cs = [cs 'Run2_cube;'];
else
    'oops 2'
    return;
end

eval( cs );


return;
    