function goforit

load( 'AAL_correlation_matrix.mat' );

path( path, 'export_fig' );
                                           

cv_z = ones( size( mu_z ) );
kk = find( abs( mu_z ) > 0 );
cv_z(kk) = sigma_z(kk) ./ abs( mu_z(kk) );

mu_z_clean = abs( mu_z );
mu_z_clean( abs( mu_z ) <= z_threshold ) = z_threshold;
mu_z_clean( cv_z        >= cv_threshold ) = z_threshold;

cmin = z_threshold;
cmax = 1;

cmap = MyColormap;

hf = figure('color','w');
hold on;
hm = pcolor( mu_z_clean );
colormap( cmap );
set(hm, 'EdgeColor',[1 1 1], 'LineWidth', 0.1);
caxis([cmin cmax]);
get(hm)
hold off;
axis 'equal';
axis 'off';

% set( hf, 'Units', 'Inches' );
% pos = get( hf, 'Position' );
% set( hf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)] );

export_fig('AAL_correlation_matrix','-pdf');
%print( hf, 'AAL_correlation_matrix', '-dpdf', '-loose', '-r0' );

%print 'AAL_correlation_matrix' -dpdf;


return;

function cmap = MyColormap

cmap = zeros(64,3);

cmap(1,1) = 1; cmap(1,2) = 1; cmap(1,3) = 1;

cmap(2:8,1)   = 0;
cmap(2:8,2)   = 0;
cmap(2:8,3)   = (8 + (2:8)) / 16.0;    % dark blue to blue

cmap(9:24,1)  = 0;
cmap(9:24,2)  = (1:16) / 16.0;        %  blue to cyan
cmap(9:24,3)  = 1;

cmap(25:40,1) = (1:16) / 16;
cmap(25:40,2) = 1;                    %  cyan to yellow
cmap(25:40,3) = (15:-1:0) / 16;

cmap(41:56,1) = 1;
cmap(41:56,2) = (15:-1:0) / 16;       %  yellow to red
cmap(41:56,3) = 0;

cmap(57:64,1) = (15:-1:8) / 16;
cmap(57:64,2) = 0;                    %  red to dark red
cmap(57:64,3) = 0;

return;