function plotFCgeneral(mu_z,sigma_z,label,plotParcelsLabels_idx,ROI)

doThePlot(mu_z,sigma_z,label,plotParcelsLabels_idx,ROI);

end

function doThePlot(mu_z,sigma_z,label,plotParcelsLabels_idx,ROI)
                        
mu_z(isinf(mu_z)) = 0;

z_threshold = 0.1889;
cv_threshold = 1;

cv_z = ones( size( mu_z ) );
kk = find( abs( mu_z ) > 0 );
cv_z(kk) = sigma_z(kk) ./ abs( mu_z(kk) );

% l = size(mu_z,1);
% b = [mu_z,zeros(l,1)];
% mu_z = [b;zeros(1,l+1)];

mu_z_clean = abs( mu_z );
mu_z_clean( abs( mu_z ) <= z_threshold ) = z_threshold;
mu_z_clean( cv_z        >= cv_threshold ) = z_threshold;

cmin = z_threshold;
cmax = 1;

% cmap = MyColormap;
cmap = MyColormapNaN;

hf = figure('color','w');
hold on;
% hm = pcolor( mu_z_clean );
hm = imagesc( mu_z_clean );
colormap( cmap );
% set(hm, 'EdgeColor',[1 1 1], 'LineWidth', 0.00000000001);
caxis([cmin cmax]);

% set( hf, 'Units', 'Inches' );
% pos = get( hf, 'Position' );
% set( hf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)] );
    
switch plotParcelsLabels_idx

    case 1 %%% plot FC MD758 Clusters labels script

        run_FC_MD758_Clusters_labels_script;

    case 2 %%% plot FC JB Clusters labels script

        run_FC_JB_Clusters_labels_script;

end

get(hm);
hold off;
axis 'equal';
axis 'off';

set(gcf, 'color', 'none');
set(gca, 'color', 'none');

set(hm,'AlphaData',mu_z_clean > cmin); 
 
export_fig(label,'-transparent');

% print( hf, label, '-dpdf');
% print( hf, label, '-depsc2');
% save2pdf(label,hf,600);

%print 'AAL_correlation_matrix' -dpdf;

close all

end

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

end

function cmap = MyColormapNaN

cmap = zeros(64,3);

cmap(1,1) = NaN; cmap(1,2) = NaN; cmap(1,3) = NaN;

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

end