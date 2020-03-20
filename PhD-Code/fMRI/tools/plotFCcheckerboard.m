function plotFCcheckerboard( log_mu, log_sigma, cmin, cmax, log_mu_threshold, applyThreshold, applyColorBar, label )

[Clusters_per_BIGROI, BIGROI_labels, BIGROI_of_Cluster] = GetInfoBIGROI();

if applyThreshold == 2
    
    log_mu(isinf(log_mu)) = 0;

    z_threshold = 0.1889;   
    cv_threshold = 1;

    cv_z = ones( size( log_mu ) );
    kk = find( abs( log_mu ) > 0 );
    cv_z(kk) = log_sigma(kk) ./ abs( log_mu(kk) );

    log_mu_clean = abs( log_mu );
    log_mu_clean( abs( log_mu ) <= z_threshold ) = nan;
    log_mu_clean( cv_z        >= cv_threshold ) = nan;

    cmin = z_threshold;
    cmax = 1;
    
    csub0 = cmin - (2/62) * ( cmax - cmin );
    csub1 = cmin - (0.5/62) * ( cmax - cmin );

elseif applyThreshold == 1
    
    log_mu(isinf(log_mu)) = nan;

    log_mu_clean = log_mu;
    log_mu_clean( log_mu < log_mu_threshold ) = nan;

%     cmin = log_mu_threshold;
%     cmax = max( log_mu_clean(:) );
    
    csub0 = cmin - (2/62) * ( cmax - cmin );
    csub1 = cmin - (0.5/62) * ( cmax - cmin );

elseif applyThreshold == 0
    
    log_mu(isinf(log_mu)) = nan;
    
    log_mu_clean = log_mu;
    
    csub0 = cmax - (32/64) * ( cmax - cmin );
    csub1 = cmax - (33/64) * ( cmax - cmin );
    
end

                                     % fill in checkerboard
[ifcl, jfcl] = meshgrid(BIGROI_of_Cluster);
log_mu_clean( isnan(log_mu_clean) & mod(ifcl,2)==mod(jfcl,2) ) = csub0;
log_mu_clean( isnan(log_mu_clean) & mod(ifcl,2)~=mod(jfcl,2) ) = csub1;

if applyThreshold == 1 || applyThreshold == 2
    
    [n,m] = size( log_mu_clean );        % extend array size by one ...
    % log_mu_extended = log_mu_threshold * ones( n+1, m+1 );
    log_mu_extended = zeros( n+1, m+1 );
    log_mu_extended(1:n, 1:m) = log_mu_clean;
    
else
   
    [n,m] = size( log_mu_clean );        % extend array size by one ...
    log_mu_extended = zeros( n+1, m+1 );
    log_mu_extended(1:n, 1:m) = log_mu_clean;
    
end

cmap = MyJetColormapForCheckerboard(applyThreshold);

hf = figure('color','w');
hold on;
hm = imagesc( log_mu_extended );
colormap( cmap );
%set(hm, 'EdgeColor',[1 1 1], 'LineWidth', 0.1);

if applyThreshold
    caxis([csub0 cmax]);
else
    caxis([cmin cmax]);
end

hold off;
axis 'equal';
set( gca, 'XTick', [1 m], 'YTick', [1 n] );

switch applyColorBar
    
    case 1
        
        hb = colorbar;
        % set(hb, 'YTick', log([1 10 100 1000]), 'YTickLabel', {'0'; '1'; '2'; '3'} );

end

get(hm);
hold off;
axis 'equal';
axis 'off';
 
export_fig(label,'-pdf');

close all

switch applyColorBar
    
    case 2
        
        hf = figure('color','w');
        cmap = MyJetColormapForCheckerboard(applyThreshold);
        hold on;
        colormap( cmap );
        
        if applyThreshold
            caxis([csub0 cmax]);
        else
            caxis([cmin cmax]);
        end

        hb = colorbar;
        % set(hb, 'YTick', log([1 10 100 1000]), 'YTickLabel', {'0'; '1'; '2'; '3'} );
        
        export_fig(strcat(label,'-colorbar'),'-pdf');

end

return;

