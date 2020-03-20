function ShowConnectivityMatrixWithCheckerboard_example( log_mu, log_mu_threshold )

[~, ~, BIGROI_of_Cluster] = GetInfoBIGROI();


log_mu_clean = log_mu;
log_mu_clean( log_mu < log_mu_threshold ) = nan;

cmin = log_mu_threshold;
cmax = max( log_mu_clean(:) );
csub0 = cmin - (2/62) * ( cmax - cmin );
csub1 = cmin - (0.5/62) * ( cmax - cmin );

                                     % fill in checkerboard
[ifcl, jfcl] = meshgrid(BIGROI_of_Cluster);
log_mu_clean( isnan(log_mu_clean) & mod(ifcl,2)==mod(jfcl,2) ) = csub0;
log_mu_clean( isnan(log_mu_clean) & mod(ifcl,2)~=mod(jfcl,2) ) = csub1;


[n,m] = size( log_mu_clean );        % extend array size by one ...
log_mu_extended = log_mu_threshold * ones( n+1, m+1 );
log_mu_extended(1:n, 1:m) = log_mu_clean;

cmap = MyJetColormapForCheckerboard(1);

figure;
hold on;
h = imagesc( log_mu_extended );
colormap( cmap );
%set(h, 'EdgeColor',[1 1 1], 'LineWidth', 0.1);
caxis([csub0 cmax]);

hold off;
axis 'equal';
set( gca, 'XTick', [1 m], 'YTick', [1 n] );

hb = colorbar;
set(hb, 'YTick', log([1 10 100 1000]), 'YTickLabel', {'0'; '1'; '2'; '3'} );

%log_mu_clean( log_mu_clean == log_mu_threshold ) = 0;

return;

