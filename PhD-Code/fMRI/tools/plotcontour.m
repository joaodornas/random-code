function plotcontour(label, bintype, pairtype, x_bin, y_bin, D_xy)

f = figure;
r_threshold = 0.186728;  % corresponds to p-value = 0.01 for N=150
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );

hold on;
h = fill(z_threshold*[-1 1 1 -1], [0 0 0.55 0.55], [0.9 0.9 0.9] );
% h = fill(z_threshold*[-1 1 1 -1], [0 0 0.35 0.35], [0.9 0.9 0.9] ); %%
% STANDARD FILL
set(h,'EdgeColor', [1 1 1]);
contour( x_bin, y_bin, D_xy', 'LineWidth', 1 );
colormap('jet');
%colorbar;
plot([0 1],[0 1.0],'k--');
% plot([0 1],[0 0.1],'k--');

axis 'square';
xlabel('\mu_z', 'FontSize', 20);
ylabel('\sigma_z', 'FontSize', 20);
xlim([-0.2 1.0]); %% STANDARD INTERVAL
% xlim([-0.05 0.4]);
% ylim([0 0.35]);
ylim([0 0.55]); %% STANDARD INTERVAL
% set( gca, 'XTick', [0 0.4 0.8], 'YTick', [0 0.1 0.2 0.3 0.4 0.5 0.6], 'FontSize', 20 );
% ylim([0 20000]);
set( gca, 'XTick', [0 0.4 0.8], 'YTick', [0 0.1 0.2 0.3 0.4 0.5], 'FontSize', 20 );
% set( gca, 'XTick', [0 0.4 0.8], 'YTick', [0 0.1 0.2 0.3], 'FontSize', 20 ); %% STANDARD TICK
title([label '-' pairtype]);

print(f,'-depsc',strcat('Contour-',label,'-',bintype, '-', pairtype, '.eps'));

return;