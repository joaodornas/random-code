
% voxels_corr_z = (1/2)*log((1 + voxels_corr)./(1 - voxels_corr));

mu = squeeze(nanmean(voxels_corr_rest,1));
sigma = squeeze(nanstd(voxels_corr_rest));

mu_z = (1/2)*log((1 + mu)./(1 - mu));
sigma_z = (1/2)*log((1 + sigma)./(1 - sigma));

iTotal=10;

for i=1:iTotal

    % label = strcat('FC-20-Voxels-',int2str(i));
    label = 'FC-20-Voxels-Total';
    
    %plotFCgeneral(mu_z((1+(i-1)*240):(240+(i-1)*240),(1+(i-1)*240):(240+(i-1)*240)),sigma_z((1+(i-1)*240):(240+(i-1)*240),(1+(i-1)*240):(240+(i-1)*240)),label);

    z_threshold = 0.1889;
 
    cmin = 0;
    cmax = 1;

    cmap = colormap('jet');
    cmap(1,:) = [1 1 1];
    
    hf = figure('color','w');
    hold on;
    % hm = pcolor( mu_z((1+(i-1)*240):(240+(i-1)*240),(1+(i-1)*240):(240+(i-1)*240)) );
    hm = imagesc( mu_z );
    colormap( cmap );
    % set(hm, 'EdgeColor',[1 1 1], 'LineWidth', 0.00000000001);
    caxis([cmin cmax]);
    get(hm)
    hold off;
    axis 'equal';
    axis 'off';

    export_fig(label,'-pdf');

end