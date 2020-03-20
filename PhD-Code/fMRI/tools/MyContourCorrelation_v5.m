function MyContourCorrelation_v5(X,Y,label,compute_bins)

X(isinf(X)) = NaN;
Y(isinf(Y)) = NaN;

X(find(X==0)) = NaN;
Y(find(Y==0)) = NaN;

mu_z = X;
sigma_z = Y;
nbin = 5;
bintype = 'Quantile';
% bintype = 'Linear';

[mu_med, sigma_med, D_mu_sigma] = fillmusigmabins( mu_z, sigma_z, nbin, bintype );

x_mid = mu_med; y_mid = sigma_med; D_xy = D_mu_sigma;

save(strcat('Contour-',label,'-',bintype,'-bins.mat'),'x_mid','y_mid','D_xy','label');

plotcontour(label, bintype, 'Mean-STD', x_mid, y_mid, D_xy);

close all


end
