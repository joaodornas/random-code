function IM_total = TotalInformationContent( mu_z, sigma_z, z_threshold, cv_threshold )

mu_z(isinf(mu_z)) = NaN;
sigma_z(isinf(sigma_z)) = NaN;

cv_z = zeros( size( mu_z ) );                % coefficient of variation (lots of spurious values for small means ...

kk = find( mu_z ~= 0 );
N_all_values = numel( kk );
cv_z(kk) = sigma_z(kk) ./ abs( mu_z(kk) );   


[cv_z_sorted, idx] = sort( cv_z, 'descend' );  % sort cv's in descending order

mu_z_sorted = mu_z( idx );

sigma_z_sorted = sigma_z( idx );


kk = find( mu_z_sorted > z_threshold & cv_z_sorted < cv_threshold );  % find significant and consistent values

N_significant_and_consistent_values = numel(kk)

N_all_values

fraction = N_significant_and_consistent_values / N_all_values


%% obtain significant and consistent z-values


mu_z_contributing = mu_z_sorted(kk);

rho_contributing = tanh( mu_z_contributing);


%% Mutual information of a bivariate normal distribution with correlation coefficient rho

% IM = - 0.5 * log2( 1 - rho^2 );

IM_contribution = -0.5 * log2( 1 - rho_contributing.^2 );



%% total information content equals expected IC (per observation), times fraction of significant & consistent observations, times number of observations

IM_total = nansum( IM_contribution(:) );

['Total information content ' num2str( IM_total/8000 ) ' kBytes']

return;