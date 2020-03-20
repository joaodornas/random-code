function fraction = FractionSignificantConsistent( mu_z, sigma_z, z_threshold, cv_threshold )

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


figure;
subplot(1,3,1);
hist( mu_z_sorted(kk) );
xlabel('\mu_z');

subplot(1,3,2);
hist( sigma_z_sorted(kk) );
xlabel('\sigma_z');

subplot(1,3,3);
hist( cv_z_sorted(kk) );
xlabel('c_v');

return;