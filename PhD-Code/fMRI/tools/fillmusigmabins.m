function [mu_med, sigma_med, D_mu_sigma] = FillMuSigmaBins( mu_z, sigma_z, nbin, bintype )


if     strcmp( bintype, 'Quantile' )
    
    [mu_bin, mu_med, mu_del]          = MakeQuantileBins( mu_z, nbin );   % x_bin[ 1 x (nbin+1) ], x_med[ 1 x nbin ]
    
    [sigma_bin, sigma_med, sigma_del] = MakeQuantileBins( sigma_z, nbin );
    
elseif strcmp( bintype, 'Linear' )
    
    [mu_bin, mu_med, mu_del]          = MakeLinearBins( mu_z, nbin );   % x_bin[ 1 x (nbin+1) ], x_med[ 1 x nbin ]
    
    [sigma_bin, sigma_med, sigma_del] = MakeLinearBins( sigma_z, nbin );
    
else
    'bintype not recognized'
    return;
end

D_mu_sigma = zeros( nbin, nbin );

for i = 1 : nbin
    
    
    ii = find( mu_bin(i) < mu_z & mu_z <= mu_bin(i+1) );
    
    
    if ~isempty( ii )
        
        for j = 1 : nbin
            
            kk = find( sigma_bin(j) < sigma_z(ii) & sigma_z(ii) <= sigma_bin(j+1) );
            
            if ~isempty(kk)
                
                D_mu_sigma( i, j ) = numel(kk);
                
            end
        end
        
    else
        ['empty bin ' num2str(i)]
    end
end


for i = 1 : nbin                                  % convert to densities
    
    D_mu_sigma( i, : ) = D_mu_sigma( i, : ) / mu_del(i);
    D_mu_sigma( :, i ) = D_mu_sigma( :, i ) / sigma_del(i);  
     
end
    
return;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_bin, x_med, x_del] = MakeLinearBins( x, nbin ) 

    x_med = zeros(1,nbin);

    x_bin = linspace(min(x), max(x), nbin+1 );               % linear bins
    
    for i=1 : nbin                                  
        x_med(i) = median( x( x_bin(i) <= x & x <= x_bin(i+1) ) );                          % bin medians
    end
    
    x_del = x_bin(2:end) - x_bin(1:end-1);
return;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_bin, x_med, x_del] = MakeQuantileBins( x, nbin ) 

    x_med = zeros(1,nbin);

    x_bin = zeros( 1, nbin+1 );
    x_bin_inner     = quantile( x, nbin-3 );           % quantile bins
    x_bin(3:nbin-1) = x_bin_inner;
    
    
    kk = find( x < x_bin(3) );                         % first two bins    
    x_bin(2)        = median( x(kk) );   
    x_bin(1)        = min(x(kk));
    
    kk = find( x_bin(nbin-1) < x );                    % last two bins
    x_bin(nbin)     = median( x(kk) );  
    x_bin(nbin+1)   = max(x(kk));  
    
    for i=1 : nbin                                  
        x_med(i)    = median( x( x_bin(i) <= x & x <= x_bin(i+1) ) );                          % bin medians
    end

    x_del = x_bin(2:end) - x_bin(1:end-1);
    
return;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
