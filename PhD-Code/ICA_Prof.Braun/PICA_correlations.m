function PICA_correlations

clear all;
close all;
fs = 14;

load('Low-High-Subject-1-FC-ICA-High-Rest-Run-1.mat');

who 

%% evaluate high attention run

NoA = 212;
NpA = 218;

size( high_run_1.rho )
NoA + NpA

OOA = high_run_1.rho(1:NoA,1:NoA);

PPA = high_run_1.rho(NoA+1:NoA+NpA,NoA+1:NoA+NpA);

OPA = high_run_1.rho(1:NoA, NoA+1:NoA+NpA); 

%% cumulative sum of absolute correlations

NA = min( [NoA NpA] );  

SOOA = zeros(1,NA);
SPPA = zeros(1,NA);
SOPA = zeros(1,NA);

for nsc = 2:NA    % number of components to be summed
    
    for ic = 2 : nsc   % loop over components in i-direction
        for jc = 1 : ic-1   % loop over componnets in j-direction
            
            SOOA(nsc) = SOOA(nsc) + abs( OOA(ic,jc) );
            SPPA(nsc) = SPPA(nsc) + abs( PPA(ic,jc) );
            SOPA(nsc) = SOPA(nsc) + abs( OPA(ic,jc) );
            
        end
    end
    
    SOOA(nsc) = 2 * SOOA(nsc) / (nsc*nsc);   % normalization
    SPPA(nsc) = 2 * SPPA(nsc) / (nsc*nsc);   % normalization
    SOPA(nsc) = 2 * SOPA(nsc) / (nsc*nsc);   % normalization
end
    
figure;
subplot(1,2,1);
hold on;
plot(1:NA, SOOA, 'k', 'LineWidth', 2 );
plot(1:NA, SPPA, 'b', 'LineWidth', 2 );
plot(1:NA, SOPA, 'r', 'LineWidth', 2 );
hold off;
h = legend('O-O', 'P-P', 'O-P');
set(h,'FontSize', fs, 'Location', 'NorthEast' );
axis 'square';
axis([0 200 0 0.5]);
title('high attention', 'FontSize', fs);
xlabel( 'component no', 'FontSize', fs);
ylabel( 'mean abs corr', 'FontSize', fs );


%% evaluate resting run

NoR = 231;
NpR = 250;

size( rest_run_1.rho )
NoR + NpR

OOR = rest_run_1.rho(1:NoR,1:NoR);

PPR = rest_run_1.rho(NoR+1:NoR+NpR,NoR+1:NoR+NpR);

OPR = rest_run_1.rho(1:NoR, NoR+1:NoR+NpR); 

%% cumulative sum of absolute correlations

NR = min( [NoR NpR] );  

SOOR = zeros(1,NR);
SPPR = zeros(1,NR);
SOPR = zeros(1,NR);

for nsc = 2:NR    % number of components to be summed
    
    for ic = 2 : nsc   % loop over components in i-direction
        for jc = 1 : ic-1   % loop over componnets in j-direction
            
            SOOR(nsc) = SOOR(nsc) + abs( OOR(ic,jc) );
            SPPR(nsc) = SPPR(nsc) + abs( PPR(ic,jc) );
            SOPR(nsc) = SOPR(nsc) + abs( OPR(ic,jc) );
            
        end
    end
    
    SOOR(nsc) = 2 * SOOR(nsc) / (nsc*nsc);   % normalization
    SPPR(nsc) = 2 * SPPR(nsc) / (nsc*nsc);   % normalization
    SOPR(nsc) = 2 * SOPR(nsc) / (nsc*nsc);   % normalization
end
    

subplot(1,2,2);
hold on;
plot(1:NR, SOOR, 'k', 'LineWidth', 2 );
plot(1:NR, SPPR, 'b', 'LineWidth', 2 );
plot(1:NR, SOPR, 'r', 'LineWidth', 2 );
hold off;
h = legend('O-O', 'P-P', 'O-P');
set(h,'FontSize', fs, 'Location', 'NorthEast' );
axis 'square';
axis([0 200 0 0.5]);
title('resting state', 'FontSize', fs);
xlabel( 'component no', 'FontSize', fs);
ylabel( 'mean abs corr', 'FontSize', fs );

suptitle('Inter- vs Intra-lobe correlations');

print 'inter_vs_intra' -depsc2;

%%

return;