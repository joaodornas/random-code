function PICA_correlations

clear all;
close all;
fs = 14;

Ncomponent = 20;

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

ShowFirstAttComponents( OOA, PPA, OPA, Ncomponent );

%% evaluate resting run

NoR = 231;
NpR = 250;

size( rest_run_1.rho )
NoR + NpR

OOR = rest_run_1.rho(1:NoR,1:NoR);

PPR = rest_run_1.rho(NoR+1:NoR+NpR,NoR+1:NoR+NpR);

OPR = rest_run_1.rho(1:NoR, NoR+1:NoR+NpR); 

ShowFirstRestComponents( OOR, PPR, OPR, Ncomponent );

%% cumulative mean of absolute correlations

NA = min( [NoA NpA] );  

SOOA = zeros(1,NA);
SPPA = zeros(1,NA);
SOPA = zeros(1,NA);

for nsc = 2:NA    % number of components to be summed
    
    soo = abs( OOA(1:nsc,1:nsc) );
    spp = abs( PPA(1:nsc,1:nsc) );
    sop = abs( OPA(1:nsc,1:nsc) );
    
    
    SOOA(nsc) = ( sum(soo(:)) - sum(diag(soo)) ) / (nsc*nsc-nsc);
    SPPA(nsc) = ( sum(spp(:)) - sum(diag(spp)) ) / (nsc*nsc-nsc);
    SOPA(nsc) = ( sum(sop(:))                  ) / (nsc*nsc);
end
    
figure(3);
subplot(1,2,1);
hold on;
plot(1:NA, SOOA, 'k', 'LineWidth', 2 );
plot(1:NA, SPPA, 'b', 'LineWidth', 2 );
plot(1:NA, SOPA, 'r', 'LineWidth', 2 );
hold off;
h = legend('O-O', 'P-P', 'O-P');
set(h,'FontSize', fs, 'Location', 'NorthEast' );
hold on;
plot(1:NA, SOOA, 'ko', 'LineWidth', 2 );
plot(1:NA, SPPA, 'bo', 'LineWidth', 2 );
plot(1:NA, SOPA, 'ro', 'LineWidth', 2 );
hold off;
axis 'square';
axis([0 Ncomponent 0 1]);
title('high attention', 'FontSize', fs);
xlabel( 'component no', 'FontSize', fs);
ylabel( 'mean abs corr', 'FontSize', fs );




%% NON-CUMULATIVE sum of absolute correlations

NR = min( [NoR NpR] );  

SOOR = zeros(1,NR);
SPPR = zeros(1,NR);
SOPR = zeros(1,NR);

for nsc = 2:NR    % number of components to be summed
    
    soo = abs( OOR(1:nsc,1:nsc) );
    spp = abs( PPR(1:nsc,1:nsc) );
    sop = abs( OPR(1:nsc,1:nsc) );
    
    
    SOOR(nsc) = ( sum(soo(:)) - sum(diag(soo)) ) / (nsc*nsc-nsc);
    SPPR(nsc) = ( sum(spp(:)) - sum(diag(spp)) ) / (nsc*nsc-nsc);
    SOPR(nsc) = ( sum(sop(:))                  ) / (nsc*nsc);
end
    
subplot(1,2,2);
hold on;
plot(1:NR, SOOR, 'k', 'LineWidth', 2 );
plot(1:NR, SPPR, 'b', 'LineWidth', 2 );
plot(1:NR, SOPR, 'r', 'LineWidth', 2 );
hold off;
h = legend('O-O', 'P-P', 'O-P');
set(h,'FontSize', fs, 'Location', 'NorthEast' );
hold on;
plot(1:NR, SOOR, 'ko', 'LineWidth', 2 );
plot(1:NR, SPPR, 'bo', 'LineWidth', 2 );
plot(1:NR, SOPR, 'ro', 'LineWidth', 2 );
hold off;
axis 'square';
axis([0 Ncomponent 0 1]);
title('resting state', 'FontSize', fs);
xlabel( 'component no', 'FontSize', fs);
ylabel( 'mean abs corr', 'FontSize', fs );

suptitle('Inter- vs Intra-lobe correlations');

print 'inter_vs_intra_2' -depsc2;




return;

%%
function ShowFirstAttComponents(  OOA, PPA, OPA, Ncomponent )

figure(1);
subplot(1,3,1);
hold on;
caxis([-0.5 0.5]);
pcolor( OOA(1:Ncomponent,1:Ncomponent) );
hold off;
axis 'equal'; 
axis 'off';
title('OOA');

subplot(1,3,2);
hold on;
caxis([-0.5 0.5]);
pcolor( PPA(1:Ncomponent,1:Ncomponent) );
hold off;
axis 'equal'; 
axis 'off';
title('PPA');


subplot(1,3,3);
hold on;
caxis([-0.5 0.5]);
pcolor( OPA(1:Ncomponent,1:Ncomponent) );
hold off;
axis 'square'; 
axis 'off';
title('OPA');

return;

%%
function ShowFirstRestComponents(  OOR, PPR, OPR, Ncomponent )

figure(2);
subplot(1,3,1);
hold on;
caxis([-0.5 0.5]);
pcolor( OOR(1:Ncomponent,1:Ncomponent) );
hold off;
axis 'equal'; 
axis 'off';
title('OOR');

subplot(1,3,2);
hold on;
caxis([-0.5 0.5]);
pcolor( PPR(1:Ncomponent,1:Ncomponent) );
hold off;
axis 'equal'; 
axis 'off';
title('PPR');


subplot(1,3,3);
hold on;
caxis([-0.5 0.5]);
pcolor( OPR(1:Ncomponent,1:Ncomponent) );
hold off;
axis 'square'; 
axis 'off';
title('OPR');

return;

