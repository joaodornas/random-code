function HT_check_sanity

clear all;
close all;

fs = 14;


%% run Hilbert transforms on sanity-check data


% retrieve values Ci and Di to file 

load( 'sanity.mat', 'Ci', 'Di', 'Tevent' );

% time axis

Tstart = 1;
Tend   = length(Ci);
Ti     = 1:Tend;

Tend

% number of transition events

Nevent = length(Tevent);

Nevent

%% Hilbert transform

ASCi = hilbert( Ci );  % analytical signal

ASDi = hilbert( Di );

RCi = abs( ASCi );

RDi = abs( ASDi );

PCi = atan2( imag( ASCi ), real( ASCi ) );

PDi = atan2( imag( ASDi ), real( ASDi ) );


figure;
subplot(2,1,1);

hold on;
plot( Ti(1:10000), Ci(1:10000), 'r', 'LineWidth', 2 );

plot( Ti(1:10000), Di(1:10000), 'b', 'LineWidth', 2 );
hold off;
h = legend( 'C', 'D');
set( h, 'FontSize', fs );

xlabel('time steps', 'FontSize', fs );
ylabel('signal', 'FontSize', fs );

subplot(2,1,2);

hold on;
plot( Ti(1:10000), PCi(1:10000), 'r', 'LineWidth', 2 );

plot( Ti(1:10000), PDi(1:10000), 'b', 'LineWidth', 2 );
hold off;
h = legend( 'C', 'D');
set( h, 'FontSize', fs );

xlabel('time steps', 'FontSize', fs );
ylabel('HT \Phi', 'FontSize', fs );


%% remove steps

dPCi = diff( PCi );

dPDi = diff( PDi );

% histogram of step sizes BEFORE removal

ShowStepSizes( dPCi, dPDi );
title( 'before', 'FontSize', fs );

% remove steps larger pi and smaller -pi

kk = find( dPCi < -pi );
dPCi(kk) = dPCi(kk) + 2*pi;

jj = find( dPCi >  pi );
dPCi(jj) = dPCi(jj) - 2*pi;

NCstep = length(kk) + length(jj);

kk = find( dPDi < -pi );
dPDi(kk) = dPDi(kk) + 2*pi;

jj = find( dPDi >  pi );
dPDi(jj) = dPDi(jj) - 2*pi;

NDstep = length(kk) + length(jj);

% histogram of step sizes AFTER removal

ShowStepSizes( dPCi, dPDi );
title( 'after', 'FontSize', fs );

% compare number of removed jumps with number of events

[Nevent NCstep NDstep]



%% subtract average rate of phase increase (detrending)

mdPCi = mean( dPCi );
mdPDi = mean( dPDi );


% reconstitute detrended signals

rdtPCi = cumsum( [0 dPCi - mdPCi] );

rdtPDi = cumsum( [0 dPDi - mdPDi] );


% note that mean is not zero!  only final value is guaranteed to be zero
% ...


% reconstitute full signals

rPCi = rdtPCi + Ti * mdPCi;

rPDi = rdtPDi + Ti * mdPDi;


figure;

subplot(2,1,1);

hold on;
plot( Ti, rdtPCi, 'r', 'LineWidth', 1 );

plot( Ti, rdtPDi, 'b', 'LineWidth', 1 );

h = legend( 'C', 'D');
set( h, 'FontSize', fs );

plot( Ti(Tevent), rdtPCi(Tevent), 'r.' );
plot( Ti(Tevent), rdtPDi(Tevent), 'b.' );

hold off;

xlabel('time steps', 'FontSize', fs );
ylabel('detrended HT \Phi', 'FontSize', fs );

subplot(2,1,2);

hold on;
plot( Ti, rPCi, 'r', 'LineWidth', 1 );

plot( Ti, rPDi, 'b', 'LineWidth', 1 );

h = legend( 'C', 'D');
set( h, 'FontSize', fs );

plot( Ti(Tevent), rPCi(Tevent), 'r.' );
plot( Ti(Tevent), rPDi(Tevent), 'b.' );

hold off;

xlabel('time steps', 'FontSize', fs );
ylabel('HT \Phi', 'FontSize', fs );


%% plot orbits around selected transition events

Nselect = 20;   % must be smaller than Nevent

idx = randperm(Nevent);          

Tselect = Tevent( idx(1:Nselect) );   % choose Nselect of Nevent events

fullwindow = 51;
halfwindow = 25;


figure;

subplot(1,2,1);

hold on;
for i=1:Nselect
    
    ievent = Tselect(i);
    
    iwindow = ievent+[-halfwindow:halfwindow];
    
    rdtPC = rdtPCi(iwindow) - rdtPCi(ievent);
    rdtPD = rdtPDi(iwindow) - rdtPDi(ievent);
    
    plot( rdtPC, rdtPD, 'r', 'LineWidth', 1 );
end

plot( 0, 0, 'k.' );

hold off;

axis 'equal';
axis 'square';

xlabel('\Delta \Phi_C', 'FontSize', fs );
ylabel('\Delta \Phi_D', 'FontSize', fs );

subplot(1,2,1);

hold on;
for i=1:Nselect
    
    ievent = Tselect(i);
    
    iwindow = ievent+[-halfwindow:halfwindow];
    
    rdtPC = rdtPCi(iwindow) - rdtPCi(ievent);
    rdtPD = rdtPDi(iwindow) - rdtPDi(ievent);
    
    plot( rdtPC, rdtPD, 'r', 'LineWidth', 1 );
end

plot( 0, 0, 'k.' );

hold off;

axis 'square';

xlabel('\Delta \Phi_C', 'FontSize', fs );
ylabel('\Delta \Phi_D', 'FontSize', fs );

title('detrended', 'FontSize', fs );

subplot(1,2,2);

hold on;
for i=1:Nselect
    
    ievent = Tselect(i);
    
    iwindow = ievent+[-halfwindow:halfwindow];
    
    rPC = rPCi(iwindow) - rPCi(ievent);
    rPD = rPDi(iwindow) - rPDi(ievent);
    
    plot( rPC, rPD, 'r', 'LineWidth', 1 );
end

plot( 0, 0, 'k.' );

hold off;

axis 'equal';
axis 'square';

xlabel('\Delta \Phi_C', 'FontSize', fs );
ylabel('\Delta \Phi_D', 'FontSize', fs );

title('reconstructed', 'FontSize', fs );


return;



%% function defs

function ShowStepSizes( dPCi, dPDi )

fs = 14;

figure;

X = linspace(-5*pi/2,5*pi/2,81);

dX = X(2) - X(1);

midx = find( abs(X)== min(abs(X)) );

[NC,X] = hist( dPCi, X );

[ND,X] = hist( dPDi, X );

% convert to density

PC = NC / ( length(dPCi) * dX );
PD = ND / ( length(dPDi) * dX );

hold on;

h = bar( X, PC );

set(h, 'FaceColor', [0 0 1]);

h = bar( X, PD );

set(h, 'FaceColor', [1 0 0]);


hold off;

h = legend( 'C', 'D');
set( h, 'FontSize', fs );

axis([X(1) X(end) 0 0.1]);

xlabel('\Delta \Phi', 'FontSize', fs );
ylabel('density', 'FontSize', fs );

return;