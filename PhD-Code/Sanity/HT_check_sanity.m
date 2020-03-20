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

% difference values

dCi = [diff(Ci) 0];


dDi = [diff(Di) 0];

%% scatter plots of raw signal (D vs dC)

figure;

hold on;

plot( Di, dCi, 'b.' );

hold off;

axis 'square';

xlabel('D', 'FontSize', fs );
ylabel('\Delta C', 'FontSize', fs );



%% Hilbert transform

ASCi = hilbert( Ci );  % analytical signal

ASDi = hilbert( Di );

RCi = abs( ASCi );

RDi = abs( ASDi );

PCi = atan2( imag( ASCi ), real( ASCi ) );

PDi = atan2( imag( ASDi ), real( ASDi ) );




%% histogram of phase differences

dPi = PCi - PDi;

% keep to interval -pi to pi
kk = find( dPi > pi );
dPi(kk) = dPi(kk) - 2*pi;

kk = find( dPi < -pi );
dPi(kk) = dPi(kk) + 2*pi;


X = linspace(-5*pi/2,5*pi/2,81);

dX = X(2) - X(1);

midx = find( abs(X)== min(abs(X)) );

[NdP,X] = hist( dPi, X );

% convert to density

dP = NdP / ( length(dPi) * dX );

figure;

hold on;

h = bar( X, dP );

set(h, 'FaceColor', [1 0 0]);


hold off;

axis([X(1) X(end) 0 0.1]);

xlabel('\Phi_C - \Phi_D', 'FontSize', fs );
ylabel('density', 'FontSize', fs );



%% remove steps

dPCi = diff( PCi );

dPDi = diff( PDi );

% histogram of step sizes BEFORE removal

ShowStepSizes( dPCi, dPDi );
title( '\Delta\Phi BEFORE removing jumps', 'FontSize', fs );

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
title( '\Delta\Phi AFTER removing jumps', 'FontSize', fs );

% compare number of removed jumps with number of events

[Nevent NCstep NDstep]


%% compare original signal, Hilbert phases, and Hilbert frequencies (without jumps)

figure;
subplot(3,1,1);

ii = 50000:60000;
hold on;
plot( Ti(ii), Ci(ii), 'r', 'LineWidth', 2 );

plot( Ti(ii), Di(ii), 'b', 'LineWidth', 2 );
hold off;
h = legend( 'C', 'D');
set( h, 'FontSize', fs );

xlabel('time steps', 'FontSize', fs );
ylabel('signal', 'FontSize', fs );



subplot(3,1,2);

hold on;
plot( Ti(ii), PCi(ii), 'r', 'LineWidth', 2 );

plot( Ti(ii), PDi(ii), 'b', 'LineWidth', 2 );
hold off;
h = legend( 'C', 'D');
set( h, 'FontSize', fs );

xlabel('time steps', 'FontSize', fs );
ylabel('HT \Phi', 'FontSize', fs );


subplot(3,1,3);

hold on;
plot( Ti(ii), dPDi(ii), 'b.', 'LineWidth', 1 );

plot( Ti(ii), dPCi(ii), 'r.', 'LineWidth', 1 );
hold off;
h = legend( 'C', 'D');
set( h, 'FontSize', fs );

axis( [Ti(ii(1)) Ti(ii(end)) -0.01 0.01] );

xlabel('time steps', 'FontSize', fs );
ylabel('HT \Phi', 'FontSize', fs );


suptitle('D causes C: Signals, HT phases, HT frequencies');


print 'HTilluA' -dpdf;

%% scatter plots of unwrapped HT frequencies

figure;

hold on;

plot( dPCi, dPDi, 'b.' );  % frequencies (difference per timestep)

hold off;

axis 'square';

xlabel('\Delta\Phi_C', 'FontSize', fs );
ylabel('\Delta\Phi_D', 'FontSize', fs );

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

suptitle( 'HT \Phi without jumps, detrended and not' );


%% plot orbits around selected transition events


idxup   = find( Di(Tevent+1) - Di(Tevent-1) ==  2 );  % index to Di, Ci

idxdown = find( Di(Tevent+1) - Di(Tevent-1) == -2 );  % index to Di, Ci

Nselect = 200;   % must be smaller than Nevent

idx = randperm(length(idxup));  % index to idxup       

TselectUP   = Tevent( idxup( idx(1:Nselect) ) );   % choose Nselect UP events

TselectDOWN = Tevent( idxdown( idx(1:Nselect) ) ); % choose Nselect DOWN events

fullwindow = 51;
halfwindow = 25;

rt = -halfwindow:halfwindow;
%

figure;

subplot(1,2,1);

hold on;
for i=1:Nselect
    
    ievent = TselectUP(i);
    
    plot3( dPCi(ievent+rt), dPDi(ievent+rt), rt, 'r', 'LineWidth', 1 );
    
    plot3( dPCi(ievent), dPDi(ievent), 0, 'k.');
end
hold off;

axis( [-pi/8 pi/8 -pi/8 pi/8 -halfwindow halfwindow] );

view( 50, -15 );

xlabel('\Delta \Phi_C', 'FontSize', fs );
ylabel('\Delta \Phi_D', 'FontSize', fs );
zlabel('t', 'FontSize', fs );

title(' UP', 'FontSize', fs );

subplot(1,2,2);

hold on;
for i=1:Nselect
    
    ievent = TselectDOWN(i);
    
    plot3( dPCi(ievent+rt), dPDi(ievent+rt), rt, 'g', 'LineWidth', 1 );
    
    plot3( dPCi(ievent), dPDi(ievent), 0, 'k.');
end

hold off;

axis( [-pi/8 pi/8 -pi/8 pi/8 -halfwindow halfwindow] );

view( 50, -15 );

xlabel('\Delta \Phi_C', 'FontSize', fs );
ylabel('\Delta \Phi_D', 'FontSize', fs );
zlabel('t', 'FontSize', fs );

title('DOWN', 'FontSize', fs );






suptitle( 'Trajectories around reversal (black dot)' );
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

axis([X(1) X(end) 0 0.01]);

xlabel('\Delta \Phi', 'FontSize', fs );
ylabel('density', 'FontSize', fs );

return;