
clear all;
close all;

fs = 16;

path( path, './Sources.Matlab');
path( path, './Settings');
%% generate trajectories

MOT = GenerateBraunianWalkTrajectories;

xi = MOT.xi;
yi = MOT.yi;

sizeX = 400;
sizeY = 400;
repR = 50;

nFrames = size(xi,2);
sFrames = 5;   % sampling of frames
nBalls = size(xi,1);

%startStr = input('First Frame (i.e. 1): ','s');
%finishStr = input(strcat('Last Frame (max. ',int2str(nFrames),'): '),'s');

%start = str2num(startStr);
%finish = str2num(finishStr);

start = 1;
finish = 300;   % max 600 (for 10s trial duration)

%% plot trajectories

color{1} = 'r';
color{2} = 'b';
color{3} = 'g';
color{4} = 'k';
color{5} = 'y';

balls = 1:nBalls;

BallColors = MOT.BallColors;
AllRGB = MOT.Settings.MOT.Target.AllRGB;
nColors = size(AllRGB,1);

for b=1:nBalls/2

    colorBalls(b) = 1;
    
end
for b=nBalls/2+1:nBalls
    
    colorBalls(b) = 2;
end
 
figure;

hold on;
plot( [-0.5 0.5 0.5 -0.5 -0.5] * sizeX, [-0.5 -0.5 0.5 0.5 -0.5]*sizeY, 'k', 'LineWidth', 2 );
plot( [-0.5 0.5 0.5 -0.5 -0.5] * (sizeX-2*repR), [-0.5 -0.5 0.5 0.5 -0.5]*(sizeY-2*repR), 'k--', 'LineWidth', 2 );

for b=1:length(balls)
    
   plot(xi(balls(b),start:finish),yi(balls(b),start:finish),color{colorBalls(b)});
    
end
hold off;

axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

axis 'equal';
axis 'off';

print 'untitled1' -depsc2;

%% compute pairwise distances to find number of nearest neighbors with same color

Nneighbor = min( 12, nBalls );  % number of neighbors to be considered

Nsigma = 100;

Ns2 = 2*Nsigma*Nsigma;

Wsame = nan(nBalls, ceil(nFrames/sFrames) );  % sum of exponential weights of same balls
Wall  = nan(nBalls, ceil(nFrames/sFrames) );  % sum of exponential weights of all balls
Wfrac = nan(nBalls, ceil(nFrames/sFrames) );  % sum of same weights / sum of all weights
Wprob = nan(nBalls, ceil(nFrames/sFrames) );  % probability density of above fraction

for kf=1:ceil(nFrames/sFrames)
    
    jf = (kf-1)*sFrames+1;
    
    Nexp2 = nan(nBalls,nBalls);
    
    for ib1=1:nBalls-1  % compute pairwise distances
        
        x1 = xi(ib1,jf);
        y1 = yi(ib1,jf);
        
        for ib2=ib1+1:nBalls
            
            x2 = xi(ib2,jf);
            y2 = yi(ib2,jf);
            
            dsqr = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
            Nexp2(ib1, ib2) = exp( -dsqr / Ns2 );
            Nexp2(ib2, ib1) = exp( -dsqr / Ns2 );
            
        end
        
    end
    
    

    for ib1=1:nBalls;  % examine neighborhood
        
        Cball = colorBalls( ib1 );       % color of central ball
        
        [nexp2, nidx] = sort(Nexp2(ib1,:));   % sort neighbors by distance
        
        Cneigh = colorBalls( nidx );  % get colors of neighbors, same order as ndis2
        
        sidx = find( Cneigh == Cball );  % index to ndis2 of same color neighbors
        
        Wsame( ib1, kf ) = nansum( nexp2(sidx) );  % weighted sum of same
        
        Wall( ib1, kf )  = nansum( nexp2 );  % weighted sum of same
        % Nsame( ib1, kf ) = length( find( Cneigh == Cball ) );  % number of same balls
        
    end
    
end



kk = find( Wall ~= 0 );

Wfrac(kk) = Wsame(kk) ./ Wall(kk);


%% show distributions

[nsame,esame] = hist(Wsame(:),21);

fsame = nsame / sum(nsame);

dsame = fsame / (esame(2) - esame(1));

[nall,eall] = hist(Wall(:),21);

fall = nall / sum(nall);

dall = fall / (eall(2) - eall(1));

[nfrac,efrac] = hist(Wfrac(:),21);

ffrac = nfrac / sum(nfrac);

dfrac = ffrac / (eall(2) - eall(1));

figure;
subplot(3,1,1);
hold on;
bar( esame, dsame ); 
hold off;
xlabel('E same', 'FontSize', fs );
ylabel('density', 'FontSize', fs );

subplot(3,1,2);
hold on;
bar( eall, dall ); 
hold off;
xlabel('E all', 'FontSize', fs );
ylabel('density', 'FontSize', fs );

subplot(3,1,3);
hold on;
bar( efrac, dfrac ); 
hold off;
xlabel('E fraction', 'FontSize', fs );
ylabel('density', 'FontSize', fs );

%% save probabilities

Wprob = nan(size(Wfrac));

kk = find( Wfrac <= 0.5*(efrac(1)+efrac(2)) );
Wprob(kk) = efrac(1);

for ii=2:length(efrac)-1
    elo = 0.5*(efrac(ii-1)+efrac(ii));
    ehi = 0.5*(efrac(ii)+efrac(ii+1));
    kk = find( elo < Wfrac & Wfrac <= ehi );
    Wprob(kk) = efrac(ii);
end

kk = find( Wfrac > 0.5*(efrac(end-1)+efrac(end)) );
Wprob(kk) = efrac(end);


%% check the results

figure;

hold on;
for ib=1:nBalls
    
    plot( xi(ib,1), yi(ib,1), 'o', 'Color', color{colorBalls(ib)} );
    
    text( 10+xi(ib,1), yi(ib,1), [num2str(Wfrac(ib,1),'%4.2f')], 'FontSize', 12 );
    
end

hold off;

axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

axis 'equal';
axis 'off';


%% compute Entropy of distribution of counts in windows of 30 frames

%  also look for balls with complete homogeneity or heterogeneity

nWindow = 30/sFrames;

Hcolor = [];
Lfrac = [];
Hfrac = [];
count = 1;


for kf=1:ceil(nFrames/sFrames)-nWindow+1
    
    jf = (kf-1)*sFrames+1;  % index to frames
    rf = kf:kf+nWindow-1;   % index to sampled frames
    
    Rfrac = Wfrac(:,rf);   % same fraction of all balls at all sampled frames in window
    Rprob = Wprob(:,rf);
    
    % overall count probability of same neighbors
    
    Hcolor(count) = nansum( -Rprob(:) .* log( Rprob(:) ) );
    
    % homogeneity of neighbors of individual ball with most uniform
    % neighborhood
    
    arf = sort( mean( Rfrac, 2 ) );  % average weight fraction of each ball over sampled frames in window
    
    Lfrac(count) = arf(1);   % lowest fraction of all balls
    Hfrac(count) = arf(end);  % highest fraction of all balls
    
    count = count+1;
    
end

figure;
subplot(2,1,1);
hold on;
plot( 1:ceil(nFrames/sFrames)-nWindow+1, Hcolor, 'r.' );
hold off;
xlabel('frames', 'FontSize', fs );
ylabel('Entropy', 'FontSize', fs );

subplot(2,1,2);
hold on;
plot( 1:ceil(nFrames/sFrames)-nWindow+1, Lfrac, 'r.' );
plot( 1:ceil(nFrames/sFrames)-nWindow+1, Hfrac, 'b.' );
hold off;
xlabel('frames', 'FontSize', fs );
ylabel('W fraction', 'FontSize', fs );

print 'untitled2' -depsc2;

kfmin = find( Hcolor == min(Hcolor ) );   % smallest entropy over all sampled frames

kfmax = find( Hcolor == max(Hcolor ) );   % highest entropy over all sampled frames

dummy = find( Hfrac == max(Hfrac ) );     % highest fraction over all sampled frames
hnmax = dummy(1);

dummy = find( Lfrac == min(Lfrac ) );     % lowest fraction over all sampled frames
lnmin = dummy(end);

%% check the results


figure;

subplot(2,2,1);

kfrange = kfmin:kfmin+nWindow-1;  % index to sampled frames;
jfmin = (kfmin-1)*sFrames+1;
jfrange = (kfrange-1)*sFrames+1;  % index to original frames;

hold on;
for ib=1:nBalls
    
    plot( xi(ib,jfrange), yi(ib,jfrange), 'o', 'Color', color{colorBalls(ib)} );
    
    text( 10+xi(ib,jfmin), yi(ib,jfmin), [num2str(Wfrac(ib,kfmin),'%4.2f')], 'FontSize', 12 );
    
end


hold off;

axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

axis 'equal';
axis 'off';
title( 'Lowest entropy', 'FontSize', fs );

subplot(2,2,2);

kfrange = kfmax:kfmax+nWindow-1;  % index to sampled frames;
jfmax = (kfmax-1)*sFrames+1;
jfrange = (kfrange-1)*sFrames+1;  % index to original frames;

hold on;
for ib=1:nBalls
    
    plot( xi(ib,jfrange), yi(ib,jfrange), 'o', 'Color', color{colorBalls(ib)} );
    
    text( 10+xi(ib,jfmax), yi(ib,jfmax), [num2str(Wfrac(ib,kfmax),'%4.2f')], 'FontSize', 12 );
    
end


hold off;

axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

axis 'equal';
axis 'off';
title( 'Highest entropy', 'FontSize', fs );

subplot(2,2,3);

hnmax

kfrange = hnmax:hnmax+nWindow-1;
jfmax = (hnmax-1)*sFrames+1;
jfrange = (kfrange-1)*sFrames+1;  % index to original frames;

hold on;
for ib=1:nBalls
    
    plot( xi(ib,jfrange), yi(ib,jfrange), 'o', 'Color', color{colorBalls(ib)} );
    
    text( 10+xi(ib,jfmax), yi(ib,jfmax), [num2str(Wfrac(ib,hnmax),'%4.2f')], 'FontSize', 12 );
    
end


hold off;

axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

axis 'equal';
axis 'off';
title( 'Highest fraction', 'FontSize', fs );

subplot(2,2,4);

lnmin

kfrange = lnmin:lnmin+nWindow-1;
jfmin = (lnmin-1)*sFrames+1;
jfrange = (kfrange-1)*sFrames+1;  % index to original frames;

hold on;
for ib=1:nBalls
    
    plot( xi(ib,jfrange), yi(ib,jfrange), 'o', 'Color', color{colorBalls(ib)} );
    
    text( 10+xi(ib,jfmin), yi(ib,jfmin), [num2str(Wfrac(ib,lnmin),'%4.2f')], 'FontSize', 12 );
    
end


hold off;

axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

axis 'equal';
axis 'off';
title('Lowest fraction', 'FontSize', fs );

print 'untitled3' -depsc2;

return;
