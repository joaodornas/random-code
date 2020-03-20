
function AnalyzeTrajectoriesV1(saveAnalyze,datafilename,saveAnalyzeDataFileName)

%%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%

disp('Load data file');

load(datafilename);
    
%MOT = GenerateBraunianWalkTrajectories;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DEFINE PLOTS %%%%%%%%%%%%%%%%%%%%%

plotTrajectories = false;
plotColorFraction = false;
checkResults = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SET PARAMETERS %%%%%%%%%%%%%%%%%%%%%

% display size and repulsion radius

fs = 12;

sizeX = 400;
sizeY = 400;
repR = 50;

% select part of sequence to be analyzed

start  = 61;               % skip first second
finish = length(MOT.xi);   % 

xi = MOT.xi(:,start:finish);
yi = MOT.yi(:,start:finish);

size(xi);
size(yi);

% sampling of frames and evaluation window

nFrames = size(xi,2);
sFrames = 3;   % sampling of frames
nBalls = size(xi,1);

nWindow = 18/sFrames;  % number of sampled frames per window

hWindow = ceil( nWindow / 2);   % nWindow is even, hWindow is half ...
nWindow = 2*hWindow;

% neighborhood for proximity evaluation

Nneighbor = min( 12, nBalls );  % number of neighbors to be considered

Nsigma = 100;    % size of Gaussian neighborhood

Ns2 = 2*Nsigma*Nsigma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PLOT TRAJECTORIES %%%%%%%%%%%%%%%%%%%%%

color{1} = 'r';
color{2} = 'b';
color{3} = 'g';
color{4} = 'k';
color{5} = 'y';

BallColors = MOT.BallColors;
AllRGB = MOT.Settings.MOT.Target.AllRGB;

for i=1:size(AllRGB,1)
    
    k = 1;
    
    thiscolor = AllRGB(i,:);
    
    for j=1:size(BallColors,1)
       
        ballcolor = BallColors(j,:);
        
        if isequal(thiscolor,ballcolor)
            
            ballsPerColor(i,k) = j;
            
            k = k + 1;
            
        end
        
    end
    
end

colorBalls = zeros(1,size(ballsPerColor,2));

for j=1:size(ballsPerColor,1)
    
    balls = ballsPerColor(j,:);
    
    for b=balls

        colorBalls(b) = j;

    end
    
end
 
if plotTrajectories
    
    figure;

    hold on;
    plot( [-0.5 0.5 0.5 -0.5 -0.5] * sizeX, [-0.5 -0.5 0.5 0.5 -0.5]*sizeY, 'k', 'LineWidth', 2 );
    plot( [-0.5 0.5 0.5 -0.5 -0.5] * (sizeX-2*repR), [-0.5 -0.5 0.5 0.5 -0.5]*(sizeY-2*repR), 'k--', 'LineWidth', 2 );

    for b=1:length(balls)

       plot(xi(balls(b),1:10:end),yi(balls(b),1:10:end),color{colorBalls(b)});

    end
    hold off;

    axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

    axis 'equal';
    axis 'off';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% COMPUTE PARAMETERS  %%%%%%%%%%%%%%%%%%%%%

disp('Compute parameters');

%% compute pairwise distances to evaluate the proximity of various qualities
%
% qualities include same color, similar direction of motion, bouncing motion (consistent changes in
% relative direction) and shifting quadrant (consistent changes in quadrant)
%
% note that same color and similar motion are INSTANTANEOUS, whereas
% bounces (changes in relative motion) and shifts (changes in relative
% quadrant) are NON-INSTANTANEOUS
%
%
% in the future, we may consider introducing changes of COLOR as yet another class 
% NON-INSTANTANEOUS EVENT
%

% end results to be saved and evaluated over all sampled frames

Wcolor  = nan(nBalls, ceil(nFrames/sFrames) );  % weighted average of same COLOR in proximity
Pcolor =  nan(nBalls, ceil(nFrames/sFrames) );  % probability density of above 

Wmotion = nan(nBalls, ceil(nFrames/sFrames) );  % weighted average of similar MOTION in proximity
Pmotion = nan(nBalls, ceil(nFrames/sFrames) );  % probability density of above

Wbounce = nan(nBalls, ceil(nFrames/sFrames) );  % weighted average of BOUNCING balls in proximity
Pbounce = nan(nBalls, ceil(nFrames/sFrames) );  % probability density of above fraction

Wshift  = nan(nBalls, ceil(nFrames/sFrames) );  % weighted average of SHIFTING balls in proximity
Pshift  = nan(nBalls, ceil(nFrames/sFrames) );  % probability density of above 

% end result to be saved over all sampled frames

Ecolor = nan(ceil(nFrames/sFrames),1);    % color entropy, low value, high value
Lcolor = nan(ceil(nFrames/sFrames),1);
Hcolor = nan(ceil(nFrames/sFrames),1);

Emotion = nan(ceil(nFrames/sFrames),1);   % motion entropy, low value, high value
Lmotion = nan(ceil(nFrames/sFrames),1);
Hmotion = nan(ceil(nFrames/sFrames),1);

Ebounce = nan(ceil(nFrames/sFrames),1);   % motion change (bounce) entropy, low value, high value
Lbounce = nan(ceil(nFrames/sFrames),1);
Hbounce = nan(ceil(nFrames/sFrames),1);

Eshift = nan(ceil(nFrames/sFrames),1);   % quadrant change (shift) entropy, low value, high value
Lshift = nan(ceil(nFrames/sFrames),1);
Hshift = nan(ceil(nFrames/sFrames),1);

%% iterate over all sampled frames

% intermediate results to be saved and evaluated over nWindow sampled
% frames

Nexp2       = nan(nBalls,nBalls,nWindow);   % pairwise Gaussian factor, measure of proximity
Ncolor2     = nan(nBalls,nBalls,nWindow);  % pairwise realtive color, same=1, different=0
Nmotion2    = nan(nBalls,nBalls,nWindow);  % pairwise relative motion, velocity dot product, similar=positive, dissimilar=negative
Nquadrant2  = nan(nBalls,nBalls,nWindow);  % pairwise index of relative quadrant
    
dv12all = nan(nBalls,nBalls,ceil(nFrames/sFrames));

for kf=1:ceil(nFrames/sFrames)  % index to sampled frames
    
    jf  = (kf-1)*sFrames+1;  % index to original frames
    
    xcur = xi(:,jf);           % current position
    ycur = yi(:,jf);
    
    if kf>1                         % starting with second frame, compute relative color, relative
                                    % velocity, and relative quadrant, for each sampled frame
                                    
        vxcur = xcur - xpvs;    % current velocity
        vycur = ycur - ypvs;
        
        for ib1=1:nBalls-1  % compute pairwise distances
            
            x1  = xcur(ib1);   % current position
            y1  = ycur(ib1);
            
            vx1 = vxcur(ib1);    % current velocity
            vy1 = vycur(ib1);
            v1 = sqrt( vx1^2 + vy1^2 );
            
            clr1 = colorBalls(ib1);  % current color (may change over time in future)
            
            for ib2=ib1+1:nBalls
                
                x2  = xcur(ib2);  % current position
                y2  = ycur(ib2);
                
                vx2 = vxcur(ib2);    % current velocity
                vy2 = vycur(ib2);
                v2 = sqrt( vx2^2 + vy2^2 );
                
                clr12 = (clr1==colorBalls(ib2));   % relative color
                
                dx12 = x2-x1;    % relative position
                dy12 = y2-y1;
                
                dvx12 = vx2 - vx1;
                dvy12 = vy2 - vy1;
                
                alpha12 = acos( dot([vx1 vy1],[vx2 vy2])/(v1*v2) );
                
                %dv12 = vx1*vx2 + vy1*vy2;   % dot product of velocities (positve or negative)
                %dv12 = sqrt( dvx12^2 + dvy12^2 ) / sqrt( v1*v2 ) ;
                dv12 = v1*v2 * cos( alpha12 );
                
                dv12all(ib1,ib2,kf) = dv12;
                
                dq12 = 1 + (dx12>=0) + 2*(dy12>=0);   % index of relative quadrant
                dq21 = 1 + (-dx12>=0) + 2*(-dy12>=0);   %
                
                %                   |
                %                   |
                %          dq12 = 3 |   dq12 = 4
                %                   |
                %        ___________|______________
                %                   |
                %                   |
                %          dq12 = 1 |   dq12 = 2
                %                   |
                
                
                expfac = exp( -(dx12*dx12 + dy12*dy12) / Ns2 );    % save exponential factors
                Nexp2(ib1, ib2, 1)      = expfac;
                Nexp2(ib2, ib1, 1)      = expfac;
                
                Ncolor2(ib1, ib2, 1)    = clr12;                               % relative color: same=1, different=0
                Ncolor2(ib2, ib1, 1)    = clr12;
                
                Nmotion2(ib1, ib2, 1)   = dv12;                                % relative motion: vel dot product
                Nmotion2(ib2, ib1, 1)   = dv12;
                
                Nquadrant2(ib1, ib2, 1) = dq12;                                % save index to relative quadrant
                Nquadrant2(ib2, ib1, 1) = dq21;
                 
            end         
        end
    end
    
    xpvs = xcur;                    % save for velocity computation
    ypvs = ycur;
    
               %               
               % <---------- nWindow frames ---------->
               %
               %  W1  W2 ...                        Wn     instantaneous proximity measure
               %
               %       rel. color: 0 or 1
               %      rel. motion: -val to +val (instantaneous)
               %       degree of bounce
               %       degree of shift

    if kf >= 2      % starting with 2nd frame, evaluate instantaneous qualities (color, motion)
        
        for ib1=1:nBalls; 
            
            nexp2   = squeeze(Nexp2(ib1,:,1));      % exponential weight of neigbhors (one frame)
            
            clr2     = squeeze(Ncolor2(ib1,:,1));   % get relative color of neighbors
            
            %vel2     = squeeze(Nmotion2(ib1,:,1));  % get relative motion of neighbors
            
            Wcolor( ib1, kf )  = nansum( nexp2 .* clr2 ) / nansum( nexp2 );  % weighted sum of relative color
            
            %Wmotion( ib1, kf ) = nansum( nexp2 .* vel2 ) / nansum( nexp2 );  % weighted sum of relative motion
            
        end
       
    end
    
    if kf >= nWindow+1      % starting with nWindow + 1st frame, evaluate non-instantaneous qualities (bounce, shift)
        
        for ib1=1:nBalls; 
            
            nexp2  = squeeze(Nexp2(ib1,:,hWindow));      % exponential weight of neigbhors (middle frame)
            
            vel2pre = nanmean(Nmotion2(ib1,:,1:hWindow),3); % relative velocity in first half
            vel2pos = nanmean(Nmotion2(ib1,:,hWindow+1:nWindow),3);  % relative velocity in second half
            
            bounce2 = -vel2pre .* vel2pos;   % measure of bounce (positive to negative)
            
            %% it should be the opposite, the first half is the future, the last half is the past
            qua2pre = nanmean(Nquadrant2(ib1,:,1:hWindow),3); % relative quadrant in first half
            qua2pos = nanmean(Nquadrant2(ib1,:,hWindow+1:nWindow),3);  % relative quadrant in second half
            
            shift2  = abs((qua2pos - qua2pre)/2) + abs( mod(qua2pos,2) - mod(qua2pre,2)); % measure of shifts (zero to two)
            
            vel2     = nansum( squeeze(Nmotion2(ib1,:,1:nWindow)), 2 );  % get relative motion of neighbors
            
            %Wbounce( ib1, kf-hWindow ) = nansum( nexp2 .* bounce2) / nansum( nexp2 ); % weighted sum of bounce
            Wbounce( ib1, kf-nWindow ) = nansum( nexp2 .* bounce2) / nansum( nexp2 ); % weighted sum of bounce
               
            %Wshift( ib1, kf-hWindow )  = nansum( nexp2 .* shift2) / nansum( nexp2 ); % weighted sum of shift 
            Wshift( ib1, kf-nWindow )  = nansum( nexp2 .* shift2) / nansum( nexp2 ); % weighted sum of shift 
            
            Wmotion( ib1, kf-nWindow ) = nansum( nexp2 .* vel2' ) / nansum( nexp2 );  % weighted sum of relative motion
            
        end
            
    end
    
    Nexp2(:,:,2:end)      = Nexp2(:,:,1:end-1);               % save for Window evaluation
    Ncolor2(:,:,2:end)    = Ncolor2(:,:,1:end-1);             
    Nmotion2(:,:,2:end)   = Nmotion2(:,:,1:end-1);
    Nquadrant2(:,:,2:end) = Nquadrant2(:,:,1:end-1);
    
    if mod( kf, 1000) == 0
        disp(strcat('kf:',int2str(kf)));
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% COMPUTE DISTRIBUTIONS AND PROBABILITES %%%%%%%%%%%%%%%%%%%%%

disp('Compute Distributions and Probabilites');

%% show distributions and establish histograms

[Fcolor,  Cbins]  = ShowDistributions( Wcolor, 'Same color' );
[Fmotion, Mbins]  = ShowDistributions( Wmotion, 'Same motion' );
[Fbounce, Bbins]  = ShowDistributions( Wbounce, 'Motion change' );
[Fshift,  Sbins]  = ShowDistributions( Wshift, 'Quadrant change' );

%% save probabilities

Pcolor  = Weigth2Probability( Wcolor,  Fcolor,  Cbins );
Pmotion = Weigth2Probability( Wmotion, Fmotion, Mbins );
Pbounce = Weigth2Probability( Wbounce, Fbounce, Bbins );
Pshift  = Weigth2Probability( Wshift,  Fshift,  Sbins );

%% check the computation of color fraction

if plotColorFraction
    
    figure;

    hold on;
    for ib=1:nBalls

        plot( xi(ib,2), yi(ib,2), 'o', 'Color', color{colorBalls(ib)} );

        text( 10+xi(ib,2), yi(ib,2), [num2str(Wcolor(ib,2),'%4.2f')], 'FontSize', fs );

    end

    hold off;

    axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

    axis 'equal';
    axis 'off';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FIND EXTREMES %%%%%%%%%%%%%%%%%%%%%

disp('Find Extremes');

%% we have now computed several proximity-weighted qualities, as well as the associate probabilities

%  Wcolor, Pcolor       proximity of balls of same color

%  Wmotion, Pmotion     proximity of balls with similar motion

%  Wbounce, Pbounce     proximity of balls changing relative motion (bouncing)

%  Wshift, Pshift       proximity of balls changing relative quadrant (shifting)

% we will evaluate the results for 'sliding windows', each containing nWindow frames
%  
% exceptionally high or low values of qualities will be evaluated for the middle frame of the window

% high ensemble entropy indicates that proximity measures are highly heterogeneous
% low ensemble entropy indicates highly homogeneous proximity measures

% in the case of instantaneous qualities (color, motion) ensemble entropy will be computed over all nWindow frames

% in the case non-instantaneous qualities (shift, bounce) ensemble entropy will be computed over the middle frame

for kf=2:ceil(nFrames/sFrames)-nWindow-nWindow   % compute entropy from probabilities P, compute high and low values from weights W
    
    mf = kf+hWindow;      % index to middle frame
    %rf = kf:kf+nWindow-1;   % index to sampled frames
    rf = kf:kf+nWindow;   % index to sampled frames
    
    Pselect = Pcolor(:,rf);   % evaluate color (entire range rf)
    Wselect = Wcolor(:,rf);
    Ecolor(kf) = nansum( -Pselect(:) .* log( Pselect(:) ) );  % evaluate entropy from probabilities ...
    meanfrac = nanmean( Wselect, 2 ) ;  % mean fraction over all frames, separately for each ball, sorted
    idx_min = find(meanfrac == min(meanfrac));
    idx_max = find(meanfrac == max(meanfrac));
    Lcolor(kf) = meanfrac(idx_min(1));   % lowest value
    Hcolor(kf) = meanfrac(idx_max(end)); % highest values
    LcolorBall(kf) = idx_min(1);
    HcolorBall(kf) = idx_max(end);
    
    Pselect = Pmotion(:,rf);  % evaluate motion (entire range rf)
    Wselect = Wmotion(:,rf); 
    Emotion(kf) = nansum( -Pselect(:) .* log( Pselect(:) ) ); 
    meanfrac = nanmean( Wselect, 2 ) ;  % mean fraction over all frames, separately for each ball, sorted
    idx_min = find(meanfrac == min(meanfrac));
    idx_max = find(meanfrac == max(meanfrac));
    Lmotion(kf) = meanfrac(idx_min(1));   % lowest value
    Hmotion(kf) = meanfrac(idx_max(end)); % highest values
    LmotionBall(kf) = idx_min(1);
    HmotionBall(kf) = idx_max(end);
    
    Pselect = Pbounce(:,mf);   % evaluate bounce  (middle frame mf only)
    Wselect = Wbounce(:,mf);
    Ebounce(kf) = nansum( -Pselect(:) .* log( Pselect(:) ) ); 
    meanfrac = nanmean( Wselect, 2 ) ;  % mean fraction over all frames, separately for each ball, sorted
    idx_min = find(meanfrac == min(meanfrac));
    idx_max = find(meanfrac == max(meanfrac));
    Lbounce(kf) = meanfrac(idx_min(1));   % lowest value
    Hbounce(kf) = meanfrac(idx_max(end)); % highest values
    LbounceBall(kf) = idx_min(1);
    HbounceBall(kf) = idx_max(end);
    
    Pselect = Pshift(:,mf);  % evaluate shift (middle frame mf only)
    Wselect = Wshift(:,mf);
    Eshift(kf) = nansum( -Pselect(:) .* log( Pselect(:) ) ); 
    meanfrac = nanmean( Wselect, 2 ) ;  % mean fraction over all frames, separately for each ball, sorted
    idx_min = find(meanfrac == min(meanfrac));
    idx_max = find(meanfrac == max(meanfrac));
    Lshift(kf) = meanfrac(idx_min(1));   % lowest value
    Hshift(kf) = meanfrac(idx_max(end)); % highest values
    LshiftBall(kf) = idx_min(1);
    HshiftBall(kf) = idx_max(end);
    
end

%% find extremes in different analyses

[kECmin, kECmax, kVCmin, kVCmax, kVCballmin, kVCballmax] = FindExtremes( Ecolor, Lcolor, Hcolor, LcolorBall, HcolorBall, 'Color' );

[kEMmin, kEMmax, kVMmin, kVMmax, kVMballmin, kVMballmax] = FindExtremes( Emotion, Lmotion, Hmotion, LmotionBall, HmotionBall, 'Motion' );

[kEBmin, kEBmax, kVBmin, kVBmax, kVBballmin, kVBballmax] = FindExtremes( Ebounce, Lbounce, Hbounce, LbounceBall, HbounceBall, 'Bounce' );

[kEQmin, kEQmax, kVQmin, kVQmax, kVQballmin, kVQballmax] = FindExtremes( Eshift, Lshift, Hshift, LshiftBall, HshiftBall, 'Quadrant' );


if checkResults
    
    %% check the results
    CheckResults( kECmin, kECmax, kVCmin, kVCmax, Wcolor, Ecolor, sFrames, nWindow, nBalls, xi, yi, color, colorBalls, sizeX, sizeY, 'Color'   )
    CheckResults( kEMmin, kEMmax, kVMmin, kVMmax, Wmotion, Emotion, sFrames, nWindow, nBalls, xi, yi, color, colorBalls, sizeX, sizeY, 'Motion'  )
    CheckResults( kEBmin, kEBmax, kVBmin, kVBmax, Wbounce, Ebounce, sFrames, nWindow, nBalls, xi, yi, color, colorBalls, sizeX, sizeY, 'Bounce'  )
    CheckResults( kEQmin, kEQmax, kVQmin, kVQmax, Wshift, Eshift, sFrames, nWindow, nBalls, xi, yi, color, colorBalls, sizeX, sizeY, 'Quadrant'  )

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SAVE RESULTS %%%%%%%%%%%%%%%%%%%%

if saveAnalyze
    
    clear MOT
    clear xi
    clear yi
    
    save(saveAnalyzeDataFileName);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FUNCTIONS DEFINITIONS %%%%%%%%%%%%%%%%%%%%%

function [Pquality, Bins] = ShowDistributions( Wquality, titlestring )

    plotDistributions = false;
    
    fs = 12;

    [nquality,xquality] = hist(Wquality(:),21);    % bin counts, bin positions

    fquality = nquality / sum(nquality);           % count fractions

    dquality = fquality / (xquality(2) - xquality(1));   % count density

    if plotDistributions
        
        figure;

        hold on;
        bar( xquality, dquality  ); 
        hold off;
        xlabel('W fraction', 'FontSize', fs );
        ylabel('density', 'FontSize', fs );

        suptitle( titlestring );

    end
    
    Pquality = fquality;
    Bins     = xquality;

end

function Wprob = Weigth2Probability( Wfrac, Pfrac, Bins )

    % replace exponential weight fractions with probability fractions, bin by bin

    Wprob = nan(size(Wfrac));

    kk = find( Wfrac <= 0.5*(Bins(1)+Bins(2)) );       % first bin
    Wprob(kk) = Pfrac(1);

    for ii=2:length(Bins)-1
        elo = 0.5*(Bins(ii-1)+Bins(ii));               % intermediate bins
        ehi = 0.5*(Bins(ii)+Bins(ii+1));
        kk = find( elo < Wfrac & Wfrac <= ehi );
        Wprob(kk) = Pfrac(ii);
    end

    kk = find( Wfrac > 0.5*(Bins(end-1)+Bins(end)) );  % last bin
    Wprob(kk) = Pfrac(end);

end

function [kEmin, kEmax, kVmin, kVmax, kVballmin, kVballmax] = FindExtremes( Etype, Ltype, Htype, LtypeBall, HtypeBall, titlestring )

    plotExtremes = false;

    fs = 12;

    dummy = find( Etype == nanmin(Etype ) );   % smallest entropy over all sampled frames
    kEmin = dummy(1);

    dummy = find( Etype == nanmax(Etype ) );   % highest entropy over all sampled frames
    kEmax = dummy(1);

    dummy = find( Ltype == nanmin(Ltype ) );     % highest fraction over all sampled frames
    kVmin = dummy(1);

    dummy = find( Htype == nanmax(Htype ) );     % lowest fraction over all sampled frames
    kVmax = dummy(1);
    
    kVballmin = LtypeBall(kVmin);
    kVballmax = HtypeBall(kVmax);

    if plotExtremes
        
        figure;
        subplot(2,1,1);
        hold on;
        plot( 1:length(Etype), Etype, 'r.' );
        plot( kEmin, Etype(kEmin), 'ko', 'MarkerSize', 10 );
        plot( kEmax, Etype(kEmax), 'ko', 'MarkerSize', 10 );
        hold off;
        xlabel('frames', 'FontSize', fs );
        ylabel('Entropy', 'FontSize', fs );

        subplot(2,1,2);
        hold on;
        plot( 1:length(Ltype), Ltype, 'r.' );
        plot( 1:length(Htype), Htype, 'b.' );
        plot( kVmin, Ltype(kVmin), 'ko', 'MarkerSize', 10 );
        plot( kVmax, Htype(kVmax), 'ko', 'MarkerSize', 10 );
        hold off;
        xlabel('frames', 'FontSize', fs );
        ylabel('W fraction', 'FontSize', fs );

        suptitle( titlestring );
    
    end

end

function CheckResults( kemin, kemax, kvmin, kvmax, Wfrac, Hentropy, sFrames, nWindow, nBalls, xi, yi, color, colorBalls, sizeX, sizeY, titlestring  )

    fs = 12;

    figure;

    subplot(2,2,1);

    kerange = kemin:kemin+nWindow;  % index to sampled frames;
    %jemin = (kemin-1)*sFrames+1;
    jerange = (kerange-1)*sFrames+1;  % index to original frames;

    hold on;
    for ib=1:nBalls

        plot( xi(ib,jerange), yi(ib,jerange), 'o', 'Color', color{colorBalls(ib)} ); 

    end

    hold off;

    axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

    axis 'equal';
    axis 'off';
    title( strcat('Lowest entropy:',titlestring,':',num2str(Hentropy(kemin),'%4.2f')), 'FontSize', fs );

    subplot(2,2,2);

    kerange = kemax:kemax+nWindow;  % index to sampled frames;
    jemax = (kemax-1)*sFrames+1;
    jerange = (kerange-1)*sFrames+1;  % index to original frames;

    hold on;
    for ib=1:nBalls

        plot( xi(ib,jerange), yi(ib,jerange), 'o', 'Color', color{colorBalls(ib)}, 'MarkerSize', 8 );
        plot( xi(ib,jerange(1)), yi(ib,jerange(1)), 'x', 'Color', color{3}, 'MarkerSize', 3 );

    end


    hold off;

    axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

    axis 'equal';
    axis 'off';
    title( strcat('Highest entropy:',titlestring,':',[num2str(Hentropy(kemax),'%4.2f')]), 'FontSize', fs );

    subplot(2,2,3);

    kvrange = kvmax:kvmax+nWindow;
    jvmax = (kvmax-1)*sFrames+1;
    jvrange = (kvrange-1)*sFrames+1;  % index to original frames;

    hold on;
    for ib=1:nBalls

        plot( xi(ib,jvrange), yi(ib,jvrange), 'o', 'Color', color{colorBalls(ib)} );

        text( 10+xi(ib,jvmax), yi(ib,jvmax), [num2str(Wfrac(ib,kvmax),'%4.2f')], 'FontSize', fs );

    end


    hold off;

    axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

    axis 'equal';
    axis 'off';
    title( strcat('Highest fraction:',titlestring), 'FontSize', fs );

    subplot(2,2,4);

    kvrange = kvmin:kvmin+nWindow;
    jvmin = (kvmin-1)*sFrames+1;
    jvrange = (kvrange-1)*sFrames+1;  % index to original frames;

    hold on;
    for ib=1:nBalls

        plot( xi(ib,jvrange), yi(ib,jvrange), 'o', 'Color', color{colorBalls(ib)} );

        text( 10+xi(ib,jvmin), yi(ib,jvmin), [num2str(Wfrac(ib,kvmin),'%4.2f')], 'FontSize', fs );

    end

    hold off;

    axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

    axis 'equal';
    axis 'off';
    title(strcat('Lowest fraction:',titlestring), 'FontSize', fs );
    
end