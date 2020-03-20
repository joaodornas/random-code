
function MOT = GenerateBraunianWalkTrajectories

%% items moving randomly in square field, constantly changing their relative positions
  %  tracking these items engages attention effectively

  Folder.Settings= 'Settings/'; % path to settings folder
  
  ObserverName = 'trajectories';
  
  Settings= CExperimentalSettings(ObserverName, Folder.Settings);
  
  PresentationDuration_Sec = Settings.MOT.Duration + Settings.MOT.Delay;
  Oversample = Settings.MOT.Oversample;
  FrameDuration= 1./Settings.Screen.ScreenMode.hz;
  ItemsN = Settings.MOT.TotalN;
  
  %% array that describe a complete sequence
  FramesN = ceil(PresentationDuration_Sec*Oversample/FrameDuration);  % number of iteration steps
  vi = nan(ItemsN, FramesN);   % speed of item iItem at iternation ti (units of Settings.MOT.v0)
  ai = nan(ItemsN, FramesN);   % direction of item iItem at iteration ti (-pi to pi)
  oi = nan(ItemsN, FramesN);   % angular velocity of item iItem at iteration ti (-pi to pi)
  xi = nan(ItemsN, FramesN);   % x position of item iItem
  yi = nan(ItemsN, FramesN);   % y position of item iItem

  vixx = nan(ItemsN, FramesN);
  viyy = nan(ItemsN, FramesN);
  
  %% generating initial position
  %xi(1,1) = round( (Settings.MOT.Area.Size(1)-2*Settings.MOT.RepulsionRadius)*(rand - 0.5) );
  %yi(1,1) = round( (Settings.MOT.Area.Size(2)-2*Settings.MOT.RepulsionRadius)*(rand - 0.5) );
  xi(1,1) = (Settings.MOT.Area.Size(1)-2*Settings.MOT.RepulsionRadius)*(rand - 0.5) ;
  yi(1,1) = (Settings.MOT.Area.Size(2)-2*Settings.MOT.RepulsionRadius)*(rand - 0.5) ;
  iItem = 2;
  while (iItem<=ItemsN),
    % position of the next item
    %xi(iItem,1) = round( (Settings.MOT.Area.Size(1)-2*Settings.MOT.RepulsionRadius)*(rand - 0.5) );
    %yi(iItem,1) = round( (Settings.MOT.Area.Size(2)-2*Settings.MOT.RepulsionRadius)*(rand - 0.5) );
    xi(iItem,1) = (Settings.MOT.Area.Size(1)-2*Settings.MOT.RepulsionRadius)*(rand - 0.5) ;
    yi(iItem,1) = (Settings.MOT.Area.Size(2)-2*Settings.MOT.RepulsionRadius)*(rand - 0.5) ;

    %% checking that it is sufficiently far away from all other items
    D= hypot(xi(iItem, 1)-xi(1:(iItem-1), 1), yi(iItem, 1)-yi(1:(iItem-1), 1));
    if min(D)>Settings.MOT.RepulsionRadius
      iItem=iItem+1;
    end
  end

  %% generating speed
  vi(:, 1)= Settings.MOT.v0;
  ai(:, 1)= 2*pi*rand(ItemsN, 1)-pi;
  oi(:, 1)= 0;
  dvx = vi(:,1).*cos(ai(:,1));
  dvy = vi(:,1).*sin(ai(:,1));


  %% iterate trajectories
  for ti=2:FramesN,
    %% initial linear velocities and angular velocities (vir and oir)
    vix = vi(:,ti-1).*cos(ai(:,ti-1));
    viy = vi(:,ti-1).*sin(ai(:,ti-1));
    vir = sqrt( vix.^2 + viy.^2 );
    oir = oi(:,ti-1);


    %% random increments to velocities (Brownian motion accelerations)
    % OU random increments to radial velocity, ss at Settings.MOT.v0
    dvir = -(vir-Settings.MOT.v0) * FrameDuration + Settings.MOT.sigma_v * sqrt(FrameDuration) * (2*(rand(ItemsN,1)<0.5)-1);
    % OU random increments to angular velocity, ss at 0
    doir = -(oir-Settings.MOT.o0) * FrameDuration + Settings.MOT.sigma_o * sqrt(FrameDuration) * (2*(rand(ItemsN,1)<0.5)-1);

    % new linear velocities after random acceleration
    vix  = (vir + dvir).*cos(ai(:,ti-1));
    viy  = (vir + dvir).*sin(ai(:,ti-1));

    % new angular velocities after random acceleration
    oir = oir + doir;

    % rotate linear velocities according to new angular velocity
    soi = sin(oir*FrameDuration);
    coi = cos(oir*FrameDuration);

    vixp = coi.*vix - soi.*viy;
    viyp = soi.*vix + coi.*viy;

    vix = vixp;
    viy = viyp;


    %% boundary repulsion 
    r_left   = xi(:,ti-1)+Settings.MOT.Area.Size(1)/2;
    r_right  = xi(:,ti-1)-Settings.MOT.Area.Size(1)/2;
    r_top    = yi(:,ti-1)+Settings.MOT.Area.Size(2)/2;
    r_bottom = yi(:,ti-1)-Settings.MOT.Area.Size(2)/2;

    epsilon  = 1;

    kk = find( r_left < Settings.MOT.RepulsionRadius );
    if ~isempty(kk)
         vix(kk) = vix(kk) + FrameDuration * Settings.MOT.RepulsionK ./ max(epsilon, r_left(kk)) ; 
    end

    kk = find( r_right > -Settings.MOT.RepulsionRadius );
    if ~isempty(kk)
         vix(kk) = vix(kk) + FrameDuration * Settings.MOT.RepulsionK ./ min(-epsilon, r_right(kk)) ; 
    end 

    kk = find( r_top < Settings.MOT.RepulsionRadius );
    if ~isempty(kk)
         viy(kk) = viy(kk) + FrameDuration * Settings.MOT.RepulsionK ./ max(epsilon, r_top(kk)) ; 
    end

    kk = find( r_bottom > -Settings.MOT.RepulsionRadius );
    if ~isempty(kk)
         viy(kk) = viy(kk) + FrameDuration * Settings.MOT.RepulsionK ./ min(-epsilon, r_bottom(kk)); 
    end 

    %% mutual repulsion 
    for iItem=1:ItemsN
      dx = xi(:,ti-1) - xi(iItem,ti-1);
      dy = yi(:,ti-1) - yi(iItem,ti-1);
      r_near = sqrt(dx.*dx + dy.*dy);
      kk = find(r_near > 0 & r_near < Settings.MOT.RepulsionRadius);
        if ~isempty(kk)
          dv = FrameDuration * Settings.MOT.RepulsionK ./ max( epsilon, r_near(kk) );

          dxn = dx(kk)./ r_near(kk);
          dyn = dy(kk)./ r_near(kk);

          vix(kk) = vix(kk) + dv .* dxn;
          viy(kk) = viy(kk) + dv .* dyn;
        end
    end

      %% update position and velocity

      xi(:,ti) = xi(:,ti-1) + vix * FrameDuration;
      yi(:,ti) = yi(:,ti-1) + viy * FrameDuration;

      vi(:,ti) = sqrt( vix.^2 + viy.^2 );
      vixx(:,ti) = vix;
      viyy(:,ti) = viy;
      ai(:,ti) = atan2( viy, vix ); 
      oi(:,ti) = oir;
  end

  xi= xi(:, 1:Oversample:end);
  yi= yi(:, 1:Oversample:end);
  
  vixx = vixx(:,1:Oversample:end);
  viyy = viyy(:,1:Oversample:end);
  
  BallColors = RNDColors(Settings);
  
    %% Randomise colors
    function BallColors = RNDColors(Settings)
         %Generate different colors
     NumOfColors=Settings.MOT.Target.NumOfDifferentColors;% Number of color types in run
     Colors=Settings.MOT.Target.AllRGB;% All possible colors
     NumOfAllColors=size(Colors,1);
     ColorsForRun=Colors(randperm(NumOfAllColors,NumOfColors),:);
     NumOfBallsPerColor=Settings.MOT.Target.NumOfBallsPerColor;
     BallColors=repmat(ColorsForRun,NumOfBallsPerColor,1);
     BallColors=BallColors(randperm(size(BallColors,1)),:);
     % cues should be a same color
     
    end

    MOT.xi = xi;
    MOT.yi = yi;
    MOT.vi = vi;
    MOT.vix = vixx;
    MOT.viy = viyy;
    MOT.ai = ai;
    MOT.oi = oi;
    MOT.BallColors = BallColors;
    MOT.Settings.MOT = Settings.MOT;
    MOT.Settings.Screen = Settings.Screen;

end