    
%% Load Data
Folder.Settings= 'Settings/'; % path to settings folder 

%ObserverName = 'events';

SettingsFile = input('Please enter the Settings datafile name: ', 's');

Settings= CExperimentalSettings(SettingsFile, Folder.Settings);

DataFile= input('Please enter the MOT datafile name: ', 's');

iTrial = input('Please enter the Trial Number: ', 's');

data = load(strcat('MOT-',DataFile,'/MOT-',DataFile,'-T',iTrial));

xi = data.MOT.xi;
yi = data.MOT.yi;
vi = data.MOT.vi;
oi = data.MOT.oi;
vix = data.MOT.vix;
viy = data.MOT.viy;
BallColors = data.MOT.BallColors;

Settings.MOT = data.MOT.Settings.Current.MOT;
Settings.Screen = data.MOT.Settings.Current.Screen;

%% Set Parameters
%FrameRate = floor2(floor(size(xi,2)/Settings.MOT.Duration),2);
FrameRate = Settings.Screen.ScreenMode.hz;

nBalls = 3;

allBalls = Settings.MOT.TotalN;

%% Select Balls
ib = 1;
for b=1:allBalls
   
    if isequal(BallColors(b,:),Settings.MOT.Target.AllRGB(Settings.Events.SelectedColor,:))
        
        selectedBalls(ib) = b;
    
        ib = ib + 1;
        
    end
    
end

distractor = selectedBalls(nBalls+1);
allDistractors = selectedBalls(nBalls+1:end);
selectedBalls = selectedBalls(1:nBalls);

maxvi = max(max(vi(selectedBalls,:)));
maxoi = max(max(oi(selectedBalls,:)));

nFrames = size(xi,2);

%% Area Events
for n=1:nFrames
   
    areaBalls(n) = tri_area([xi(selectedBalls(1),n) yi(selectedBalls(1),n)],[xi(selectedBalls(2),n) yi(selectedBalls(2),n)],[xi(selectedBalls(3),n) yi(selectedBalls(3),n)]);
    
end

maxArea = max(areaBalls);
eventsArea = zeros(1,Settings.MOT.Duration*2);
hs = 1;
for s=1:Settings.MOT.Duration
    
    startWindow = (1 + FrameRate*(s-1));
    endWindow = (FrameRate*(s-1) + (FrameRate/2));
    
    if sum(areaBalls(startWindow:endWindow) < Settings.Events.AreaThreshold*maxArea) ~= 0

        eventsArea(hs) = 1;

    end
    
    hs = hs + 1;
    
    startWindow = ((FrameRate/2) + 1 + FrameRate*(s-1));
    endWindow = (FrameRate*s);
    
    if sum(areaBalls(startWindow:endWindow) < Settings.Events.AreaThreshold*maxArea) ~= 0

        eventsArea(hs) = 1;

    end

    hs = hs + 1;
    
end

%% Displacement Events

displacement = zeros(nBalls,Settings.MOT.Duration*2);
eventsBallDistance = zeros(nBalls,Settings.MOT.Duration*2);

hs = 1;
for s=1:Settings.MOT.Duration
    
    for b=1:nBalls
        
        firstFrame = (1 + FrameRate*(s-1));
        lastFrame = (FrameRate*(s-1) + FrameRate/2);
        ball = selectedBalls(b);
        displacement(b,hs) = hypot(xi(ball,firstFrame) - xi(ball,lastFrame),yi(ball,firstFrame) - yi(ball,lastFrame));
        
        if displacement(b,hs) > Settings.Events.DisplacementThreshold*Settings.MOT.Area.Size(1)
            
            eventsBallDistance(b,hs) = 1;
            
        end
        
    end
    
    hs = hs + 1;
    
    for b=1:nBalls
        
        firstFrame = (FrameRate/2 + 1 + FrameRate*(s-1));
        lastFrame = (FrameRate*s);
        ball = selectedBalls(b);
        displacement(b,hs) = hypot(xi(ball,firstFrame) - xi(ball,lastFrame),yi(ball,firstFrame) - yi(ball,lastFrame));
        
        if displacement(b,hs) > Settings.Events.DisplacementThreshold*Settings.MOT.Area.Size(1)
            
            eventsBallDistance(b,hs) = 1;
            
        end
        
    end
    
    hs = hs + 1;
    
end

%% Distance Between All Balls

distanceAllBalls = zeros(allBalls,allBalls,nFrames);
for n=1:nFrames
    
    for b=1:allBalls
        
        distanceAllBalls(:,b,n) = hypot(xi(b,n) - xi(:,n),yi(b,n) - yi(:,n));
        
    end
    
end

%% Color Similarity Events

near = zeros(allBalls,nBalls,nFrames);
eventsBallColorMatch = zeros(nBalls,Settings.MOT.Duration*2);

for b=1:nBalls
    
   for n=1:nFrames
       
        near(:,b,n) = distanceAllBalls(:,selectedBalls(b),n) < Settings.Events.SameColorProximityThreshold;
      
        if (n == 1); hs = 1; elseif (mod(n,FrameRate/2) == 1); hs = hs + 1; end
 
        if ( (mod(n,FrameRate/2) == 0) && (n < FrameRate*Settings.MOT.Duration) )
        
          idx_alwaysNear = 1;
          
          for nn=1:size(near,1)
              
%               if ~any(near(nn,b,(n - FrameRate/2 + 1):n) - 1)
                allTimesNear = squeeze(near(nn,b,(n - FrameRate/2 + 1):n));
            
            if sum(allTimesNear) >= Settings.Events.SameColorProximityFrameThreshold
                  
                alwaysNear(idx_alwaysNear) = nn;
                
                idx_alwaysNear = idx_alwaysNear + 1;
                
            end
         
          end
          
          alwaysNear(alwaysNear == selectedBalls(b)) = [];
          
          if Settings.Events.SameColorTargetWithSameColorDistractor
              
              for g=1:length(selectedBalls)
                  
                  alwaysNear(alwaysNear == selectedBalls(g)) = [];
                  
              end
              
          end

          if ~isempty(alwaysNear)
             
              for nn=1:length(alwaysNear)

                if isequal(BallColors(alwaysNear(nn),:),BallColors(selectedBalls(b),:))

                    eventsBallColorMatch(b,hs) = 1;

                end
                
              end

          end
        
        end
        
   end
    
end

%% Mixing Events

if Settings.Events.MixingAllBalls; numBalls = allBalls; else numBalls = nBalls; end

windowSize = FrameRate/Settings.Events.MixingWindowsPerSecond;

if ~Settings.Events.MixingNoOverlappingWindows; step = 2; else step = 1; end

overlapOffset = windowSize/step;
   
near = zeros(allBalls,allBalls,nFrames);
sweep = zeros(allBalls,allBalls,Settings.MOT.Duration*Settings.Events.MixingWindowsPerSecond*step);
eventsBallMixing = zeros(allBalls,allBalls,Settings.MOT.Duration*Settings.Events.MixingWindowsPerSecond*step);
       
nearMixing(:,:,:) = distanceAllBalls(:,:,:) < Settings.Events.MixingProximityDistanceThreshold;
   
for b=1:numBalls
    
   if Settings.Events.MixingAllBalls; actualBall = b; else actualBall = selectedBalls(b); end
   
   for n=1:(nFrames - (1-Settings.Events.MixingNoOverlappingWindows)*overlapOffset) 
      
      nextWindowCondition = mod(n,(windowSize - (1-Settings.Events.MixingNoOverlappingWindows)*overlapOffset));
 
      if (n == 1); window = 1; elseif (nextWindowCondition == 1); window = window + 1; end
      
      if ( nextWindowCondition == 0) && (n <= FrameRate*Settings.MOT.Duration) 
      
          idx_alwaysNear = 1;
          
          alwaysNear = [];
          
          for nn=1:size(nearMixing,1)

              allTimesNear = squeeze(nearMixing(nn,actualBall,(n - (windowSize - (1-Settings.Events.MixingNoOverlappingWindows)*overlapOffset) + 1):(n + (1-Settings.Events.MixingNoOverlappingWindows)*overlapOffset)));
            
              if (sum(allTimesNear) >= Settings.Events.MixingProximityFrameThreshold) && (nn ~= actualBall)

                alwaysNear(idx_alwaysNear) = nn;
                
                idx_alwaysNear = idx_alwaysNear + 1;
                
              end
         
          end
          
          if ~isempty(alwaysNear)
             
              for nn=1:length(alwaysNear)
                  
                  iiFrame = 1;
                  angle = zeros(1,FrameRate/Settings.Events.MixingWindowsPerSecond);
                  
                  for iFrame=(n - (windowSize - (1-Settings.Events.MixingNoOverlappingWindows)*overlapOffset) + 1):(n + (1-Settings.Events.MixingNoOverlappingWindows)*overlapOffset)
                      
                      origin = [xi(actualBall,iFrame) yi(actualBall,iFrame)];
                      
                      relativeNear = [xi(alwaysNear(nn),iFrame) yi(alwaysNear(nn),iFrame)] - origin;
                      
                      angle(iiFrame) = radtodeg( atan ( relativeNear(2) / relativeNear(1) ) );
                      
                      iiFrame = iiFrame + 1;
                      
                  end
                  
                  sweep(alwaysNear(nn),actualBall,window) = sum( abs( diff( angle ) ) );
                  
                  if sweep(alwaysNear(nn),actualBall,window) > Settings.Events.MixingSweepThreshold
                  
                      eventsBallMixing(alwaysNear(nn),actualBall,window) = 1;
                      
                  end
                  
              end
              
          end
          
      end
      
   end
   
end

clear data;
clear n;
clear nn;
clear b;
clear hs;

save('eventsBalls.mat');
    