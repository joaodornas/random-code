
DataFile= input('Please enter the Parameters file name: ', 's');

load(DataFile);

Folder.Settings = 'Settings/';

SettingsFile = input('Please enter the Events Settings file name: ', 's');

SettingsEvents = CExperimentalSettings(SettingsFile, Folder.Settings);

Settings.Events = SettingsEvents.Events;
Settings.Screen = SettingsEvents.Screen;
Settings.Keyboard = SettingsEvents.Keyboard;

nBalls = 3;

%% Area Events

eventsArea = zeros(1,nFrames);
areaComb = 0;
for i=1:size(allCombSelectedBalls,1)
    
    if sum(ismember(selectedBalls(1:nBalls),allCombSelectedBalls(20,:))) == nBalls
        
        areaComb = i;
        break
        
    end
    
end

for s=1:nFrames

    if areaBalls(areaComb,s) < Settings.Events.AreaThreshold*maxArea(areaComb) 

        eventsArea(s) = 1;

    end
    
end

%% Displacement Events

eventsBallDisplacement = zeros(nBalls,nFrames);

for s=1:(nFrames-FrameRate*Settings.Events.DisplacementWindowInSeconds)
    
    for b=1:nBalls
        
        if  maxDisplacement(selectedBalls(b),s) > Settings.Events.DisplacementThreshold*Settings.MOT.Area.Size(1)
            
            eventsBallDisplacement(b,s) = 1;
            
        end
        
    end
    
end


%% Color Similarity Events

eventsBallColorMatch = zeros(nBalls,nFrames);
nMatches = zeros(nBalls,nFrames);

for b=1:nBalls
    
   for n=1:(nFrames - Settings.Events.SameColorWindowInSeconds*FrameRate)
       
       for bb=1:allBalls
           
           if ballColorMatch(selectedBalls(b),bb,n)
               
               nMatches(b,n) = nMatches(b,n) + sum(ballColorMatch(selectedBalls(b),bb,n:n+Settings.Events.SameColorWindowInSeconds*FrameRate));
               
           end
           
       end
       
       if nMatches(b,n) > Settings.Events.SameColorProximityNumberOfBalls
           
           eventsBallColorMatch(b,n) = 1;
           
       end
       
   end
    
end

%% Mixing Events

eventsBallMixing = zeros(nBalls,allBalls,nFrames);
  
targetsDistance(1:nBalls,1:allBalls,1:nFrames) = distanceAllBalls(selectedBalls(1:nBalls),1:allBalls,1:nFrames);

nearThreshold = (targetsDistance > Settings.MOT.RepulsionRadius) & (targetsDistance < 2*Settings.MOT.RepulsionRadius);

sameColor = zeros(nBalls,nFrames);
diffColor = zeros(nBalls,nFrames);

for b=1:nBalls
   
   for bb=1:allBalls
      
       for n=1:(nFrames - Settings.Events.MixingMovingAroundThresholdTime*FrameRate)
           
%            if sweep(selectedBalls(b),bb,n) > Settings.Events.MixingSweepThreshold
%                
%                selectedBall = [xi(selectedBalls(b),n), yi(selectedBalls(b),n)];
%                allBall = [xi(bb,n), yi(bb,n)];
%                X = [selectedBall;allBall];
%                euclidean = pdist(X);
%                
%                if euclidean < Settings.Events.MixingProximityDistanceThreshold
%                
%                     eventsBallMixing(b,bb,n) = 1;
%                     
%                end
%                
%            end
           
            if selectedBalls(b) ~= bb
                
                if all(squeeze(nearThreshold(b,bb,n:(n+Settings.Events.MixingMovingAroundThresholdTime*FrameRate))),1)

                    eventsBallMixing(b,bb,n) = 1;

                end
            
            end

       end
       
   end
   
end

for n=2:(nFrames - Settings.Events.MixingMovingAroundThresholdTime*FrameRate)
    
    for b=1:nBalls
        
        for bb=1:allBalls
            
            if eventsBallMixing(b,bb,n)
                
                if isequal(BallColors(selectedBalls(b),:),BallColors(bb,:))
                    
                    sameColor(b,n) = sameColor(b,n) + 1;
                    
                else
                    
                    diffColor(b,n) = diffColor(b,n) + 1;
                    
                end
                
            elseif eventsBallMixing(b,bb,n-1) && ~eventsBallMixing(b,bb,n)
                
                if isequal(BallColors(selectedBalls(b),:),BallColors(bb,:))

                    sameColor(b,n:(n+Settings.Events.MixingMovingAroundThresholdTime*FrameRate)) = sameColor(b,n:(n+Settings.Events.MixingMovingAroundThresholdTime*FrameRate)) + 1;

                else

                    diffColor(b,n:(n+Settings.Events.MixingMovingAroundThresholdTime*FrameRate)) = diffColor(b,n:(n+Settings.Events.MixingMovingAroundThresholdTime*FrameRate)) + 1;

                end
            
            end
            
        end
        
    end
    
end
   
    