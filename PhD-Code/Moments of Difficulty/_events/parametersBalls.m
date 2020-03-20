    
%% Load Data
Folder.Settings= 'Settings/'; % path to settings folder 

SettingsFile = input('Please enter the Events Settings file name: ', 's');

SettingsEvents= CExperimentalSettings(SettingsFile, Folder.Settings);

MOTDataFile= input('Please enter the MOT datafile name: ', 's');

data = load(MOTDataFile);

Settings.MOT = data.MOT.Settings.MOT;
Settings.Screen = data.MOT.Settings.Screen;
Settings.Events = SettingsEvents.Events;
Settings.Keyboard = SettingsEvents.Keyboard;

xi = data.MOT.xi;
yi = data.MOT.yi;
vi = data.MOT.vi;
oi = data.MOT.oi;
vix = data.MOT.vix;
viy = data.MOT.viy;
BallColors = data.MOT.BallColors;

%% Set Parameters
FrameRate = Settings.Screen.ScreenMode.hz;

allBalls = Settings.MOT.TotalN;

%% Select Balls
ib = 1;
for b=1:allBalls
   
    if isequal(BallColors(b,:),Settings.MOT.Target.AllRGB(Settings.Events.SelectedColor,:))
        
        selectedBalls(ib) = b;
    
        ib = ib + 1;
        
    end
    
end

maxvi = max(max(vi(selectedBalls,:)));
maxoi = max(max(oi(selectedBalls,:)));

nFrames = size(xi,2);

allCombSelectedBalls = combnk(selectedBalls,3);

%% Area Events
for l=1:size(allCombSelectedBalls,1)

    for n=1:nFrames
    
        areaBalls(l,n) = tri_area([xi(allCombSelectedBalls(l,1),n) yi(allCombSelectedBalls(l,1),n)],[xi(allCombSelectedBalls(l,2),n) yi(allCombSelectedBalls(l,2),n)],[xi(allCombSelectedBalls(l,3),n) yi(allCombSelectedBalls(l,3),n)]);
        
    end
    
    maxArea(l) = max(areaBalls(l,:));
    
end

%% Displacement Events

displacement = zeros(allBalls,nFrames,FrameRate*Settings.Events.DisplacementWindowInSeconds);

for s=1:(nFrames-FrameRate*Settings.Events.DisplacementWindowInSeconds)
    
    for b=1:allBalls
        
        for f=1:FrameRate*Settings.Events.DisplacementWindowInSeconds
            
            firstFrame = s;
            lastFrame = s+f;
            displacement(b,s,f) = hypot(xi(b,firstFrame) - xi(b,lastFrame),yi(b,firstFrame) - yi(b,lastFrame));
            
        end
        
        maxDisplacement(b,s) = max(displacement(b,s,1:end));
        endDisplacement(b,s) = find(displacement(b,s,1:end) == maxDisplacement(b,s));
        
    end

end

%% Distance Between All Balls

distanceAllBalls = zeros(allBalls,allBalls,nFrames);
for n=1:nFrames
    
    for b=1:allBalls
        
        distanceAllBalls(:,b,n) = hypot(xi(b,n) - xi(:,n),yi(b,n) - yi(:,n));
        
    end
    
end

%% Color Similarity Events

nearBallColorMatch = zeros(allBalls,allBalls,nFrames);

nearBallColorMatch(:,:,:) = distanceAllBalls(:,:,:) < Settings.Events.SameColorProximityThreshold;
ballColorMatch = nearBallColorMatch;

for b=1:allBalls
    
   for n=1:nFrames
       
        idx = find(nearBallColorMatch(:,b,n));
        
        idx(idx == b) = [];
        
        if ~isempty(idx)
            
            for k=1:length(idx)
        
                if ~isequal(BallColors(idx(k),:),BallColors(b,:))

                    ballColorMatch(idx(k),b,n) = 0;

                end
                
            end
            
        end
        
   end
    
end

%% Mixing Events

nearBallMixing = zeros(allBalls,allBalls,nFrames);

nearBallMixing(:,:,:) = distanceAllBalls(:,:,:) < Settings.Events.MixingProximityDistanceThreshold;

angle = zeros(allBalls,allBalls,nFrames);
sweep = zeros(allBalls,allBalls,nFrames);
endSweep = zeros(allBalls,allBalls,nFrames);

for b=1:allBalls
    
    for bb=1:allBalls
        
        for n=1:nFrames 

            origin = [xi(b,n) yi(b,n)];
                      
            relativeNear = [xi(bb,n) yi(bb,n)] - origin;
                      
            angle(b,bb,n) = radtodeg( atan ( relativeNear(2) / relativeNear(1) ) );
                      
        end
        
        for n=1:(nFrames - Settings.Events.MixingWindowInSeconds*FrameRate)
            
           for s=(n+1):(n + Settings.Events.MixingWindowInSeconds*FrameRate)
               
               this_sweep = abs(sum(diff(angle(b,bb,n:s))));
               
               if this_sweep > sweep(b,bb,n)
                   
                    sweep(b,bb,n) = this_sweep;
                    endSweep(b,bb,n) = s;
               
               end
               
           end
            
        end
    
    end
    
end

clear data;
clear n;
clear idx;
clear b;
clear bb;
clear origin;
clear relativeNear;
clear firstFrame;
clear lastFrame;
clear ib;
clear f;
clear iTrial;
clear k;
clear l;
clear s;
clear this_sweep;

save(strcat('param-',MOTDataFile,'.mat'));
    