
%clear all

data = load('MOT.mat');

xi = data.MOT.xi;
yi = data.MOT.yi;
vi = data.MOT.vi;
oi = data.MOT.oi;
vix = data.MOT.vixx;
viy = data.MOT.viyy;
BallColors = data.MOT.BallColors;

Settings.MOT.Duration= 60; %% in seconds

FrameRate = floor(size(xi,2)/Settings.MOT.Duration);

nBalls = 3;
allBalls = 16;

Settings.MOT.Target.AllRGB=[ 200 30 30; 250 150 30; 30 150 30; 50 50 255 ];
Settings.MOT.RepulsionRadius = 50;   % radius of repulsive interactions % pix
DirectionDifference = 10;
VelocityPercentage = 0.9;
AreaThreshold = 0.1;
idx_selectedColor = 1;

ib = 1;
for b=1:allBalls
   
    if isequal(BallColors(b,:),Settings.MOT.Target.AllRGB(idx_selectedColor,:))
        
        selectedBalls(ib) = b;
    
        ib = ib + 1;
        
    end
    
end

selectedBalls = selectedBalls(1:nBalls);

nFrames = size(data.MOT.xi,2);

viBalls(1:nBalls,1:nFrames) = vi([selectedBalls(1) selectedBalls(2) selectedBalls(3)],1:nFrames);

for n=1:nFrames
   
    areaBalls(n) = tri_area([xi(selectedBalls(1),n) yi(selectedBalls(1),n)],[xi(selectedBalls(2),n) yi(selectedBalls(2),n)],[xi(selectedBalls(3),n) yi(selectedBalls(3),n)]);
    
end

for n=1:nFrames
    
    for b=1:nBalls
        
        distance(:,b,n) = hypot(xi(selectedBalls(b),n) - xi(:,n),yi(selectedBalls(b),n) - yi(:,n));
        
    end
    
end

for b=1:nBalls
    
   for n=1:nFrames
      
        near = find(distance(:,b,n) < Settings.MOT.RepulsionRadius);    
       
        ballColorMatch(b,n) = 0;
        ballDirectionMatch(b,n) = 0;
        
        near(near == selectedBalls(b)) = [];
        
        vball = [vix(selectedBalls(b),n) viy(selectedBalls(b),n)];
        
        for nn=1:length(near)
           
            if isequal(BallColors(near(nn),:),BallColors(selectedBalls(b),:))
                
                ballColorMatch(b,n) = 1;
            
            end
            
            vnear = [vix(near(nn),n) viy(near(nn),n)];
            
            teta = radtodeg( acos ( dot( vball, vnear ) / ( norm(vball)*norm(vnear) ) ) );
            
            if teta < DirectionDifference
                
                ballDirectionMatch(b,n) = 1;
                
            end
            
        end
        
   end
    
end
    
similarityBalls(1:nBalls,:) = (ballDirectionMatch(1:nBalls,:) == 1) & (ballColorMatch(1:nBalls,:) == 1);

maxvi = max(max(viBalls(:,:)));
maxArea = max(areaBalls);

eventsBalls = zeros(nBalls,Settings.MOT.Duration);
eventsAllBalls = zeros(1,Settings.MOT.Duration);

for f=1:Settings.MOT.Duration
    
    for b=1:nBalls
    
        if sum(viBalls(b,(1 + FrameRate*(f-1)):FrameRate*f) > VelocityPercentage*maxvi) ~= 0
            
            eventsBalls(b,f) = 1;
    
        end
        
        if sum(similarityBalls(b,(1 + FrameRate*(f-1)):FrameRate*f)) ~= 0
            
            eventsBalls(b,f) = 1;
            
        end
        
    end
    
    if sum(areaBalls(1 + FrameRate*(f-1):FrameRate*f) < AreaThreshold*maxArea) ~= 0
        
        eventsAllBalls(f) = 1;
        
    end
  
end



