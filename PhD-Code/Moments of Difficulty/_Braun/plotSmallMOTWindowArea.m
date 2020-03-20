
sFrames = 3;   % sampling of frames
nWindow = 18/sFrames;  % number of sampled frames per window

FrameRate = 60;
Time = 10;
fs = 12;

% display size and repulsion radius

sizeX = 400;
sizeY = 400;
repR = 50;

start  = 1;               % skip first second
finish = length(MOT.xi);   % 

xi = MOT.xi(:,start:finish);
yi = MOT.yi(:,start:finish);

AllRGB = MOT.Settings.MOT.Target.AllRGB;
BallColors = MOT.BallColors;

targetBalls = MOT.Targets;

start_ = Time*FrameRate/2;
end_ = Time*FrameRate/2 + nWindow*sFrames;

nBalls = size(xi,1);

figure;
    
set(gca,'nextplot','replacechildren');

colormap jet;  

jerange = start_:end_;  % index to original frames;

color{1} = 'r';
color{2} = 'b';
color{3} = 'k';

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

hold on;

for k=1:length(targetBalls)
    
    plot( xi(targetBalls(k),jerange(1)), yi(targetBalls(k),jerange(1)), 'h', 'Color', color{colorBalls(targetBalls(k))}, 'MarkerSize', 18 );
    colorBalls(targetBalls(k)) = 3;

end

for ib=1:nBalls

    plot( xi(ib,jerange), yi(ib,jerange), 'o', 'Color', color{colorBalls(ib)}, 'MarkerSize', 8 );
    if isempty(find(targetBalls == ib)); plot( xi(ib,jerange(1)), yi(ib,jerange(1)), 'h', 'Color', color{3}, 'MarkerSize', 18 ); end 

end

plot((-sizeX/2):(sizeX/2),0,'k');
plot(0,(-sizeY/2):(sizeY/2),'k');

hold off;

axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

axis 'equal';
axis 'off';

hold on

figure

plot(1:length(MOT.Areas),MOT.Areas);
xlabel('Window idx');
ylabel('Area');
