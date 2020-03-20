
sFrames = 3;   % sampling of frames
nWindow = 18/sFrames;  % number of sampled frames per window

FrameRate = 60;
Time = 10;
fs = 12;

% display size and repulsion radius

sizeX = 400;
sizeY = 400;
repR = 50;

%start  = 61;  % skip first second
start = 1;
finish = length(MOT.xi);   % 

xi = MOT.xi(:,start:finish);
yi = MOT.yi(:,start:finish);

AllRGB = MOT.Settings.MOT.Target.AllRGB;
BallColors = MOT.BallColors;

start_ = Time*FrameRate/2;
end_ = Time*FrameRate/2 + nWindow*sFrames;

nBalls = size(xi,1);

isEntropy = input('Is entropy?');

ball = 1;
if ~isEntropy; ball = MOT.MainTarget; end

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

if ~isEntropy; plot( xi(ball,jerange(1)), yi(ball,jerange(1)), 'h', 'Color', color{colorBalls(ball)}, 'MarkerSize', 18 ); end
if ~isEntropy; colorBalls(ball) = 3; end

for ib=1:nBalls

    plot( xi(ib,jerange), yi(ib,jerange), 'o', 'Color', color{colorBalls(ib)}, 'MarkerSize', 8 );
    if (ib ~= ball) || (isEntropy); plot( xi(ib,jerange(1)), yi(ib,jerange(1)), 'h', 'Color', color{3}, 'MarkerSize', 18 ); end 

end

plot((-sizeX/2):(sizeX/2),0,'k');
plot(0,(-sizeY/2):(sizeY/2),'k');

hold off;

axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

axis 'equal';
axis 'off';

hold on

figure

if isEntropy
    plot(1:length(MOT.Entropies),MOT.Entropies);
    xlabel('Window idx');
    ylabel('Entropy');
else
    plot(1:length(MOT.Weights),MOT.Weights);
    xlabel('Window idx');
    ylabel('Weight');
end