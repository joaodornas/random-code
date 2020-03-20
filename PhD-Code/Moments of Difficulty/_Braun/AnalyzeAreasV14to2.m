function AnalyzeAreasV14to2(saveAnalyze,datafilename,saveAnalyzeDataFileName)

%%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%

disp('Load data file');

load(datafilename);

%%% SET PARAMETERS %%%%%%%%%%%%%%%%%%%%%

% select part of sequence to be analyzed

start  = 61;               % skip first second
finish = length(MOT.xi);   % 

xi = MOT.xi(:,start:finish);
yi = MOT.yi(:,start:finish);
BallColors = MOT.BallColors;
%TargetN = MOT.Settings.MOT.TargetN;
TargetN = 4;
%AllRGB = MOT.Settings.MOT.Target.AllRGB;

% sampling of frames and evaluation window

nFrames = size(xi,2);
sFrames = 3;   % sampling of frames
nBalls = size(xi,1);

nWindow = 18/sFrames;  % number of sampled frames per window

hWindow = ceil( nWindow / 2);   % nWindow is even, hWindow is half ...
nWindow = 2*hWindow;

%%% COMPUTE PARAMETERS  %%%%%%%%%%%%%%%%%%%%%

disp('Compute parameters');

clear AllRGB
AllRGB(1,:) = MOT.Settings.MOT.Target.AllRGB(1,:);
AllRGB(2,:) = MOT.Settings.MOT.Target.AllRGB(3,:);

firstcolor = MOT.Settings.MOT.Target.AllRGB(1,:);
thirdcolor = MOT.Settings.MOT.Target.AllRGB(3,:);

secondcolor = MOT.Settings.MOT.Target.AllRGB(2,:);
fourthcolor = MOT.Settings.MOT.Target.AllRGB(4,:);

for i=1:size(BallColors,1)

    ballcolor = BallColors(i,:);

    if isequal(secondcolor,ballcolor)

        BallColors(i,:) = firstcolor;

    end

    if isequal(fourthcolor,ballcolor)

        BallColors(i,:) = thirdcolor;

    end

end

clear ballsPerColor
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

for i=1:size(AllRGB,1)
    
    selectedBalls = ballsPerColor(i,:);
    
    allColorComb(i).comb = combnk(selectedBalls,TargetN);

end

areaComb = zeros(size(AllRGB,1),size(allColorComb(1).comb,1),ceil(nFrames/sFrames));
area = NaN;

for kf=1:ceil(nFrames/sFrames)
    
    iFrame = (kf-1)*sFrames + 1;
    
    for i=1:size(AllRGB,1)

        for kcomb=1:size(allColorComb(i).comb,1)
            
            balls = allColorComb(i).comb(kcomb,:);
            
            X = zeros(1,length(balls));
            Y = zeros(1,length(balls));
            
            for kb=1:length(balls)
                
                X(kb) = xi(balls(kb),iFrame);
                Y(kb) = yi(balls(kb),iFrame);
                
            end
           
            cv = convhull(X,Y);
            area = polyarea(X(cv),Y(cv));
            
            areaComb(i,kcomb,kf) = area;
            
            if (kf == 1) || (area < minArea.area && ~isnan(area))
                
                minArea.area = area; 
                minArea.color = i;
                minArea.comb = kcomb;
                minArea.kf = kf;
                
            end
            
        end
        
    end
    
end

clear start
clear finish 
clear MOT
clear xi
clear yi
clear k
clear color
clear ballcolor
clear selectedBalls
clear area
clear iFrame
clear X
clear Y
clear balls
clear area
clear i
clear j
clear kb

if saveAnalyze
    
    save(saveAnalyzeDataFileName);
    
end

    
end