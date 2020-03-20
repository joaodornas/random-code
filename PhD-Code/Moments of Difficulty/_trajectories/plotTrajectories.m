
MOTDataFile = input('Please type the MOT datafile name: ','s');

load(MOTDataFile);

xi = MOT.xi + MOT.Settings.Screen.CenterXY(1);
yi = MOT.yi + MOT.Settings.Screen.CenterXY(2);

nFrames = size(xi,2);
nBalls = size(xi,1);

startStr = input('First Frame (i.e. 1): ','s');
finishStr = input(strcat('Last Frame (max. ',int2str(nFrames),'): '),'s');

start = str2num(startStr);
finish = str2num(finishStr);

color{1} = 'r';
color{2} = 'b';
color{3} = 'g';
color{4} = 'k';
color{5} = 'y';

balls = 1:nBalls;

BallColors = MOT.BallColors;
AllRGB = MOT.Settings.MOT.Target.AllRGB;
nColors = size(AllRGB,1);

for b=1:nBalls

    thisBallColor = BallColors(b,:);
    idx_thisBallColor = find(ismember(AllRGB,thisBallColor,'rows'));
    colorBalls(b) = idx_thisBallColor;
    
end
 
figure

for b=1:length(balls)
    
   plot(xi(balls(b),start:finish),yi(balls(b),start:finish),color{colorBalls(b)});
   
   hold on
    
end