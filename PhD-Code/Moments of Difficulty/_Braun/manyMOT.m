
firstMOT = input('first MOT:');
lastMOT = input('last MOT:');

for n=firstMOT:lastMOT
    
    disp(strcat('nMOT:',int2str(n)));
    
    tic
    
    MOT = GenerateBraunianWalkTrajectories;
    
    Time = MOT.Settings.MOT.DurationInMinutes;
    
    save(strcat('MOT-',int2str(Time),'-','m','-',int2str(n),'.mat'),'MOT');
    
    toc
    
    clear MOT
    clear Time
    
end