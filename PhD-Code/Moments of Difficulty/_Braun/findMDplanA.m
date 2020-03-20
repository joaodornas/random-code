
function findMDplanA

%%% THRESHOLDS FOR MOMENTS OF DIFFICULTY

VelocitiesThresholds(1).minArea_threshold = 5000;
VelocitiesThresholds(1).maxCollinearity_threshold = 0.1;
VelocitiesThresholds(1).minColorEntropy_threshold = 15;
VelocitiesThresholds(1).minMotionWeightSimilarity_threshold = -9.25;

VelocitiesThresholds(2).minArea_threshold = 5000;
VelocitiesThresholds(2).maxCollinearity_threshold = 0.1;
VelocitiesThresholds(2).minColorEntropy_threshold = 13;
VelocitiesThresholds(2).minMotionWeightSimilarity_threshold = -10.5;

VelocitiesThresholds(3).minArea_threshold = 5000;
VelocitiesThresholds(3).maxCollinearity_threshold = 0.1;
VelocitiesThresholds(3).minColorEntropy_threshold = 13;
VelocitiesThresholds(3).minMotionWeightSimilarity_threshold = -11.25;

VelocitiesThresholds(4).minArea_threshold = 5000;
VelocitiesThresholds(4).maxCollinearity_threshold = 0.1;
VelocitiesThresholds(4).minColorEntropy_threshold = 13;
VelocitiesThresholds(4).minMotionWeightSimilarity_threshold = -12;

VelocitiesThresholds(5).minArea_threshold = 5000;
VelocitiesThresholds(5).maxCollinearity_threshold = 0.1;
VelocitiesThresholds(5).minColorEntropy_threshold = 13;
VelocitiesThresholds(5).minMotionWeightSimilarity_threshold = -12.8;

VelocitiesThresholds(6).minArea_threshold = 5000;
VelocitiesThresholds(6).maxCollinearity_threshold = 0.1;
VelocitiesThresholds(6).minColorEntropy_threshold = 13;
VelocitiesThresholds(6).minMotionWeightSimilarity_threshold = -13.5;

VelocitiesThresholds(7).minArea_threshold = 5000;
VelocitiesThresholds(7).maxCollinearity_threshold = 0.1;
VelocitiesThresholds(7).minColorEntropy_threshold = 13;
VelocitiesThresholds(7).minMotionWeightSimilarity_threshold = -14.25;

%%% MOT Parameters

velocities = [20 25 30 35 40 45 50];
nMOT = 10;
nK = 13200;
nMD = 3;
nComb = 70;
nRGB = 2;

TargetN = 4;

prefix = 'MOT-Analyze-V2-11-m-';

%%% Load Balls and Colors


%% COMPUTE MOMENTS OF DIFFICULTY

moments_of_difficulty = zeros(nK,nMD,nComb,nRGB,nMOT,length(velocities));
balls_motion_of_difficulty = zeros(nK,TargetN,nComb,nRGB,nMOT,length(velocities));

idxVelocity = 0;

for iVelocity=velocities
   
    disp(strcat('Velocity:',int2str(iVelocity)));
    
    idxVelocity = idxVelocity + 1;
    
    for iMOT=1:nMOT
        
        disp(strcat('MOT:',int2str(iMOT)));
  
        analyzedata = load(strcat(prefix,int2str(iVelocity),'-',int2str(iMOT),'.mat'));
        analyzeareadata = load(strcat(prefix,'Area','-',int2str(iVelocity),'-',int2str(iMOT),'.mat'));
        analyzecollineardata = load(strcat(prefix,'Collinear','-',int2str(iVelocity),'-',int2str(iMOT),'.mat'));
                
        allColorComb = analyzeareadata.allColorComb;
        
        for iColor=1:nRGB
        
            for iComb=1:size(allColorComb(1).comb,1)
                
                balls = allColorComb(iColor).comb(iComb,:);
                
                Ecolor = analyzedata.Ecolor;
                Wmotion = analyzedata.Wmotion;
                
                areaComb = squeeze(analyzeareadata.areaComb(iColor,iComb,:));
                collinearComb = squeeze(analyzecollineardata.collinearComb(iColor,iComb,:));
                
                color_difficulty = Ecolor < VelocitiesThresholds(idxVelocity).minColorEntropy_threshold;
                
                motion_difficulty = Wmotion < VelocitiesThresholds(idxVelocity).minMotionWeightSimilarity_threshold;
                
                moments_of_difficulty(:,1,iComb,iColor,iMOT,idxVelocity) = color_difficulty;
                
                moments_of_difficulty(:,2,iComb,iColor,iMOT,idxVelocity) = motion_difficulty(balls(1),:) + motion_difficulty(balls(2),:) + motion_difficulty(balls(3),:) + motion_difficulty(balls(4),:);
                
                balls_motion_of_difficulty(:,1,iComb,iColor,iMOT,idxVelocity) = motion_difficulty(balls(1),:) * balls(1);
                balls_motion_of_difficulty(:,2,iComb,iColor,iMOT,idxVelocity) = motion_difficulty(balls(2),:) * balls(2);
                balls_motion_of_difficulty(:,3,iComb,iColor,iMOT,idxVelocity) = motion_difficulty(balls(3),:) * balls(3);
                balls_motion_of_difficulty(:,4,iComb,iColor,iMOT,idxVelocity) = motion_difficulty(balls(4),:) * balls(4);
                
                area_difficulty = areaComb < VelocitiesThresholds(idxVelocity).minArea_threshold;
                collinear_difficulty = collinearComb < VelocitiesThresholds(idxVelocity).maxCollinearity_threshold;
                
                cross_difficulty = area_difficulty == 1 & collinear_difficulty == 1;
                
                moments_of_difficulty(:,3,iComb,iColor,iMOT,idxVelocity) = cross_difficulty;
                
            end
            
        end
        
    end
    
end

save('MDplanA.mat','moments_of_difficulty','balls_motion_of_difficulty','VelocitiesThresholds','-v7.3');

return

function allColorComb = getColorsComb(datafilename)

    load(datafilename);
    
    TargetN = 4;
    
    BallColors = MOT.BallColors;
    
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

return



