
function manyMOTwithAnalysis

firstMOT = input('first MOT:');
lastMOT = input('last MOT:');

param = loadParam;

for n=firstMOT:lastMOT
    
    param.n = n;
    
    disp(strcat('nMOT:',int2str(param.n)));
    
    param.notAccomplished = true;
    
    while param.notAccomplished
    
        disp('Trying');
        
        param.datafilename = strcat('MOT-TrajectoriesAnalysis-',int2str(param.Time),'-',param.Timestr,'-',int2str(param.n));

        if param.newTrajectories; AnalyzeTrajectories(param.loadData,param.plotMovie,param.saveAnalyze,param.datafilename); end
        
        close all
        
        load(param.datafilename);
    
        %% getSegments
        
        nFrames = size(MOT.xi,2);
        
                    %% COLOR

                    minLcolor = Lcolor(kVCmin);
                    maxHcolor = Hcolor(kVCmax);

                    if (1 - minLcolor) > maxHcolor

                        mdC = kVCmin;
                        mdCball = kVCballmin;

                    else

                        mdC = kVCmax;
                        mdCball = kVCballmax;

                    end

                    kemin = mdC;
                    kerange = kemin:kemin+nWindow-1;  % index to sampled frames;
                    jerange = (kerange-1)*sFrames+1;  % index to original frames;

                    jmin = jerange(1);
                    jmax = jerange(end);

                    windowRange = jmin:jmax;

                    start_ = (jmin - 1 - (param.segTime/2)*param.FrameRate):(jmin - 1);
                    end_ = (jmax + 1):(jmax + 1 + (param.segTime/2)*param.FrameRate);

                    mdCrange = [start_ windowRange end_];

                    if start_(1) >= 1 && end_(end) <= nFrames

                        allunderC = true;

                    else

                        allunderC = false;

                    end


                    %% MOVEMENT

                    minLmov = Lmotion(kVMmin);
                    maxLmov = Hmotion(kVMmax);

                    if minLmov > maxLmov

                        mdM = kVMmin;
                        mdMball = kVMballmin;

                    else

                        mdM = kVMmax;
                        mdMball = kVMballmax;

                    end

                    kemin = mdM;
                    kerange = kemin:kemin+nWindow-1;  % index to sampled frames;
                    jerange = (kerange-1)*sFrames+1;  % index to original frames;

                    jmin = jerange(1);
                    jmax = jerange(end);

                    windowRange = jmin:jmax;

                    start_ = (jmin - 1 - (param.segTime/2)*param.FrameRate):(jmin - 1);
                    end_ = (jmax + 1):(jmax + 1 + (param.segTime/2)*param.FrameRate);

                    mdMrange = [start_ windowRange end_];

                    if start_(1) >= 1 && end_(end) <= nFrames

                        allunderM = true;

                    else

                        allunderM = false;

                    end

                    %% to see if they cross each other

                    allrange = [mdCrange mdMrange];

                    count = hist(allrange,min(allrange):max(allrange));

                    if sum(count > 1) == 0

                        rangedoesnotcross = true;

                    else

                        rangedoesnotcross = false;

                    end

                    if allunderC && allunderM && rangedoesnotcross

                        param.notAccomplished = false;

                        ColorMOT = MOT;

                        ColorMOT.xi = ColorMOT.xi(:,mdCrange);
                        ColorMOT.yi = ColorMOT.yi(:,mdCrange);
                        ColorMOT.vi = ColorMOT.vi(:,mdCrange);
                        ColorMOT.vix = ColorMOT.vix(:,mdCrange);
                        ColorMOT.viy = ColorMOT.viy(:,mdCrange);
                        ColorMOT.ai = ColorMOT.ai(:,mdCrange);
                        ColorMOT.oi = ColorMOT.oi(:,mdCrange);
                        
                        ColorMOT.kemin = mdC;
                        ColorMOT.criteria = 'color';
                        ColorMOT.ball = mdCball;
                        
                        MovMOT = MOT;

                        MovMOT.xi = MovMOT.xi(:,mdMrange);
                        MovMOT.yi = MovMOT.yi(:,mdMrange);
                        MovMOT.vi = MovMOT.vi(:,mdMrange);
                        MovMOT.vix = MovMOT.vix(:,mdMrange);
                        MovMOT.viy = MovMOT.viy(:,mdMrange);
                        MovMOT.ai = MovMOT.ai(:,mdMrange);
                        MovMOT.oi = MovMOT.oi(:,mdMrange);

                        MovMOT.kemin = mdM;
                        MovMOT.criteria = 'movement';
                        MovMOT.ball = mdMball;
                        
                        clear MOT
                        
                        MOT = ColorMOT;
                        save(strcat('MOT-Color-',int2str(param.segTime),'-',param.segTimestr,'-',int2str(param.n),'.mat'),'MOT');
                        clear MOT
                        
                        MOT = MovMOT;
                        save(strcat('MOT-Mov-',int2str(param.segTime),'-',param.segTimestr,'-',int2str(param.n),'.mat'),'MOT');
                        clear MOT
                        
                    else

                        if ~allunderC; disp(strcat('MOT:',int2str(param.n),':color range negative')); end

                        if ~allunderM; disp(strcat('MOT:',int2str(param.n),':movement range negative')); end

                        if ~rangedoesnotcross; disp(strcat('MOT:',int2str(param.n),':range does cross')); end

                    end
                    
                    if param.newTrajectories; save('param.mat','param'); end

        clear all
        close all
        
        newparam = load('param.mat');
        
        param = newparam.param;
        
        if param.notAccomplished && param.newTrajectories; delete(strcat(param.datafilename,'.mat')); end
    
    end
    
    clear all
    close all
        
    param = loadParam;
    
end

return

function param = loadParam

    param.nMOT = 100;
    param.Time = 5;
    param.Timestr = 'm';

    param.loadData = false;
    param.plotMovie = false;
    
    param.saveAnalyze = true;
   
    param.segTime = 10;
    param.segTimestr = 's';
    param.FrameRate = 60;
    
    param.notAccomplished = true;
    
    param.newTrajectories = true;

return