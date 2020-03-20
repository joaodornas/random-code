
eventsBalls;

if Settings.Events.ShowTrialResponse
    
    Trial = input('Indicate the Trial Number:','s');
    
end

%% initializing PTB screen 
ptbScreen= CPTBScreen(Settings.Screen);
if (Settings.Screen.FullScreen)
    HideCursor;
    ListenChar(2); %% cut-off matlab from key presses
end;

%% creating a USB device input queue
PsychHID('KbQueueCreate');
PsychHID('KbQueueStart');

CenterXY = Settings.Screen.CenterXY;

textDisplacement = 100;

ItemX = xi;
ItemY = yi;

startTime = GetSecs;

side = Settings.MOT.Area.Size(1)/2;

iFrame = 1;
hs = 1;
VBLTimestamp = 0;

eventsDisplacementRunning = false;
turnBallColorMatchOn(1:nBalls) = 0;

selectedFrames = 1:nFrames;

timer = ones(nBalls,allBalls)*1000;

if Settings.Events.ShowTrialResponse
    
    expData = load('experiment.mat');
    
    nBlocks = length(expData.Block);
    
    for i=1:nBlocks
        
        seqTrials = expData.Block{i}.Sequence_Of_Trials;
        
        idx = find(seqTrials == str2num(Trial));
        
        frame(i) = expData.Block{i}.Trial(idx).LostTrackFrame;
        
        before(i) = frame(i) - Settings.Events.ShowTrialResponseBeforeInSeconds*FrameRate;
       
        if before(i) < 0
            
            before(1) = 1;
            
        end
        
    end
    
end

runPreview = Settings.Events.ShowTrialResponseBeforeInSeconds*FrameRate + 1;
runNow = Settings.Events.ShowTrialResponseBeforeInSeconds*FrameRate + 1;
    

for iFrame=selectedFrames
    
    if (iFrame > 1) && (mod(iFrame,FrameRate) == 1); hs = hs + 1; end 
    
    %% SHOW EVENTS
    
    if Settings.Events.ShowEventsAtAll && (iFrame/FrameRate <= Settings.MOT.Duration)

        FontSize = 20;
        Screen('TextSize', ptbScreen.WinID, FontSize);

        if eventsArea(iFrame) && Settings.Events.AreaEvents.Show

               minArea = areaBalls(areaComb,iFrame);

               Color = Settings.Events.Thresholds.Color;
               Screen('DrawText',ptbScreen.WinID,strcat(num2str(round((minArea/maxArea)*100)),'%'),CenterXY(1),CenterXY(2)-side-30,Color,1);
               Color = Settings.Events.AreaEvents.Cue.Color;
               Width = Settings.Events.AreaEvents.Cue.Width;
               Screen('DrawText',ptbScreen.WinID,'Area Event',(CenterXY(1)-side), (CenterXY(2)-side)-30,Color,1);

               for b=1:length(selectedBalls)
                   thisBall = selectedBalls(b);
                   DrawCues;
               end

        end

        for b=1:nBalls

           if Settings.Events.DisplacementEvents.Show 

               if eventsBallDistance(b,iFrame)

                   if ~eventsDisplacementRunning; startPoint(1,b) = xi(selectedBalls(b),iFrame) + CenterXY(1); end
                   endPoint(1,b) = xi(selectedBalls(b),endDisplacement(selectedBalls(b),iFrame)+iFrame) + CenterXY(1); 
                   if ~eventsDisplacementRunning; startPoint(2,b) = yi(selectedBalls(b),iFrame) + CenterXY(2); end
                   endPoint(2,b) = yi(selectedBalls(b),endDisplacement(selectedBalls(b),iFrame)+iFrame) + CenterXY(2); 

                   eventsDisplacementRunning = true;
                   frameDisplaceLimit = endDisplacement(selectedBalls(b),iFrame);
                   positionDisplaceFrame = 1;

               end

               if eventsDisplacementRunning

                   Color = Settings.Events.DisplacementEvents.Cue.Color;
                   Screen('DrawText',ptbScreen.WinID,'Displacement Event',(CenterXY(1)-side), (CenterXY(2)-side) - 60,Color,1);
                   Color = Settings.Events.DisplacementEvents.Trajectory.Color;
                   Line_Width = Settings.Events.DisplacementEvents.Trajectory.Width;
                   Screen('DrawLine',ptbScreen.WinID,Color,startPoint(1,b),startPoint(2,b),endPoint(1,b),endPoint(2,b),Line_Width);
                   Color = Settings.Events.Thresholds.Color;
                   Screen('DrawText',ptbScreen.WinID,strcat(num2str(round((displacement(b,iFrame)/(Settings.MOT.Area.Size(1)))*100)),'%'),CenterXY(1)+b*100,CenterXY(2)-side-60,Color,1);
                   Color = Settings.Events.DisplacementEvents.Cue.Color;
                   Width = Settings.Events.DisplacementEvents.Cue.Width;

                   thisBall = selectedBalls(b);
                   DrawCues;

                   positionDisplaceFrame = positionDisplaceFrame + 1;

                   if positionDisplaceFrame > frameDisplaceLimit

                       eventsDisplacementRunning = false;

                   end

               end

           end

           if Settings.Events.ColorMatchEvents.Show

               if eventsBallColorMatch(b,iFrame)

                   turnBallColorMatchOn(b) = 1;
                   endFrame(b) = iFrame + Settings.Events.SameColorWindowInSeconds*FrameRate;

                   iBall = 1;

                   for bb=1:allBalls

                       if sum(ballColorMatch(selectedBalls(b),bb,iFrame:iFrame+Settings.Events.SameColorWindowInSeconds*FrameRate)) ~= 0

                           balls(iBall) = bb;
                           iBall = iBall + 1;

                       end

                   end

               end

               if turnBallColorMatchOn(b)

                    Color = Settings.Events.ColorMatchEvents.Cue.Color;
                    Width = Settings.Events.ColorMatchEvents.Cue.Width;
                    Screen('DrawText',ptbScreen.WinID,'Color Match Event',(CenterXY(1)-side), (CenterXY(2)-side) - 90,Color,1);

                    for g=1:length(balls)
                        thisBall = balls(g);
                        DrawCues;
                    end

               end

               if iFrame > endFrame(b)

                   turnBallColorMatchOn(b) = 0;

               end

           end 

            if Settings.Events.MixingEvents.Show
                
                for bb=1:allBalls
                    
                    if eventsBallMixing(b,bb,iFrame); timer(b,bb) = 0; end
                    
                    if timer(b,bb) < Settings.Events.MixingMovingAroundThresholdTime*FrameRate
                        
                        startPoint = [(xi(selectedBalls(b),iFrame) + CenterXY(1)) (yi(selectedBalls(b),iFrame) + CenterXY(2))];
                        endPoint = [(xi(bb,iFrame) + CenterXY(1)) (yi(bb,iFrame) + CenterXY(2))];
                        Color = Settings.Events.MixingEvents.Link.Color;
                        Line_Width = Settings.Events.MixingEvents.Link.Width;
                        Screen('DrawLine',ptbScreen.WinID,Color,startPoint(1),startPoint(2),endPoint(1),endPoint(2),Line_Width);
                        
                        timer(b,bb) = timer(b,bb) + 1;
                        
                        thisBall = selectedBalls(b);

                        Color = Settings.Events.MixingEvents.Cue.Color;
                        Width = Settings.Events.MixingEvents.Cue.Width;
                        Screen('DrawText',ptbScreen.WinID,'Mixing Event',(CenterXY(1)-side), (CenterXY(2)-side) - 120,Color,1);
                        DrawCues;

                        Color = Settings.Events.MixingEvents.Cue.Color;
                        Screen('DrawText',ptbScreen.WinID,strcat(int2str(b),':'),(CenterXY(1)-side + 40 + (b-1)*90), (CenterXY(2)-side) - 80,Color,1);
                        Color = [255 255 255];
                        Screen('DrawText',ptbScreen.WinID,int2str(sameColor(b,iFrame)),(CenterXY(1)-side + 40 + 30 + (b-1)*90), (CenterXY(2)-side) - 80,Color,1);
                        Color = [0 0 0];
                        Screen('DrawText',ptbScreen.WinID,int2str(diffColor(b,iFrame)),(CenterXY(1)-side + 40 + 50 + (b-1)*90), (CenterXY(2)-side) - 80,Color,1);
 
                    end
                    
                end
    
            end 

        end
        
    end
    
    if Settings.Events.ShowTrialResponse
        
        if ~isempty(find(before == iFrame,1))

            runPreview = 0;

        end

        if runPreview < Settings.Events.ShowTrialResponseBeforeInSeconds*FrameRate

            runPreview = runPreview + 1;
            Color = [0 0 255];
            Screen('DrawText',ptbScreen.WinID,'Preview',contourRect(3)+10, contourRect(4)+10,Color,1);

        end

        if ~isempty(find(frame == iFrame,1))

            runNow = 0;

        end

        if runNow < Settings.Events.ShowTrialResponseBeforeInSeconds*FrameRate

            runNow = runNow + 1;
            Color = [255 0 0];
            Screen('DrawText',ptbScreen.WinID,'LOST',contourRect(3)+30, contourRect(4)+30,Color,1);

        end
    
    end
       
    if Settings.Events.Distractor.Show

        distractors = selectedBalls(nBalls+1:end);
        Color = Settings.Events.Distractor.Cue.Color;
        Width = Settings.Events.Distractor.Cue.Width;
        DrawDistractors;

    end

    DrawBalls;

    if Settings.Events.ShowRepulsionRadiusForAll
        
        balls = 1:Settings.MOT.TotalN;
        
    else
        
        balls = selectedBalls;
        
    end
    
    DrawRepulsionRadius;

    FontSize = 10;
    Screen('TextSize', ptbScreen.WinID, FontSize);
    Color = Settings.MOT.Target.Label.Color;
    for b=1:nBalls
        Screen('DrawText',ptbScreen.WinID,int2str(b),ItemX(selectedBalls(b), iFrame)+CenterXY(1)+Settings.MOT.Target.Diameter,ItemY(selectedBalls(b), iFrame)+CenterXY(2)+Settings.MOT.Target.Diameter/2,Color,1);
    end

    contourRect = [(CenterXY(1)-side) (CenterXY(2)-side) (CenterXY(1)+side) (CenterXY(2)+side)];
    Screen('FrameRect', ptbScreen.WinID, Settings.MOT.Area.Contour.Color, contourRect, Settings.MOT.Area.Contour.Width);

    %timeElapsed = GetSecs - hsTimeStarted;
    Screen('DrawText',ptbScreen.WinID,strcat(int2str(hs)),contourRect(3), contourRect(4),Settings.Events.Timer.Color,1);

    if Settings.Events.ShowInSlowMotion
        when = VBLTimestamp + Settings.Events.SlowMotionTime;
    else
        %when = VBLTimestamp + 1/FrameRate;
        when = 0;
    end

    [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, Beampos] = Screen('Flip', ptbScreen.WinID, when );

    Flip(iFrame,1) = hs;
    Flip(iFrame,2) = VBLTimestamp;
    Flip(iFrame,3) = StimulusOnsetTime;
    Flip(iFrame,4) = FlipTimestamp;
    Flip(iFrame,5) = Missed;
    Flip(iFrame,6) = Beampos;

    [keyIsDown, firstKeyPressTimes, ~, ~, ~]=PsychHID('KbQueueCheck');

    escapeKeyTime = firstKeyPressTimes(Settings.Keyboard.Escape);
    upKeyTime = firstKeyPressTimes(Settings.Keyboard.Up);

    if keyIsDown || (iFrame/FrameRate > Settings.MOT.Duration)

        if escapeKeyTime ~= 0
        
            break
            
        end
        
        if upKeyTime
            
            iFrame = iFrame + FrameRate;
            hs = hs + 1;
            
        end

    end
    
end

endTime = GetSecs - startTime;

%% closing PTB
ListenChar(0); %% re-enable matlab to see key presses
if (exist('ptbScreen', 'var'))
  ptbScreen.Close;
  ShowCursor;
end;

disp(strcat('Total Time:',num2str(endTime)));
