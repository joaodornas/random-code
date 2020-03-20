
%% EYE TRACKER
%% starting recording  
if (Settings.Eyelink.Use)
  Eyelink('StartRecording');
  WaitSecs(0.1); % record a few samples before we actually start displaying
  Eyelink('Message', 'KEY_EVENT TrialOnset');
end;

%% DEFINE CONTROL VARIABLES
hasNotStarted = true;
hasNotShowedStartedMessage = true;
StillHaveThingsToRun = true;
hasPresentedCueTargets = false;
hasWaitingTimeCueTargetsPassed = false;
startMovingBalls = false;
presentTheProbe = false;
hasPresentedTheProbe = false;
showFeedBackNow = false;
hasShowedFeedBack = false;

ItemX = MOT.xi;
ItemY = MOT.yi;
BallColors = MOT.BallColors;
%MOT.FramesN = size(MOT.xi,2) - Settings.Screen.ScreenMode.hz*1;
MOT.FramesN = size(MOT.xi,2);

targets = randperm(Settings.MOT.Target.NumOfBallsPerColor);
targets = targets(1:Settings.MOT.TargetN);

TargetColor = Settings.MOT.Target.AllRGB(Settings.MOT.Target.SelectedColor,:);

if all_trials_w_conditions{iTrial,3} == 100
    
    frameInterval = 1:MOT.FramesN;

else
   
    startFirstHalf = 1;
    endFirstHalf = (MOT.FramesN - Settings.MOT.Segments.FramesInWindow)/2;
    
    startLastHalf = (MOT.FramesN - Settings.MOT.Segments.FramesInWindow)/2 + Settings.MOT.Segments.FramesInWindow + 1;
    endLastHalf = MOT.FramesN;
    
    FrameRate = Settings.Screen.ScreenMode.hz;
    
    percentage = all_trials_w_conditions{iTrial,3}/100;
    finalTimeDuration = percentage * Settings.MOT.Segments.Duration;
    totalFrames = finalTimeDuration*FrameRate;
    
    startFrameInterval = linspace(1,endFirstHalf,totalFrames/2);
    endFrameInterval = linspace(startLastHalf,MOT.FramesN,totalFrames/2);
    
    frameInterval = [startFrameInterval endFrameInterval(2:end)];
    
    frameInterval = round(frameInterval);
    
end

iiTarget = 1;
iiotherTarget = 1;
i = 1;
while i <= Settings.MOT.TotalN
   
    nextColor = MOT.BallColors(i,:);
    
    if isequal(TargetColor,nextColor)
        
        if iiTarget <= Settings.MOT.TargetN
        
            Targets(iiTarget) = i;
        
            iiTarget = iiTarget + 1;
            
        else
            
            otherTargets(iiotherTarget) = i;
            
            iiotherTarget = iiotherTarget + 1;
            
        end
        
    end
    
    i = i + 1;
    
end

Log.Block{end}.Trial(iTrial).Targets = Targets;

while (StillHaveThingsToRun)

    if hasNotStarted
        
        if hasNotShowedStartedMessage
            
            %% message before the onset

            ptbScreen.ShowCenteredMessage({Settings.Message.Trial.StartTrialInstrunctionMsg}, [255 255 255]);
            [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, Beampos] = Screen('Flip', ptbScreen.WinID);
            FlipLogLabel = Settings.Message.FlipLogLabel.StartTrialInstrunctionMsg;
            LogFlip;

            if (Settings.Eyelink.Use)
                Eyelink('Message', 'KEY_EVENT TrialMessage');
            end;
            
            hasNotShowedStartedMessage = false;
            
        end
        
 
        GotTTL= false;
        %% waiting for TR
        while (~GotTTL)
          
             ReadKeys;
                
             if (tSmallKeyTime) || (tCapitalKeyTime)
                
                GotTTL = true;
                    
             end
          
        end
        
        hasNotStarted = false;
       
        TTLLogLabel = Settings.Message.TTLLogLabel.StartTrial;
        LogTTL;
        
    elseif ~hasWaitingTimeCueTargetsPassed
        
        if ~hasPresentedCueTargets
            
            timePresentedCueTargets = GetSecs;
            
            %% presenting cued targets
            Width = Settings.MOT.Cue.Width;
            Color = Settings.MOT.Cue.RGB;
            iFrame = 1;
            for i=1:Settings.MOT.TargetN
                thisBall = Targets(i);
                DrawCues;
            end
            DrawBalls;
            Screen('FrameOval', ptbScreen.WinID, Settings.Fixation.Color, [Settings.Screen.CenterXY(1)-Settings.Fixation.Size/2 Settings.Screen.CenterXY(2)-Settings.Fixation.Size/2 Settings.Screen.CenterXY(1)+Settings.Fixation.Size/2 Settings.Screen.CenterXY(2)+Settings.Fixation.Size/2]);
            [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, Beampos] = Screen('Flip', ptbScreen.WinID);
            FlipLogLabel = Settings.Message.FlipLogLabel.CueTargets;
            LogFlip;
            
            if (Settings.Eyelink.Use)
              Eyelink('Message', 'KEY_EVENT CueTargets');
            end;

            hasPresentedCueTargets = true;
            
        end
        
        timeNow = GetSecs;
        
        if (timeNow - timePresentedCueTargets) > Settings.Schedule.TargetCue
            
            hasWaitingTimeCueTargetsPassed = true;
            
            startMovingBalls = true;
            
        end
        
        ReadKeys;
       
        if (tSmallKeyTime) || (tCapitalKeyTime)
              
            TTLLogLabel = Settings.Message.TTLLogLabel.CueTargets;
            LogTTL;
            
        end
        
    elseif startMovingBalls

        %% moving targets around
        if (Settings.Eyelink.Use) 
          Eyelink('Message', 'KEY_EVENT StartTracking');
        end;

        TrackingStart= GetSecs;
        
        %for iFrame= 1:MOT.FramesN,
        
        for iFrame = frameInterval
            
            DrawBalls;
            Screen('FrameOval', ptbScreen.WinID, Settings.Fixation.Color, [Settings.Screen.CenterXY(1)-Settings.Fixation.Size/2 Settings.Screen.CenterXY(2)-Settings.Fixation.Size/2 Settings.Screen.CenterXY(1)+Settings.Fixation.Size/2 Settings.Screen.CenterXY(2)+Settings.Fixation.Size/2]);
            [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, Beampos] = Screen('Flip', ptbScreen.WinID);
            FlipLogLabel = strcat(Settings.Message.FlipLogLabel.MOTFrame,int2str(iFrame));
            LogFlip;

            ReadKeys;
       
            if (tSmallKeyTime) || (tCapitalKeyTime)
              
                TTLLogLabel = strcat(Settings.Message.FlipLogLabel.MOTFrame,int2str(iFrame));
                LogTTL;
                
            end
            
        end
        
        startMovingBalls = false;
        presentTheProbe = true;
        
    elseif presentTheProbe

        if ~hasPresentedTheProbe
            
            %% presenting the probe
            Width = Settings.MOT.Cue.Width;
            Color = Settings.MOT.Probe.RGB;
            iFrame = MOT.FramesN;
            newTargets = Targets;
            if ~Trial.ProbeTarget(iTrial)
                newTarget = otherTargets(randi(length(otherTargets)));
                newTargets = Targets;
                newTargets(randi(length(Targets))) = newTarget;
            end
            for i=1:Settings.MOT.TargetN
                thisBall = newTargets(i);
                DrawCues;
            end
            DrawBalls;
            Screen('FrameOval', ptbScreen.WinID, Settings.Fixation.Color, [Settings.Screen.CenterXY(1)-Settings.Fixation.Size/2 Settings.Screen.CenterXY(2)-Settings.Fixation.Size/2 Settings.Screen.CenterXY(1)+Settings.Fixation.Size/2 Settings.Screen.CenterXY(2)+Settings.Fixation.Size/2]);
            [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, Beampos] = Screen('Flip', ptbScreen.WinID);
            FlipLogLabel = Settings.Message.FlipLogLabel.PresentProbe;
            LogFlip;

            if (Settings.Eyelink.Use)
              Eyelink('Message', 'KEY_EVENT Probe');
            end;
            
            hasPresentedTheProbe = true;
            
        end

        %% getting the response
        GotTheAnswer= false;
        ResponseStart= GetSecs;
        while (GetSecs-ResponseStart< Settings.Schedule.ProbeCue)
            
          if (~GotTheAnswer)
              
     ReadKeys;
    
            if (keyIsDown)
 
              if (RighKeyTime) || (ResponseStick_YellowKeyTime ~= 0)
                  
                Log.Block{end}.Response(iTrial)= 1;
                Log.Block{end}.Correct(iTrial)= Trial.ProbeTarget(iTrial)==Log.Block{end}.Response(iTrial);
                Log.Block{end}.RT(iTrial)= GetSecs-ResponseStart;
                
              end
                
              if (LeftKeyTime) || (ResponseStick_GreenKeyTime)
                  
                Log.Block{end}.Response(iTrial)= 0;
                Log.Block{end}.Correct(iTrial)= Trial.ProbeTarget(iTrial)==Log.Block{end}.Response(iTrial);
                Log.Block{end}.RT(iTrial)= GetSecs-ResponseStart;
                
             end;
             
             if (tSmallKeyTime) || (tCapitalKeyTime)
                 
                 TTLLogLabel = Settings.Message.TTLLogLabel.PresentProbe;
                 LogTTL;

             end
             
           end;
           
          end
          
        end
        
        presentTheProbe = false;
        showFeedBackNow = true;
        
    elseif showFeedBackNow
        
        if ~hasShowedFeedBack
            
            %% feedback
            timeShowedFeedBack = GetSecs;
            Width = Settings.MOT.Cue.Width;
            Color = Settings.MOT.Cue.RGB;
            iFrame = MOT.FramesN;
            for i=1:Settings.MOT.TargetN
                thisBall = Targets(i);
                DrawCues;
            end
            DrawBalls;
            Screen('FrameOval', ptbScreen.WinID, Settings.Fixation.Color, [Settings.Screen.CenterXY(1)-Settings.Fixation.Size/2 Settings.Screen.CenterXY(2)-Settings.Fixation.Size/2 Settings.Screen.CenterXY(1)+Settings.Fixation.Size/2 Settings.Screen.CenterXY(2)+Settings.Fixation.Size/2]);
            [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, Beampos] = Screen('Flip', ptbScreen.WinID);
            FlipLogLabel = Settings.Message.FlipLogLabel.PresentFeedback;
            LogFlip;

            if (Settings.Eyelink.Use)
              Eyelink('Message', 'KEY_EVENT Feedback');
            end;
            
            hasShowedFeedBack = true;
            
        end
        
        %WaitSecs(Settings.Schedule.Feedback);
        
        timeNow = GetSecs;
        
        if (timeNow - timeShowedFeedBack) > Settings.Schedule.Feedback
            
            showFeedBackNow = false;
            StillHaveThingsToRun = false;
            
        end
        
        ReadKeys;

        if (tSmallKeyTime) || (tCapitalKeyTime)
                
            TTLLogLabel = Settings.Message.TTLLogLabel.PresentFeedback;
        	LogTTL;

        end
        
    else
            
        StillHaveThingsToRun = false;
            
    end
        
end;


