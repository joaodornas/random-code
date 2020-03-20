%% getting current block settings
Settings.SettingsForBlock(iBlock);
 
%% initializing logger
Log.NewBlock({'MeasuredFPS', 'ProbeIsTarget', 'Response', 'RT', 'Correct', 'Trials_Selected', 'nCondition', 'all_trials_w_conditions'});
Log.Block{end}.MeasuredFPS= ptbScreen.Mode.hz;

CenterXY = Settings.Screen.CenterXY;

%% randomizing trials

if Settings.MOT.Segments.SpeedUp
    
    speeds = Settings.MOT.Segments.Speeds;
    nSpeeds = length(speeds);
    
else
    
    speeds = 100;
    nSpeeds = length(speeds);
    
end

nCondition = Settings.MOT.Segments.MomOfDiff.Color.Use + Settings.MOT.Segments.MomOfDiff.Movement.Use;

trials_selected = randi([1 Settings.MOT.Segments.NumberOfSegments],1,Settings.Trials*nCondition*nSpeeds);

all_trials_w_conditions = cell(Settings.Trials*nCondition*nSpeeds,3);

all_trials_w_conditions(:,1) = num2cell(trials_selected);

if Settings.MOT.Segments.MomOfDiff.Color.Use && ~Settings.MOT.Segments.MomOfDiff.Movement.Use

    for s=1:nSpeeds
        
        all_trials_w_conditions{((s-1)*Settings.Trials*nCondition+1):s*Settings.Trials*nCondition,2} = Settings.MOT.Segments.MomOfDiff.Color.String;
        all_trials_w_conditions{((s-1)*Settings.Trials*nCondition+1):s*Settings.Trials*nCondition,3} = speeds(s);
        
    end
    
elseif Settings.MOT.Segments.MomOfDiff.Movement.Use && ~Settings.MOT.Segments.MomOfDiff.Color.Use

    for s=1:nSpeeds
        
        all_trials_w_conditions{((s-1)*Settings.Trials*nCondition+1):s*Settings.Trials*nCondition,2} = Settings.MOT.Segments.MomOfDiff.Movement.String;
        all_trials_w_conditions{((s-1)*Settings.Trials*nCondition+1):s*Settings.Trials*nCondition,3} = speeds(s);
        
    end
    
elseif Settings.MOT.Segments.MomOfDiff.Movement.Use && Settings.MOT.Segments.MomOfDiff.Color.Use

    for s=1:nSpeeds
        
        stepSpeed = (s-1)*Settings.Trials*nCondition;
        
        all_trials_w_conditions{1+stepSpeed:Settings.Trials+stepSpeed,2} = Settings.MOT.Segments.MomOfDiff.Color.String;
        all_trials_w_conditions{1+stepSpeed:Settings.Trials+stepSpeed,3} = speeds(s);
 
        all_trials_w_conditions{Settings.Trials+1+stepSpeed:Settings.Trials*nCondition+stepSpeed,2} = Settings.MOT.Segments.MomOfDiff.Movement.String;
        all_trials_w_conditions{Settings.Trials+1+stepSpeed:Settings.Trials*nCondition+stepSpeed,3} = speeds(s);
    
    end
    
end

all_trials_w_conditions = all_trials_w_conditions(randperm(Settings.Trials*nCondition*nSpeeds),:);

Trial.ProbeTarget = randi([0 1],1,Settings.Trials*nCondition*nSpeeds);

%% INITIALIZE FLIP AND TRIGGER TIMESTAMP
%% book-keeping
Log.Block{end}.ttlTime = cell(50000,5);
iTTL= 1;
Log.Block{end}.vblTime= cell(50000,8);
iVBL= 1;

%% logging what we know, prepare arrays for what we do not know
Log.Block{end}.ProbeIsTarget = Trial.ProbeTarget;
Log.Block{end}.Response = nan(1, Settings.Trials);
Log.Block{end}.RT = nan(1, Settings.Trials);
Log.Block{end}.Correct = nan(1, Settings.Trials);
Log.Block{end}.Trials_Selected = trials_selected;
Log.Block{end}.nCondition = nCondition;
Log.Block{end}.all_trials_w_conditions = all_trials_w_conditions;

%% Calibrate the eye tracker
if (Settings.Eyelink.Use)
  EyelinkDoTrackerSetup(EyelinkWindowHandle);
end;

%% SHOW START BLOCK MSG
notNextBlock = true;
hasShowedWelcomeMessage = false;
itIsNotTimeToGoAhead = true;
waitForTR = false;

while notNextBlock

    %% FIRST BLOCK MSG
    if iBlock == 1

        if ~hasShowedWelcomeMessage
            
            %% finally, welcome to block message
            Strings= {sprintf(Settings.Message.Block.Start, iBlock, Settings.BlockN), Settings.Message.Block.EnterInput};
            ptbScreen.ShowCenteredMessage(Strings, [255 255 255]);
            [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, Beampos] = Screen('Flip', ptbScreen.WinID);
            FlipLogLabel = Settings.Message.FlipLogLabel.StartFirstBlock;
            
            LogFlip;
            
            hasShowedWelcomeMessage = true;
            
        end
        
        while itIsNotTimeToGoAhead
            
            ReadKeys;
            
            if (keyIsDown)
                
              if (EnterKeyTime)
                  
                   waitForTR = true;
                    
              end
              
            end
            
            if waitForTR
                
                %% waiting for TTL ('t' coming from keyboard)
                Strings= {sprintf(Settings.Message.Block.AboutToStart,Settings.fMRI.TR)};
                ptbScreen.ShowCenteredMessage(Strings, [255 255 255]);
                [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, Beampos] = Screen('Flip', ptbScreen.WinID);
                FlipLogLabel = Settings.Message.FlipLogLabel.AboutToStart;
                LogFlip;
                
                GotTTL= false;
                
                while (~GotTTL)
                  ReadKeys;
                  if (tSmallKeyTime) || (tCapitalKeyTime)
                      GotTTL = true;
                  end
                end;

                TTLLogLabel = Settings.Message.TTLLogLabel.BlockAboutToStart;
                LogTTL;
                
                itIsNotTimeToGoAhead = false;
                    
                notNextBlock = false;

            end
            
        end

    %% INTER-BLOCK MSG
    else

        if ~hasShowedWelcomeMessage
            
            if Settings.NoPauseBtwBlocks
                
                Strings = {sprintf(Settings.Message.Block.Start, iBlock, Settings.BlockN)};
                
            else
                
                Strings= {sprintf(Settings.Message.Block.Start, iBlock, Settings.BlockN), sprintf(Settings.Message.Block.PauseBtwBlocks,Settings.PauseBtwBlocksInSeconds)};
            
            end
            
            %% welcome to next block message
            ptbScreen.ShowCenteredMessage(Strings, [255 255 255]);
            [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, Beampos] = Screen('Flip', ptbScreen.WinID);
            FlipLogLabel = Settings.Message.FlipLogLabel.InterBlockMsg;
            LogFlip;
            
            hasShowedWelcomeMessage = false;
            
        end
        
        startTime = GetSecs;
         
        while itIsNotTimeToGoAhead
            
            if Settings.NoPauseBtwBlocks
            
                GotTTL = false;
                
                while ~GotTTL

                    ReadKeys;
                    if (tSmallKeyTime) || (tCapitalKeyTime)
                        GotTTL = true;
                    end
                   
                    if GotTTL
  
                        TTLLogLabel = Settings.Message.TTLLogLabel.BlockAboutToStart;
                        LogTTL;
                        
                        itIsNotTimeToGoAhead = false;
                        notNextBlock = false;
                
                    end
                    
                end
                
            else
                
                timeNow = GetSecs;
            
                if (timeNow - startTime) > Settings.PauseBtwBlocksInSeconds
            
                    itIsNotTimeToGoAhead = false;
                    notNextBlock = false;
                
                end

                GotTTL = false;
                ReadKeys;
                
                if (tSmallKeyTime) || (tCapitalKeyTime)
                
                    GotTTL = true;
                    
                end
                
                if GotTTL
                
                    TTLLogLabel = Settings.Message.TTLLogLabel.PauseBtwBlocks;
                    LogTTL;

                end
                
            end
            
        end

    end

end
