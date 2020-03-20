

%% general settings 
Folder.Settings= 'Settings/'; % path to settings folder 

ObserverName = 'prepareWalk-main';

Settings= CExperimentalSettings(ObserverName, Folder.Settings);

%% initializing PTB screen 
ptbScreen= CPTBScreen(Settings.Screen);
 if (Settings.Screen.FullScreen)
    HideCursor;
    ListenChar(2); %% cut-off matlab from key presses
  end;

%% creating a USB device input queue
PsychHID('KbQueueCreate', Settings.fMRI.hidDeviceID);
PsychHID('KbQueueStart', Settings.fMRI.hidDeviceID);

CenterXY = Settings.Screen.CenterXY;

%% LOAD POSITIONS
data = load(strcat('MOT','.mat'));

ItemX = data.MOT.xi;
ItemY = data.MOT.yi;
BallColors = data.MOT.BallColors;

clear data;
active = true;

startTime = GetSecs;

plotBallsInfo;

DrawCues = true;

iFrame = 1;
iSecond = 1;
while active

    DrawBalls;
    
    if DrawCues & (iSecond < Settings.MOT.Duration)
        
        if (iSecond == 1) || (mod(iFrame,FrameRate) == 0)

            shouldWeDrawCues = eventsBalls(:,iSecond);
            shouldWeDrawAllCues = eventsAllBalls(iSecond);

            iSecond = iSecond + 1;

        end

        if shouldWeDrawAllCues; shouldWeDrawCues(:) = shouldWeDrawAllCues; end

        for b=1:nBalls

            if shouldWeDrawCues(b)

                H = [ItemX(selectedBalls(b), iFrame)+CenterXY(1)-Settings.MOT.Cue.Width/2,...
                    ItemY(selectedBalls(b), iFrame)+CenterXY(2)-Settings.MOT.Cue.Width/2,...
                    ItemX(selectedBalls(b), iFrame)+CenterXY(1)+Settings.MOT.Cue.Width/2,...
                    ItemY(selectedBalls(b), iFrame)+CenterXY(2)+Settings.MOT.Cue.Width/2];
                Screen('FillRect', ptbScreen.WinID, Settings.MOT.Cue.RGB, H');

            end

        end
        
    end
    
    Screen('Flip', ptbScreen.WinID);
    
    [keyIsDown, firstKeyPressTimes, ~, ~, ~]=PsychHID('KbQueueCheck',Settings.fMRI.hidDeviceID);
    
    if keyIsDown || (iFrame == size(ItemX,2))
        
        active = false;
        
    end
    
    iFrame = iFrame + 1;
    
end

endTime = GetSecs - startTime;

%% closing PTB
ListenChar(0); %% re-enable matlab to see key presses
if (exist('ptbScreen', 'var'))
  ptbScreen.Close;
  ShowCursor;
end;

disp(strcat('Total Time:',num2str(endTime)));
