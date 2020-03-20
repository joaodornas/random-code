
DataFile= input('Please enter the MOT datafile name: ', 's');

data = load(DataFile);

%% general settings 
Folder.Settings= 'Settings/'; % path to settings folder 

ObserverName = 'showBalls';

Settings= CExperimentalSettings(ObserverName, Folder.Settings);

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

ItemX = data.MOT.xi;
ItemY = data.MOT.yi;
BallColors = data.MOT.BallColors;

Settings.MOT = data.MOT.Settings.MOT;

startTime = GetSecs;

iFrame = 1;
VBLTimestamp = 0;

if Settings.Display.ShowInSlowMotion

    when = Settings.Display.SlowMotionTime;

else
    
    when = 0;
    
end

if Settings.Movie.Record
    writerObj = VideoWriter('newfile.avi');
    writerObj.FrameRate = Settings.Movie.FrameRate;
    open(writerObj);
end

for iFrame=1:size(ItemX,2)

    DrawBalls;

    VBLTimestamp = Screen('Flip', ptbScreen.WinID, VBLTimestamp + when );
    
    if Settings.Movie.Record
        
        frame = Screen('GetImage', ptbScreen.WinID);
        writeVideo(writerObj,frame);
        
    end

    [keyIsDown, firstKeyPressTimes, ~, ~, ~]=PsychHID('KbQueueCheck');

    if keyIsDown

        break

    end

end

if Settings.Movie.Record; close(writerObj); end

endTime = GetSecs - startTime;

%% closing PTB
ListenChar(0); %% re-enable matlab to see key presses
if (exist('ptbScreen', 'var'))
  ptbScreen.Close;
  ShowCursor;
end;

disp(strcat('Total Time:',num2str(endTime)));
