
DataFile= input('Please enter the AnalyzeTrajectories datafile name: ', 's');

load(DataFile);

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

ItemX = MOT.xi;
ItemY = MOT.yi;
BallColors = MOT.BallColors;

Settings.MOT = MOT.Settings.MOT;
side = Settings.MOT.Area.Size(1)/2;

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

Color = [255 0 0];
hs = 1;
WindowPace = nWindow*sFrames;

for iFrame=1:size(ItemX,2)
    
    if (iFrame > 1) && (mod(iFrame,WindowPace) == 1); hs = hs + 1; end 

    DrawBalls;
    
    if Settings.Difficulty.Color.Entropy.Show
    
        if (iFrame >= (kECmin-1)*sFrames+1) && (iFrame <= (kECmin+nWindow-2)*sFrames+1) 
            
            Screen('DrawText',ptbScreen.WinID,'COLOR:min.Entropy',CenterXY(1)+side,CenterXY(2)+side/2,Color,1);
                   
        end
        
        if (iFrame >= (kECmax-1)*sFrames+1) && (iFrame <= (kECmax+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'COLOR:max.Entropy',CenterXY(1)+side,CenterXY(2)+side/2,Color,1);
            
        end
    
    end
    
    if Settings.Difficulty.Color.Weights.Show
    
        if (iFrame >= (kVCmin-1)*sFrames+1) && (iFrame <= (kVCmin+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'COLOR:min.Weigth',CenterXY(1)+side,CenterXY(2)+side/2,Color,1);
            
        end
        
        if (iFrame >= (kVCmax-1)*sFrames+1) && (iFrame <= (kVCmax+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'COLOR:max.Weigth',CenterXY(1)+side,CenterXY(2)+side/2,Color,1);
            
        end
    
    end
    
    if Settings.Difficulty.Movement.Entropy.Show

        if (iFrame >= (kEMmin-1)*sFrames+1) && (iFrame <= (kEMmin+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'MOVEMENT:min.Entropy',CenterXY(1)+side,CenterXY(2)-side/2,Color,1);
            
        end
        
        if (iFrame >= (kEMmax-1)*sFrames+1) && (iFrame <= (kEMmax+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'MOVEMENT:max.Entropy',CenterXY(1)+side,CenterXY(2)-side/2,Color,1);
            
        end
    
    end

    if Settings.Difficulty.Movement.Weights.Show

        if (iFrame >= (kVMmin-1)*sFrames+1) && (iFrame <= (kVMmin+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'MOVEMENT:min.Weigth',CenterXY(1)+side,CenterXY(2)-side/2,Color,1);
            
        end
        
        if (iFrame >= (kVMmax-1)*sFrames+1) && (iFrame <= (kVMmax+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'MOVEMENT:max.Weigth',CenterXY(1)+side,CenterXY(2)-side/2,Color,1);
            
        end
    
    end

    if Settings.Difficulty.Bounce.Entropy.Show

        if (iFrame >= (kEBmin-1)*sFrames+1) && (iFrame <= (kEBmin+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'BOUNCE:min.Entropy',CenterXY(1)-2*side,CenterXY(2)+side/2,Color,1);
            
        end
        
        if (iFrame >= (kEBmax-1)*sFrames+1) && (iFrame <= (kEBmax+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'BOUNCE:max.Entropy',CenterXY(1)-2*side,CenterXY(2)+side/2,Color,1);
            
        end
    
    end

    if Settings.Difficulty.Bounce.Weights.Show

        if (iFrame >= (kVBmin-1)*sFrames+1) && (iFrame <= (kVBmin+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'BOUNCE:min.Weight',CenterXY(1)-2*side,CenterXY(2)+side/2,Color,1);
            
        end
        
        if (iFrame >= (kVBmax-1)*sFrames+1) && (iFrame <= (kVBmax+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'BOUNCE:max.Weight',CenterXY(1)-2*side,CenterXY(2)+side/2,Color,1);
            
        end
    
    end

    if Settings.Difficulty.Quadrant.Entropy.Show

        if (iFrame >= (kEQmin-1)*sFrames+1) && (iFrame <= (kEQmin+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'QUADRANT:min.Entropy',CenterXY(1)-2*side,CenterXY(2)-side/2,Color,1);
            
        end
        
        if (iFrame >= (kEQmax-1)*sFrames+1) && (iFrame <= (kEQmax+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'QUADRANT:max.Entropy',CenterXY(1)-2*side,CenterXY(2)-side/2,Color,1);

        end
        
    end

    if Settings.Difficulty.Quadrant.Weights.Show
    
        if (iFrame >= (kVQmin-1)*sFrames+1) && (iFrame <= (kVQmin+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'QUADRANT:min.Weigth',CenterXY(1)-2*side,CenterXY(2)-side/2,Color,1);
            
        end
        
        if (iFrame >= (kVQmax-1)*sFrames+1) && (iFrame <= (kVQmax+nWindow-2)*sFrames+1)
            
            Screen('DrawText',ptbScreen.WinID,'QUADRANT:max.Weigth',CenterXY(1)-2*side,CenterXY(2)-side/2,Color,1);
            
        end
    
    end
    
    Screen('DrawText',ptbScreen.WinID,strcat(int2str(hs)),CenterXY(1)+side, CenterXY(2)+side,Color,1);

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
