
%% general settings 
Folder.Settings= 'Settings/'; % path to settings folder 
Folder.Events = 'Events/';

DataFile= input('Please enter the Paramters file name: ', 's');

load(strcat(Folder.Events,DataFile));

SettingsFile= input('Please enter the Settings file name: ', 's');

SettingsEvents = CExperimentalSettings(SettingsFile, Folder.Settings);

Settings.Events = SettingsEvents.Events;
Settings.Screen = SettingsEvents.Screen;

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
hs = 0;
VBLTimestamp = 0;
for hs=allHSEvents
            
    for iFrame=((hs-1)*FrameRate/2 + 1):(hs*FrameRate/2)
        
        if Settings.Events.ShowEventsAtAll && (iFrame/FrameRate <= Settings.MOT.Duration)

            if (iFrame == 1) || (mod(iFrame,FrameRate/2) == 1)

                startShowArea = true;

                startShowDisplacement = eventsBallDistance(1:nBalls,hs);

            end

            FontSize = 20;
            Screen('TextSize', ptbScreen.WinID, FontSize);

            for b=1:nBalls

               thisBall = selectedBalls(b);

               if eventsArea(hs) && Settings.Events.AreaEvents.Show

                   if startShowArea
                       minArea = min(areaBalls(iFrame:(iFrame+FrameRate/2)));
                       startShowArea = false;
                   end

                   Color = Settings.Events.Thresholds.Color;
                   Screen('DrawText',ptbScreen.WinID,strcat(num2str(round((minArea/maxArea)*100)),'%'),CenterXY(1),CenterXY(2)-side-30,Color,1);
                   Color = Settings.Events.AreaEvents.Cue.Color;
                   Width = Settings.Events.AreaEvents.Cue.Width;
                   Screen('DrawText',ptbScreen.WinID,'Area Event',(CenterXY(1)-side), (CenterXY(2)-side)-30,Color,1);
                   DrawCues;

               end

               if eventsBallDistance(b,hs) && Settings.Events.DisplacementEvents.Show

                   if startShowDisplacement(b) 
                       startPoint(1,b) = xi(selectedBalls(b),iFrame) + CenterXY(1);
                       endPoint(1,b) = xi(selectedBalls(b),iFrame+FrameRate/2) + CenterXY(1); 
                       startPoint(2,b) = yi(selectedBalls(b),iFrame) + CenterXY(2);
                       endPoint(2,b) = yi(selectedBalls(b),iFrame+FrameRate/2) + CenterXY(2); 
                       startShowDisplacement(b) = 0;
                   end

                   Color = Settings.Events.DisplacementEvents.Cue.Color;
                   Screen('DrawText',ptbScreen.WinID,'Displacement Event',(CenterXY(1)-side), (CenterXY(2)-side) - 60,Color,1);
                   Color = Settings.Events.DisplacementEvents.Trajectory.Color;
                   Line_Width = Settings.Events.DisplacementEvents.Trajectory.Width;
                   Screen('DrawLine',ptbScreen.WinID,Color,startPoint(1,b),startPoint(2,b),endPoint(1,b),endPoint(2,b),Line_Width);
                   Color = Settings.Events.Thresholds.Color;
                   Screen('DrawText',ptbScreen.WinID,strcat(num2str(round((displacement(b,hs)/(Settings.MOT.Area.Size(1)))*100)),'%'),CenterXY(1)+b*100,CenterXY(2)-side-60,Color,1);
                   Color = Settings.Events.DisplacementEvents.Cue.Color;
                   Width = Settings.Events.DisplacementEvents.Cue.Width;
                   DrawCues;

               end

               if eventsBallColorMatch(b,hs) && Settings.Events.ColorMatchEvents.Show
                  Color = Settings.Events.ColorMatchEvents.Cue.Color;
                  Width = Settings.Events.ColorMatchEvents.Cue.Width;
                  Screen('DrawText',ptbScreen.WinID,'Color Match Event',(CenterXY(1)-side), (CenterXY(2)-side) - 90,Color,1);
                  DrawCues;
               end 

            end
            
            if Settings.Events.MixingAllBalls; mixingBalls = 1:size(eventsBallMixing,1); else mixingBalls=selectedBalls; end
            
            for b=mixingBalls
                
                thisBall = b;
            
                 if sum(eventsBallMixing(:,b,hs))~=0 && Settings.Events.MixingEvents.Show

                      idx_mixedBalls = find(eventsBallMixing(:,b,hs));

                      for k=1:length(idx_mixedBalls)

                          if Settings.Events.MixingAllBalls && (Settings.Events.MixingEvents.TopEvents ~= 0)

                              if sweep(idx_mixedBalls(k),b,hs) == max_hs_sorted(find(max_hs_idx == hs))

                                  doIt = true;

                              else

                                  doIt = false;

                              end

                          else

                              doIt = true;

                          end

                          if doIt

                              Color = Settings.Events.MixingEvents.Cue.Color;
                              Width = Settings.Events.MixingEvents.Cue.Width;
                              Screen('DrawText',ptbScreen.WinID,'Mixing Event',(CenterXY(1)-side), (CenterXY(2)-side) - 120,Color,1);
                              DrawCues;

                              startPoint = [(xi(b,iFrame) + CenterXY(1)) (yi(b,iFrame) + CenterXY(2))];
                              endPoint = [(xi(idx_mixedBalls(k),iFrame) + CenterXY(1)) (yi(idx_mixedBalls(k),iFrame) + CenterXY(2))];
                              Color = Settings.Events.MixingEvents.Link.Color;
                              Line_Width = Settings.Events.MixingEvents.Link.Width;
                              Screen('DrawLine',ptbScreen.WinID,Color,startPoint(1),startPoint(2),endPoint(1),endPoint(2),Line_Width);

                          end

                      end

                 end 

            end
            
        end
        
        if Settings.Events.Distractor.Show

            thisBall = distractor;
            Color = Settings.Events.Distractor.Cue.Color;
            Width = Settings.Events.Distractor.Cue.Width;
            DrawDistractors;
            
        end

        DrawBalls;

        if Settings.Events.ColorMatchEvents.Show || Settings.Events.MixingEvents.Show
            
            if Settings.Events.Distractor.Show
                
                balls = [selectedBalls distractor];
                
            else
                
                balls = selectedBalls;
                
            end
            
            DrawRepulsionRadius;
            
        end

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

        if keyIsDown || (iFrame/FrameRate > Settings.MOT.Duration)

            break

        end
           
    end
    
    if (escapeKeyTime ~= 0)
            
       break
            
    end
    
end

% n = 1;
% for i=1/2:1/2:Settings.MOT.Duration
%    
%     hs_frames_duration(n) = Flip(i*FrameRate + 1,2) - Flip((i-1/2)*FrameRate + 1,2);
%     
%     n = n + 1;
%     
% end

endTime = GetSecs - startTime;

%% closing PTB
ListenChar(0); %% re-enable matlab to see key presses
if (exist('ptbScreen', 'var'))
  ptbScreen.Close;
  ShowCursor;
end;

disp(strcat('Total Time:',num2str(endTime)));
