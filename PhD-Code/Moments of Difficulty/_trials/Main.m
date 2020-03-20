clear all;
CaughtException= [];

%% general settings 
Folder.Settings= 'Settings/'; % path to settings folder 
Folder.Results= 'Results/';   % path to results folder

%% asking for user name
% ObserverName= 'testAngle'; %has to be the same as the name of the setting-file (.exp) 
if (~exist('ObserverName', 'var') || isempty(ObserverName))
  ObserverName= input('Please enter Settings file: ', 's');
end;

%% initializing book-keeping utilities
try 
  Settings= CExperimentalSettings(ObserverName, Folder.Settings);
catch exception
  if (~strcmp(exception.identifier, 'CExperimentalSettings:ErrorReadingFile'))
    rethrow(exception);
  end;
end;

%%
try
%% generating log file name
  if (Settings.Eyelink.Use) %% || Settings.LogEachSessionSeparately)
    Folder.Results= 'Results.Eye/';   % path to results folder
    CurrentTime= clock;
    LogFileName= [MOTFileName sprintf('_%d_%02d_%02d_%02d_%02d', CurrentTime(1:5))]
  else
    Folder.Results= 'Results/';   % path to results folder
    LogFileName= ObserverName;
  end;
  Log= CLog(LogFileName, Folder.Results);
  
   %% initializing eye tracker (if necessary)
  if (Settings.Eyelink.Use)
    if (EyelinkInit()~= 1)
     fprintf('Failed to initialize Eyelink\n');
      return;
    end;
  end;

  %% initializing PTB screen 
  ptbScreen= CPTBScreen(Settings.Screen);
  if (Settings.Screen.FullScreen)
    HideCursor;
    ListenChar(2); %% cut-off matlab from key presses
  end;
  
  %% setting eyetracker
  if (Settings.Eyelink.Use)
    EyelinkWindowHandle= EyelinkInitDefaults(ptbScreen.WinID);
    Eyelink('Command', Settings.Eyelink.FileEventFilter);
    Eyelink('Command', Settings.Eyelink.LinkSampleData);
    EyelinkEDFFile= [ObserverName '.edf'];
    LocalEDFFile= [Folder.Results LogFileName '.edf'];
    Eyelink('Openfile', EyelinkEDFFile);
    EyelinkDataSaved= false;
    Eyelink_TrialsBeforeCalibration= 0;
  end;    
  
%% INITIALIZE USB DEVICE
%% creating a USB device input queue
PsychHID('KbQueueCreate', Settings.fMRI.hidDeviceID);
PsychHID('KbQueueStart', Settings.fMRI.hidDeviceID);
  
  %% going through blocks
  for iBlock= 1:Settings.BlockN, 
      
    iTrial = 0;  
      
    %% preparing block
    PrepareBlock;
 
    %% running block
    for iTrial= 1:Settings.Trials*nCondition*nSpeeds
        
      dataMOT = load(strcat(char(all_trials_w_conditions{iTrial,2}),'-',int2str(Settings.MOT.Segments.Duration),'-',Settings.MOT.Segments.DurationStr,'-',int2str(all_trials_w_conditions{iTrial,1})));  

      MOT = dataMOT.MOT;
      
      TargetN = Settings.MOT.TargetN;
      Segments = Settings.MOT.Segments;

      Settings.MOT = dataMOT.MOT.Settings.MOT;
      
      Settings.MOT.TargetN = TargetN;
      Settings.MOT.Segments = Segments;
      
      RunTrial;
      
    end;
    
    %% saving collected data
    Log.Save(Settings.Current);
    
  end;
catch exception
  if (~strcmp(exception.identifier, 'Experiment:Abort'))
    CaughtException= exception;
  end;
end;

%% FINISH EXPERIMENT MSG
%% saying goodbye
Strings= {Settings.Message.ExperimentFinish};
ptbScreen.ShowCenteredMessage(Strings, [255 255 255]);
Screen('Flip', ptbScreen.WinID);
pause(2);

%% closing PTB
ListenChar(0); %% re-enable matlab to see key presses
if (exist('ptbScreen', 'var'))
  ptbScreen.Close;
  ShowCursor;
end;

%% saving eyelink data and getting file to local location
if (Settings.Eyelink.Use)
  % download data file 
  try
    Eyelink('CloseFile');
    fprintf('Receiving data file ''%s''\n', LocalEDFFile);
    status=Eyelink('ReceiveFile', [], LocalEDFFile);
    if status > 0
      fprintf('ReceiveFile status %d\n', status);
    end
    if 2==exist(LocalEDFFile, 'file')
      fprintf('Data file ''%s'' can be found in ''%s''\n', LocalEDFFile, pwd );
    end
  catch
    fprintf('Problem receiving data file ''%s''\n', LocalEDFFile);
  end
  Eyelink('ShutDown');
end;

%% rethrowing exception if we have one
if (~isempty(CaughtException))
  rethrow(CaughtException);
end