
[keyIsDown, firstKeyPressTimes, ~, ~, ~]=PsychHID('KbQueueCheck', Settings.fMRI.hidDeviceID);

EscapeKeyTime = firstKeyPressTimes(Settings.Keyboard.Escape) ~= 0;
SpaceKeyTime = firstKeyPressTimes(Settings.Keyboard.Space) ~= 0;
EnterKeyTime = firstKeyPressTimes(Settings.Keyboard.Enter) ~= 0;
LeftKeyTime = firstKeyPressTimes(Settings.Keyboard.Left) ~= 0;
RighKeyTime = firstKeyPressTimes(Settings.Keyboard.Right) ~= 0;
UpKeyTime = firstKeyPressTimes(Settings.Keyboard.Up) ~= 0;
DownKeyTime = firstKeyPressTimes(Settings.Keyboard.Down) ~= 0;
tSmallKeyTime = firstKeyPressTimes(Settings.Keyboard.tSmall) ~= 0;
tCapitalKeyTime = firstKeyPressTimes(Settings.Keyboard.tCapital) ~= 0;
ResponseStick_TopKeyTime = firstKeyPressTimes(Settings.Keyboard.ResponseStick_Top) ~= 0;
ResponseStick_YellowKeyTime = firstKeyPressTimes(Settings.Keyboard.ResponseStick_Yellow) ~= 0;
ResponseStick_GreenKeyTime = firstKeyPressTimes(Settings.Keyboard.ResponseStick_Green) ~= 0;
ResponseStick_RedKeyTime = firstKeyPressTimes(Settings.Keyboard.ResponseStick_Red) ~= 0;
ResponseStick_BlueKeyTime = firstKeyPressTimes(Settings.Keyboard.ResponseStick_Blue) ~= 0;

if tSmallKeyTime
    
    ttlTime = firstKeyPressTimes(Settings.Keyboard.tSmall);
    
elseif tCapitalKeyTime
   
    ttlTime = firstKeyPressTimes(Settings.Keyboard.tCapital);
    
end

if (EscapeKeyTime)
                  
	throw(MException('Experiment:Abort', 'Aborted'));
    
end