

%% Draw Each Ball

for b=1:length(balls)

  H= [ItemX(balls(b), iFrame)+CenterXY(1)-Settings.MOT.RepulsionRadius,...
      ItemY(balls(b), iFrame)+CenterXY(2)-Settings.MOT.RepulsionRadius,...
      ItemX(balls(b), iFrame)+CenterXY(1)+Settings.MOT.RepulsionRadius,...
      ItemY(balls(b), iFrame)+CenterXY(2)+Settings.MOT.RepulsionRadius];

  Screen('FrameOval', ptbScreen.WinID, Settings.Events.RepulsionRadius.Color, H');
  
end

