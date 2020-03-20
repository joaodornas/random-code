

%% Draw Each Ball
for iN=1:Settings.MOT.TotalN

  H= [ItemX(iN, iFrame)+CenterXY(1)-Settings.MOT.Target.Diameter/2,...
      ItemY(iN, iFrame)+CenterXY(2)-Settings.MOT.Target.Diameter/2,...
      ItemX(iN, iFrame)+CenterXY(1)+Settings.MOT.Target.Diameter/2,...
      ItemY(iN, iFrame)+CenterXY(2)+Settings.MOT.Target.Diameter/2];

  Screen('FillOval', ptbScreen.WinID, BallColors(iN,:), H');

end
