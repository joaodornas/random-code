

%% Draw Each Ball

H = [ItemX(thisBall, iFrame)+CenterXY(1) - Width/2,...
    ItemY(thisBall, iFrame)+CenterXY(2) - Width/2,...
    ItemX(thisBall, iFrame)+CenterXY(1) + Width/2,...
    ItemY(thisBall, iFrame)+CenterXY(2) + Width/2];

Screen('FrameRect', ptbScreen.WinID, Color, H');
