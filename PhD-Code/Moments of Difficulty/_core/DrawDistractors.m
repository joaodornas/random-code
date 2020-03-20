

%% Draw Each Ball

for b=1:length(distractors)

    H = [ItemX(distractors(b), iFrame)+CenterXY(1) - Width/2,...
        ItemY(distractors(b), iFrame)+CenterXY(2) - Width/2,...
        ItemX(distractors(b), iFrame)+CenterXY(1) + Width/2,...
        ItemY(distractors(b), iFrame)+CenterXY(2) + Width/2];

    Screen('FrameOval', ptbScreen.WinID, Color, H');

end