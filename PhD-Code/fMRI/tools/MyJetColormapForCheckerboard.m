function cmap = MyJetColormapForCheckerboard(applyThreshold)

if applyThreshold == 1 || applyThreshold == 2
    
    cmap = zeros(64,3);
    
    cmap(1,1) = 1; cmap(1,2) = 1; cmap(1,3) = 1;   % leave first entry white!
    cmap(2,1) = 0.9; cmap(2,2) = 0.9; cmap(2,3) = 0.9;   % leave second entry light grey!

    cmap(3:8,1)   = 0;
    cmap(3:8,2)   = 0;
    cmap(3:8,3)   = (8 + (3:8)) / 16.0;    % dark blue to blue

    cmap(9:24,1)  = 0;
    cmap(9:24,2)  = (1:16) / 16.0;        %  blue to cyan
    cmap(9:24,3)  = 1;

    cmap(25:40,1) = (1:16) / 16;
    cmap(25:40,2) = 1;                    %  cyan to yellow
    cmap(25:40,3) = (15:-1:0) / 16;

    cmap(41:56,1) = 1;
    cmap(41:56,2) = (15:-1:0) / 16;       %  yellow to red
    cmap(41:56,3) = 0;

    cmap(57:64,1) = (15:-1:8) / 16;
    cmap(57:64,2) = 0;                    %  red to dark red
    cmap(57:64,3) = 0;

else
    
    cmap = colormap('jet');

    cmap(32,1) = 1; cmap(32,2) = 1; cmap(32,3) = 1;   % leave first entry white!
    cmap(33,1) = 0.9; cmap(33,2) = 0.9; cmap(33,3) = 0.9;   % leave second entry light grey!
    
end

return;