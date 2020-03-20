
%% SAVE FLIP INFORMATION ON LOG

Log.Block{end}.vblTime{iVBL,1} = FlipLogLabel;
Log.Block{end}.vblTime{iVBL,2} = VBLTimestamp;
Log.Block{end}.vblTime{iVBL,3} = StimulusOnsetTime;
Log.Block{end}.vblTime{iVBL,4} = FlipTimestamp;
Log.Block{end}.vblTime{iVBL,5} = Missed;
Log.Block{end}.vblTime{iVBL,6} = Beampos;

if iVBL == 1
    
    Log.Block{end}.vblTime{iVBL,7} = 0;
    Log.Block{end}.vblTime{iVBL,8} = 0;
    
elseif iVBL > 1
    
    Log.Block{end}.vblTime{iVBL,7} = Log.Block{end}.vblTime{iVBL,2} - Log.Block{end}.vblTime{iVBL-1,2};
    Log.Block{end}.vblTime{iVBL,8} = 1/Log.Block{end}.vblTime{iVBL,7};
    
end

iVBL = iVBL + 1;
