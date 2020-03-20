
%% SAVE TTL INFORMATION ON LOG

clockTimeNow = clock;

Log.Block{end}.ttlTime{iTTL,1} = iTrial;
Log.Block{end}.ttlTime{iTTL,2} = TTLLogLabel;
Log.Block{end}.ttlTime{iTTL,3} = ttlTime;
Log.Block{end}.ttlTime{iTTL,4} = clockTimeNow;

if iTTL == 1
    
    Log.Block{end}.ttlTime{iTTL,5} = 0;
    
elseif iTTL > 1
    
    Log.Block{end}.ttlTime{iTTL,5} = Log.Block{end}.ttlTime{iTTL,3} - Log.Block{end}.ttlTime{iTTL-1,3};
    
end

iTTL= iTTL+1;
