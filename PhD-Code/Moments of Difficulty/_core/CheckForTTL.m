function [ttlGotIt, ttlOnThisTime, ttlTime, clockTime]= CheckForTTL(hidID, tSmall, tCapital, ttlOn, AbortKey)
  
  [keyIsDown, firstKeyPressTimes, ~, ~, ~]=PsychHID('KbQueueCheck', hidID);
  
  ttlTime = firstKeyPressTimes(tSmall);
  clockTime = clock;
  
  ttlGotIt= false;
  ttlOnThisTime= false;
  %ttlTime= nan;
  if (keyIsDown)
    if (firstKeyPressTimes(AbortKey)~=0)
      throw(MException('Experiment:Abort', 'Aborted by user'))
    end;
    if (firstKeyPressTimes(tSmall)~=0 || firstKeyPressTimes(tCapital)~=0)
      if (ttlOn==false)
        ttlGotIt= true;
      end;
      ttlOnThisTime= true;
    end;
  end;

    