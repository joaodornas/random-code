
if iTR < 10
        
    iTR_number = strcat('0000',int2str(iTR));

elseif iTR >= 10 && iTR < 100

    iTR_number = strcat('000',int2str(iTR));

else

    iTR_number = strcat('00',int2str(iTR));

end

if new_iTR < 10
        
    new_iTR_number = strcat('0000',int2str(new_iTR));

elseif new_iTR >= 10 && new_iTR < 100

    new_iTR_number = strcat('000',int2str(new_iTR));

elseif new_iTR >= 100 && new_iTR < 1000

    new_iTR_number = strcat('00',int2str(new_iTR));
    
else
    
    new_iTR_number = strcat('0',int2str(new_iTR));

end

if iTR_FSL < 10
        
    iTR_number_FSL = strcat('000',int2str(iTR_FSL));

elseif iTR_FSL >= 10 && iTR_FSL < 100

    iTR_number_FSL = strcat('00',int2str(iTR_FSL));

else

    iTR_number_FSL = strcat('0',int2str(iTR_FSL));

end