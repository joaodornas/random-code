
nTR = 331;

for iTR=1:nTR
   
    iiTR = iTR - 1;
    
    if iiTR < 10
        
        idx = strcat('000',int2str(iiTR));
    
    elseif (iiTR >= 10) && (iiTR < 100)
        
        idx = strcat('00',int2str(iiTR));
        
    elseif iiTR >= 100
        
        idx = strcat('0',int2str(iiTR));
        
    end
    
    mat(iTR).mc = load(strcat('MAT_',idx));
    
end

save('mat.mat','mat');