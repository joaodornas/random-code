function scalesZ = getscalesZ

scalesZ = [];


XC1 = double(1000);
YC1 = double(1000);
XC2 = double(200);
YC2 = double(200);

f = double(100);
sx = double(2);
sy = double(2);

dc = f * sqrt( (sx^2)*(XC1 - XC2)^2 + (sy^2)*(YC1 - YC2)^2 );

dia = 1;

for i=1:199
    
    if i+1 < 10
        index = strcat('000', int2str(i+1));
    elseif i+1 < 100
        index = strcat('00', int2str(i+1));
    elseif i+1 < 1000
        index = strcat('0', int2str(i+1));
    end

    imgP = imread(strcat(index,'.png'));
    
    [LP1, LP2, VP1] = corners(imgP,3);    
    
    [XP, YP] = getcorners(LP1, LP2, VP1);
    
    if length(XP) < 2
        
        di = dia;
        
    else
     
        di = distanceXY(XP,YP);
        
        if (di == 1 || di < (0.8*dia))
            
            di = dia;
            
        end
         
    end
        
    scalesZ(i) = dc/di;
    
    dia = di;
    
end

maxall = max(scalesZ);

scalesZ = scalesZ / maxall;

end

