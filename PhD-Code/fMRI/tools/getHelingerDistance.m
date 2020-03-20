function distance = getHellingerDistance(dist1,dist2)

X = dist1;

Y = dist2;
        
XY = X.*Y;
sqrtXY = sqrt(XY);
        
distance = sqrt( 1 - sum( sqrtXY ) );
        
end 

