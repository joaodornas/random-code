function d = distanceXY(X,Y)


l = length(X);
 

if l >= 2
    
    [minX, minXind] = min(X);
    [maxX, maxXind] = max(X);

    [minY, minYind] = min(Y);
    [maxY, maxYind] = max(Y);

    if abs(maxX - minX) > abs(maxY - minY)

        X1 = minX;
        X2 = maxX;

        Y1 = Y(minXind);
        Y2 = Y(maxXind);

    elseif abs(maxX - minX) < abs(maxY - minY)

        Y1 = minY;
        Y2 = maxY;

        X1 = X(minYind);
        X2 = X(maxYind);

    end

    d = sqrt((X1 - X2)^2 + (Y1 - Y2)^2);

else
    
    d = 1;
    
end

end

