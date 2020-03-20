function h = subtractPixelCondition(nConditions,nBmp,pixel)


for k=1:nConditions
    
    for i=1:(nBmp+1)
    
        img = imread(strcat('Condition',int2str(k),'-',int2str(i-1),'.bmp'));
        imgsub = img - pixel;
        imgsub(imgsub<0) = 1;
        
        imwrite(imgsub,strcat('Condition',int2str(k),'-',int2str(i-1),'.bmp'));
    end
    
end



end

