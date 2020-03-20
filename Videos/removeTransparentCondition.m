function h = removeTransparentCondition(nConditions,nBmp)


for k=1:nConditions
    
    for i=1:(nBmp+1)
    
        img = imread(strcat('Condition',int2str(k),'-',int2str(i-1),'.bmp'));

        img(img==1) = 0;
        
        imwrite(img,strcat('Condition',int2str(k),'-',int2str(i-1),'.bmp'));
    
    end
    
end


end

