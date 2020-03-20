function h = vectorsSameBMP(width,height,background,nFrames)

bitmap = ones(height,width,'uint8');

bitmap = bitmap.*background;

k = 1;

while (k < (nFrames + 1))

    X = reshape(bitmap.',[],1);

    X = X.';
    
    if k < 10
        frameindex = strcat('000', int2str(k));
    elseif k < 100
        frameindex = strcat('00', int2str(k));
    elseif k < 1000
        frameindex = strcat('0', int2str(k));
    else
        frameindex = int2str(k);
    end

    path = strcat('bitmap-',int2str(width),'-',int2str(background),'_', frameindex);

    dlmwrite(strcat(path, '.txt'), X, 'delimiter', '\t', '-append', 'newline', 'pc');
    
    %imwrite(bitmap,strcat(path,'.bmp'));
    
    k = k + 1;

end


end