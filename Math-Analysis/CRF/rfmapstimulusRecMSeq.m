function h = rfmapstimulusRecMSeq(width,height,side,spot_radius,background,powerVal,CC)

nCentersX = side;
nCentersY = side;

N = 2^powerVal - 1;
p_steps = 2^powerVal / (nCentersX * nCentersY);

stimulusMap = ones(height,width,'uint8');

             writerObj = VideoWriter(strcat('stimulusMap-MSeq-side-',int2str(side),'-spot_radius-',int2str(spot_radius),'-powerVal-',int2str(powerVal),'.avi'));
             writerObj.FrameRate = 30;
             open(writerObj);

msequence = mseq(2,powerVal);

for m=1:N
    
    stimulusMap(1:side*CC,1:side*CC) = background;
    
    for k=1:nCentersX

      for i=1:nCentersY
          
           D = p_steps*k + p_steps*i*nCentersY;
           
           %msequence = mseq(2,powerVal,D + m - 1);
           
           index = D + m - 1;
           
           if index > size(msequence,1)
              
               index = mod(index,size(msequence,1));
               
               if index == 0
                   
                   index = size(msequence,1);
                   
               end
               
           end
           
           if msequence(index) == 1

               color = 255;

           else

               color = 1;

           end

           Seq(k,i,m) = msequence(index);
           
           x1 = (k - 1) * 2 * spot_radius * CC + 1;
           x2 = ( (k - 1) * 2 * spot_radius + 2 * spot_radius ) * CC;

           y1 = (i - 1) * 2 * spot_radius * CC + 1;
           y2 = ( (i - 1) * 2 * spot_radius + 2 * spot_radius ) * CC;  

           if (x1 > 0) && (y1 > 0) 

               stimulusMap(y1:y2,x1:x2) = color;

           end
           
        end

    end
    
    writeVideo(writerObj,stimulusMap);
    
%     stimulusMapReverse = fliplr(stimulusMap);
    
    X = reshape(stimulusMap.',[],1);

    X = X.';

%     XR = reshape(stimulusMapReverse.',[],1);
% 
%     XR = XR.';

    if m < 10
        frameindex = strcat('0000', int2str(m));
    elseif m < 100
        frameindex = strcat('000', int2str(m));
    elseif m < 1000
        frameindex = strcat('00', int2str(m));
    elseif m < 10000
        frameindex = strcat('0', int2str(m));
    else
        frameindex = int2str(m);
    end

    path = strcat('spot-a-forward-w',int2str(width),'-s',int2str(side),'-r',int2str(spot_radius),'-b',int2str(background),'-p',int2str(powerVal),'_', frameindex);

    dlmwrite(strcat(path, '.txt'), X, 'delimiter', '\t', '-append', 'newline', 'pc');

    %imwrite(stimulusMap,strcat(path,'.bmp'));

%     pathr = strcat('spot-b-backward-w',int2str(width),'-s',int2str(side),'-r',int2str(spot_radius),'-b',int2str(background),'-p',int2str(powerVal),'_', frameindex);
% 
%     dlmwrite(strcat(pathr, '.txt'), XR, 'delimiter', '\t', '-append', 'newline', 'pc');

    %imwrite(stimulusMapReverse,strcat(pathr,'.bmp'));
    
end

close(writerObj);

save('Sequence','Seq');

end