function h = rfmapstimulusRec(width,height,side,spot_radius,background)

%nCentersX = width / (2 * spot_radius);
%nCentersY = height / (2 * spot_radius);

nCentersX = side / (2 * spot_radius);
nCentersY = side / (2 * spot_radius);

nMaps = 2 * nCentersX * nCentersY;

centers = zeros(nMaps,3);

c = 1;

for j=1:2
    
    for k=1:nCentersX

        for i=1:nCentersY

            centers(c,1) = k * 2 * spot_radius - spot_radius + (width - side)/2;

            centers(c,2) = i * 2 * spot_radius - spot_radius + (height - side)/2;

            centers(c,3) = j;

            c = c + 1;

        end

    end

end

stimulusMap(1:nMaps) = struct('map',ones(height,width,'uint8'));

for i=1:nMaps
    
    stimulusMap(i).map = stimulusMap(i).map.*background;
 
end

randOrder = randperm(nMaps);

randCenters  = centers(randOrder,:);

for i=1:nMaps
    
   if randCenters(i,3) == 1
       
       color = 1;
       
   else
       
       color = 255;
       
   end
   
   x1 = randCenters(i,1) - spot_radius;
   x2 = randCenters(i,1) + spot_radius;
   
   y1 = randCenters(i,2) - spot_radius;
   y2 = randCenters(i,2) + spot_radius;  
   
   if (x1 > 0) && (y1 > 0) 
       
       stimulusMap(i).map(y1:y2,x1:x2) = color;
      
   end
   
end

stimulusMapReverse = fliplr(stimulusMap);

k = 1;

while (k < (nMaps + 1))

    X = reshape(stimulusMap(k).map.',[],1);

    X = X.';
    
    XR = reshape(stimulusMapReverse(k).map.',[],1);

    XR = XR.';

    if k < 10
        frameindex = strcat('000', int2str(k));
    elseif k < 100
        frameindex = strcat('00', int2str(k));
    elseif k < 1000
        frameindex = strcat('0', int2str(k));
    else
        frameindex = int2str(k);
    end

    path = strcat('spot-a-forward-',int2str(width),'-',int2str(spot_radius),'-',int2str(background),'_', frameindex);

    dlmwrite(strcat(path, '.txt'), X, 'delimiter', '\t', '-append', 'newline', 'pc');
    
    %imwrite(stimulusMap(k).map,strcat(path,'.bmp'));
    
    pathr = strcat('spot-b-backward-',int2str(width),'-',int2str(spot_radius),'-',int2str(background),'_', frameindex);

    dlmwrite(strcat(pathr, '.txt'), XR, 'delimiter', '\t', '-append', 'newline', 'pc');
    
    %imwrite(stimulusMapReverse(k).map,strcat(pathr,'.bmp'));

    k = k + 1;

end


end