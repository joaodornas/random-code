function videoInsertYellowCircleCRF(video_name,nFrames,X,Y,invert)

radius = 68;

videoin = mmreader(video_name);

vidHeight = videoin.Height;
vidWidth = videoin.Width;
FrameRate = videoin.FrameRate;

mov(1:nFrames) = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),'colormap',[]);

for k = 1 : nFrames
    
    mov(k).cdata = read(videoin, k);
    
end

if invert == 1
    
   movinverted = fliplr(mov);
    
end

videoout = VideoWriter(strcat(video_name(1:end-4),'-yellow-mark-CRF.mp4'));
videoout.FrameRate = FrameRate;
open(videoout);

if invert == 1

    videooutinv = VideoWriter(strcat(video_name(1:end-4),'-yellow-mark-CRF-inverted.mp4'));
    videooutinv.FrameRate = FrameRate;
    open(videooutinv);
    
end

for k = 1 : nFrames

    yellow = uint8([255 255 0]);

    circles = int32([X Y radius;X Y radius-1;X Y radius-2;X Y radius-3;X Y radius-4;X Y radius-5;X Y radius-6;X Y radius-7]);

    shapeInserter = vision.ShapeInserter('Shape','Circles','BorderColor','Custom','CustomBorderColor',yellow);

    frame = step(shapeInserter,mov(k).cdata,circles);
    
    writeVideo(videoout,frame);
    
    if invert == 1
    
        frameinverted = step(shapeInserter,movinverted(k).cdata,circles);
    
        writeVideo(videooutinv,frameinverted);
        
    end
    
end

close(videoout);

if invert == 1
    
    close(videooutinv);
    
end

end