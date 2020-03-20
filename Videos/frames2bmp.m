function frames2bmp(video_name,X,Y,radius)

videoin = mmreader(video_name);

nFrames = 300;
vidHeight = videoin.Height;
vidWidth = videoin.Width;

% Preallocate movie structure.
mov(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
           'colormap', []);

     videoout = VideoWriter('video-RF');
     
     videoout.FrameRate = 30;
     
     open(videoout);

% Read one frame at a time.
for k = 1 : nFrames
    
    mov(k).cdata = read(videoin, k);
    
    if k < 10
       frameindex = strcat('00', int2str(k));
    elseif k < 100
       frameindex = strcat('0', int2str(k));
    else
       frameindex = int2str(k);
    end

    yellow = uint8([255 255 0]);
    %circles = int32([X Y radius;X Y radius-1;X Y radius-2;X Y radius-3;X Y radius-4;X Y radius-5;X Y radius-6;X Y radius-7;X Y radius-8;X Y radius-9;X Y radius-10;X Y radius-11;X Y radius-12;X Y radius-13;X Y radius-14]);
    circles = int32([X Y radius;X Y radius-1;X Y radius-2;X Y radius-3;X Y radius-4;X Y radius-5;X Y radius-6;X Y radius-7]); %;X Y radius-8;X Y radius-9;X Y radius-10;X Y radius-11;X Y radius-12;X Y radius-13;X Y radius-14]);
    shapeInserter = vision.ShapeInserter('Shape','Circles','BorderColor','Custom','CustomBorderColor',yellow);
    frame = step(shapeInserter, mov(k).cdata, circles);
    %imwrite(frame,['frame-' frameindex '.jpg']);

    writeVideo(videoout, frame);
    
end

close(videoout);

end


