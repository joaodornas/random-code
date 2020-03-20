function h = movie2gray(video_name)

videoin = mmreader(video_name);

nFrames = videoin.NumberOfFrames;

writerobj = VideoWriter('grayscale');
writerobj.FrameRate = videoin.FrameRate;
open(writerobj);

for k = 1 : nFrames
    
    frame = read(videoin,k);
    
    framegray = rgb2gray(frame);
    
    writeVideo(writerobj,framegray);
    
end

close(writerobj);




end

