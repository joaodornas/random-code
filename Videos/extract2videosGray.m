function h = extract2videosGray(video_name,start,nFrames)

videoin = mmreader(video_name);

writerobj720 = VideoWriter('8bit-720x480.mp4');
writerobj720.FrameRate = videoin.FrameRate;
open(writerobj720);

writerobj1024 = VideoWriter('8bit-1024x768.mp4');
writerobj1024.FrameRate = videoin.FrameRate;
open(writerobj1024);

for k=start:(nFrames + start)
    
   framergb = read(videoin,k);
   
   framegray = rgb2gray(framergb);
   
   framegrayscaled = imresize(framegray, [768 1024]);
   
   writeVideo(writerobj720,framegray);
   
   writeVideo(writerobj1024,framegrayscaled);
    
end

close(writerobj720);

close(writerobj1024);


end

