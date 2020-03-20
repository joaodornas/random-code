function h = videocrop(video_name,line_interval,column_interval)


videoin = mmreader(video_name);

nFrames = videoin.NumberOfFrames;

writerobj = VideoWriter('cropvideo');
writerobj.FrameRate = videoin.FrameRate;
open(writerobj);

for i=1:nFrames
    
    framergb = read(videoin,i);
    
    framecrop = framergb(line_interval,column_interval);
    
    writeVideo(writerobj, framecrop);
    
end

close(writerobj);



end