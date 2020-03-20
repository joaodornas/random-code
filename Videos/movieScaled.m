function h = movieScaled(video_name,width,height)


videoin = mmreader(video_name);

nFrames = videoin.numberOfFrames;

writerobj = VideoWriter('new_scale');
writerobj.FrameRate = videoin.FrameRate;
open(writerobj);

for k = 1 : nFrames
    
    frame = read(videoin,k);
    
    frame = imresize(frame,[height width]);
    
    writeVideo(writerobj, frame);
    
end

close(writerobj);



end

