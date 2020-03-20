function h = bpm2movie(FrameRate,start,finnish)


writerobj = VideoWriter('framescolados');
writerobj.FrameRate = FrameRate;
open(writerobj);

for k = start : finnish
        
    if k < 10
       frameindex = strcat('00', int2str(k));
    elseif k < 100
       frameindex = strcat('0', int2str(k));
    else
       frameindex = int2str(k);
    end

    frame = imread(['frame-' frameindex '.jpg']);

    writeVideo(writerobj, frame);
    
end

close(writerobj);


end


