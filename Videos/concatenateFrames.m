function concatenateFrames

nFrames = 300;

movieFrames(1:nFrames) = struct('cdata',zeros(720,1024));

for k=1:nFrames
        
        if k < 10
            
            frame = strcat('00',int2str(k));
            
        elseif k < 100
            
            frame = strcat('0',int2str(k));
            
        else
            
            frame = int2str(k);
            
        end
            
        img = imread(strcat('frame-',frame,'.jpg'));
        
        movieFrames(k).cdata = img;
        
end
    


writerobj = VideoWriter('FinalVideo');
writerobj.FrameRate = 30;
open(writerobj);

for i=1:nFrames

    writeVideo(writerobj, movieFrames(i).cdata);
    
end

close(writerobj);


end

