videoObj = VideoWriter('entropies.avi');
FrameRate = 3;
open(videoObj);

load('entropies.mat');

nFrames = length(F);

for j=1:nFrames
    
    for i=1:FrameRate
        
        frame = F(j).cdata;

        writeVideo(videoObj,frame);  

    end
    
end

close(videoObj);