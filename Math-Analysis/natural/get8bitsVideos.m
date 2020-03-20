function get8bitsVideos(v,ext,mainPath)


videoin = VideoReader(['_MATRIZ-v' int2str(v) '-' int2str(1024) 'x' int2str(720) ext]);

writerObj = VideoWriter(strcat(mainPath,'_MATRIZ-v',int2str(v),'-','8bit','-', int2str(1024),'x',int2str(720)));
writerObj.FrameRate = 30;
open(writerObj);

nFrames = videoin.NumberOfFrames;

for i=1:nFrames
    
   frame = read(videoin,i);
   
   frame = rgb2gray(frame);
   
   writeVideo(writerObj,frame);    
    
end

close(writerObj);

end

