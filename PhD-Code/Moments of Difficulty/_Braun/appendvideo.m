mainmovie = VideoReader('entropies.avi');

appendmovie = VideoReader('appendmovie.avi');

newvideo = VideoWriter('newentropies.avi');
open(newvideo);

nFrames = get(mainmovie, 'NumberOfFrames');

for i=1:nFrames
    
    frame = read(mainmovie,i);
    writeVideo(newvideo,frame);
    
end

newframe = read(appendmovie,1);

writeVideo(mainmovie,newframe);

close(mainmovie);

close all,