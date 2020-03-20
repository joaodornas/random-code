function h = concatenateVideos(FrameRate,height,width)


video1 = mmreader('videoclip - part 1.avi');
video2 = mmreader('videoclip - part 2.avi');
video3 = mmreader('videoclip - part 3.avi');
video4 = mmreader('videoclip - part 4.avi');
video5 = mmreader('videoclip - part 5.avi');
video6 = mmreader('videoclip - part 6.avi');
video7 = mmreader('videoclip - part 7.avi');
%video8 = mmreader('videoclip - part 8.avi');

nFrames1 = video1.NumberOfFrames;
nFrames2 = video2.NumberOfFrames;
nFrames3 = video3.NumberOfFrames;
nFrames4 = video4.NumberOfFrames;
nFrames5 = video5.NumberOfFrames;
nFrames6 = video6.NumberOfFrames;
nFrames7 = video7.NumberOfFrames;
%nFrames8 = video8.NumberOfFrames;

totalFrames = nFrames1 + nFrames2 + nFrames3 + nFrames4 + nFrames5 + nFrames6 + nFrames7;

movieFrames(1:totalFrames) = struct('cdata',zeros(height,width));

j= 0;

for k=1:nFrames1
   
    j = j + 1;
    
    frame = read(video1,k);
        
    movieFrames(j).cdata = frame;
    
end

clear video1;

for k=1:nFrames2
   
    j = j + 1;
    
    frame = read(video2,k);
        
    movieFrames(j).cdata = frame;
    
end

clear video2;

for k=1:nFrames3
   
    j = j + 1;
    
    frame = read(video3,k);
        
    movieFrames(j).cdata = frame;
    
end

clear video3;

for k=1:nFrames4
   
    j = j + 1;
    
    frame = read(video4,k);
        
    movieFrames(j).cdata = frame;
    
end

clear video4;

for k=1:nFrames5
   
    j = j + 1;
    
    frame = read(video5,k);
        
    movieFrames(j).cdata = frame;
    
end

clear video5;

for k=1:nFrames6
   
    j = j + 1;
    
    frame = read(video6,k);
        
    movieFrames(j).cdata = frame;
    
end

clear video6;

for k=1:nFrames7
   
    j = j + 1;
    
    frame = read(video7,k);
        
    movieFrames(j).cdata = frame;
    
end

clear video7;

% for k=1:nFrames8
%    
%     j = j + 1;
%     
%     frame = read(video8,k);
%         
%     movieFrames(j).cdata = frame;
%     
% end
% 
% clear video8;



writerobj = VideoWriter('concatenatedVideo');
writerobj.FrameRate = FrameRate;
open(writerobj);

for i=1:totalFrames

    writeVideo(writerobj, movieFrames(i).cdata);
    
end

close(writerobj);


end

