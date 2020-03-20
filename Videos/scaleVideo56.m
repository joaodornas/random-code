function scaleVideo56

 readerobj1chimp = VideoReader(['_MATRIZ-CHIMP-v' int2str(5) '-' int2str(720) 'x' int2str(480) '.mp4'], 'tag', 'myreader');
 readerobj2chimp = VideoReader(['_MATRIZ-CHIMP-v' int2str(6) '-' int2str(720) 'x' int2str(480) '.mp4'], 'tag', 'myreader');

 vidFrames1 = read(readerobj1chimp);
 vidFrames2 = read(readerobj2chimp);
 
 numFrames1 = get(readerobj1chimp, 'NumberOfFrames');
 numFrames2 = get(readerobj2chimp, 'NumberOfFrames');
 
%  writerObj1 = VideoWriter(strcat('_MATRIZ_CHIMP-v',int2str(5),'-',int2str(1024),'x',int2str(720),'.avi'));
%  writerObj1.FrameRate = 30;
%  open(writerObj1);
 
 writerObj2 = VideoWriter(strcat('_MATRIZ_CHIMP-v',int2str(6),'-',int2str(1024),'x',int2str(720),'.avi'));
 writerObj2.FrameRate = 30;
 open(writerObj2);
 
%  for k = 1 : numFrames1
% 
%         frame1 = vidFrames1(:,:,:,k);
%         
%         newFrame1 = imresize(frame1,[720 1024]);
%         
%         writeVideo(writerObj1,newFrame1);
%         
%  end
%  
%   close(writerObj1);
  
  
 for k = 1 : numFrames2
       
     frame2 = vidFrames2(:,:,:,k);

     newFrame2 = imresize(frame2,[720 1024]);

     writeVideo(writerObj2,newFrame2);
 
 end

 close(writerObj2);
 
end

