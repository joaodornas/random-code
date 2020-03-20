   function orientation(video,radius,nFrames)
        
       readerobj = VideoReader(video);
   
       vidFrames = read(readerobj);
       
       %nFrames = video.NumberOfFrames;
       
       for k=1:nFrames
          
           frame = vidFrames(:,:,:,k);
           
           grayframe = rgb2gray(frame);
           
           [gradient or] = canny(grayframe,1);
         
           angle(k) = mean(mean(or));        
           
       end
        
       orientation.angle = angle;
       
       orientation.sinAngle = sin(angle);
       
        save(strcat(video,'-orientation-',int2str(radius)),'orientation');
        
        f = figure;
        
        plot(1:nFrames,orientation.sinAngle);
        
        video = video(1:end-1);
        
        print(f,'-depsc',strcat(video,'-orientation-',int2str(radius),'-plot')); 
        
    end