  function timeContrast(video,radius,nFrames)
           
        readerobj = VideoReader(video);
  
        vidFrames = read(readerobj);
        
        %nFrames = video.NumberOfFrames;
        
        for k=1:nFrames
           
            frame = vidFrames(:,:,:,k);
            
            pixel(k) = frame(128,128);
            
            if k ~= nFrames
            
                nextFrame = vidFrames(:,:,:,k+1);
                
                S = abs(nextFrame - frame);
                
                meanS(k) = mean(mean(mean(S)));
            
            end
            
        end
        
        timeContrast.pixelWideChange = meanS;
        timeContrast.pixelCenterChange = pixel;
        
        video = video(1:end-4);
        
        save(strcat(video,'-timeContrast-',int2str(radius)),'timeContrast');
        
        f = figure;
        
        plot(1:nFrames-1,meanS);
        
        hold on;
        
        g = figure; 
        
        plot(1:nFrames,pixel);
        
        print(f,'-depsc',strcat(video,'-timeContrast-pixelWideChange-',int2str(radius),'-plot')); 

        print(g,'-depsc',strcat(video,'-timeContrast-pixelCenterChange-',int2str(radius),'-plot')); 

    end
