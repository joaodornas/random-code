   function spatialContrast(video,radius,nFrames)

        readerobj = VideoReader(video);
   
        vidFrames = read(readerobj);
        
        %nFrames = video.NumberOfFrames;
               
        Michelson = zeros(1,nFrames);
        
        for k=1:nFrames
           
            frame = vidFrames(:,:,:,k);
            
            maxLum = double(max(max(max(frame))));
            minLum = double(min(min(min(frame))));
            
            Michelson(k) = (maxLum - minLum) / ( (maxLum + minLum)/2 ) ;
            
        end
        
        meanMichelson = mean(Michelson);
        
        spatialContrast.Michelson = Michelson;
        spatialContrast.meanMichelson = meanMichelson;
        
        save(strcat(video,'-spatialContrast-',int2str(radius)),'spatialContrast');
        
        f = figure;
        
        plot(1:nFrames,Michelson);
        
        video = video(1:end-4);
        
        print(f,'-depsc',strcat(video,'-spatialContrast-',int2str(radius),'-plot')); 
        
        
    end