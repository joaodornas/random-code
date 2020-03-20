function visualFeatures


for i=3:3

    mainPath = strcat('/Users/joaodornas/Documents/_Research/Cenas Naturais/VIDEOS/VIDEO-',int2str(i),'/');

    visualFolder = 'visualFeature';

    mkdir(mainPath,visualFolder);

    visualFeaturePath = strcat('/Users/joaodornas/Documents/_Research/Cenas Naturais/VIDEOS/VIDEO-',int2str(i),'/',visualFolder,'/');

    videoObject1 = mmreader(strcat('video_out-v',int2str(i),'-r22.avi'));
    videoObject2 = mmreader(strcat('video_out-v',int2str(i),'-r33.avi'));
    videoObject3 = mmreader(strcat('video_out-v',int2str(i),'-r34.avi'));
    videoObject4 = mmreader(strcat('video_out-v',int2str(i),'-r44.avi'));
    videoObject5 = mmreader(strcat('video_out-v',int2str(i),'-r66.avi'));
    videoObject6 = mmreader(strcat('video_out-v',int2str(i),'-r68.avi'));
    videoObject7 = mmreader(strcat('video_out-v',int2str(i),'-r88.avi'));
    videoObject8 = mmreader(strcat('video_out-v',int2str(i),'-r99.avi'));
    videoObject9 = mmreader(strcat('video_out-v',int2str(i),'-r102.avi'));
    videoObject10 = mmreader(strcat('video_out-v',int2str(i),'-r110.avi'));
    videoObject11 = mmreader(strcat('video_out-v',int2str(i),'-r132.avi'));
    videoObject12 = mmreader(strcat('video_out-v',int2str(i),'-r136.avi'));
    videoObject13 = mmreader(strcat('video_out-v',int2str(i),'-r165.avi'));
    videoObject14 = mmreader(strcat('video_out-v',int2str(i),'-r170.avi'));
    videoObject15 = mmreader(strcat('video_out-v',int2str(i),'-r176.avi'));
    videoObject16 = mmreader(strcat('video_out-v',int2str(i),'-r198.avi'));
    videoObject17 = mmreader(strcat('video_out-v',int2str(i),'-r204.avi'));
    videoObject18 = mmreader(strcat('video_out-v',int2str(i),'-r220.avi'));
    videoObject19 = mmreader(strcat('video_out-v',int2str(i),'-r264.avi'));
    videoObject20 = mmreader(strcat('video_out-v',int2str(i),'-r272.avi'));
    videoObject21 = mmreader(strcat('video_out-v',int2str(i),'-r348.avi'));
    
    spatialContrast(videoObject1,visualFeaturePath,i,22);
    timeContrast(videoObject1,visualFeaturePath,i,22);
    orientation(videoObject1,visualFeaturePath,i,22);

    spatialContrast(videoObject2,visualFeaturePath,i,33);
    timeContrast(videoObject2,visualFeaturePath,i,33);
    orientation(videoObject2,visualFeaturePath,i,33);

    spatialContrast(videoObject3,visualFeaturePath,i,34);
    timeContrast(videoObject3,visualFeaturePath,i,34);
    orientation(videoObject3,visualFeaturePath,i,34);

    spatialContrast(videoObject4,visualFeaturePath,i,44);
    timeContrast(videoObject4,visualFeaturePath,i,44);
    orientation(videoObject4,visualFeaturePath,i,44);

    spatialContrast(videoObject5,visualFeaturePath,i,66);
    timeContrast(videoObject5,visualFeaturePath,i,66);
    orientation(videoObject5,visualFeaturePath,i,66);

    spatialContrast(videoObject6,visualFeaturePath,i,68);
    timeContrast(videoObject6,visualFeaturePath,i,68);
    orientation(videoObject6,visualFeaturePath,i,68);

    spatialContrast(videoObject7,visualFeaturePath,i,88);
    timeContrast(videoObject7,visualFeaturePath,i,88);
    orientation(videoObject7,visualFeaturePath,i,88);

    spatialContrast(videoObject8,visualFeaturePath,i,99);
    timeContrast(videoObject8,visualFeaturePath,i,99);
    orientation(videoObject8,visualFeaturePath,i,99);

    spatialContrast(videoObject9,visualFeaturePath,i,102);
    timeContrast(videoObject9,visualFeaturePath,i,102);
    orientation(videoObject9,visualFeaturePath,i,102);

    spatialContrast(videoObject10,visualFeaturePath,i,110);
    timeContrast(videoObject10,visualFeaturePath,i,110);
    orientation(videoObject10,visualFeaturePath,i,110);

    spatialContrast(videoObject11,visualFeaturePath,i,132);
    timeContrast(videoObject11,visualFeaturePath,i,132);
    orientation(videoObject11,visualFeaturePath,i,132);

    spatialContrast(videoObject12,visualFeaturePath,i,136);
    timeContrast(videoObject12,visualFeaturePath,i,136);
    orientation(videoObject12,visualFeaturePath,i,136);

    spatialContrast(videoObject13,visualFeaturePath,i,165);
    timeContrast(videoObject13,visualFeaturePath,i,165);
    orientation(videoObject13,visualFeaturePath,i,165);

    spatialContrast(videoObject14,visualFeaturePath,i,170);
    timeContrast(videoObject14,visualFeaturePath,i,170);
    orientation(videoObject14,visualFeaturePath,i,170);

    spatialContrast(videoObject15,visualFeaturePath,i,176);
    timeContrast(videoObject15,visualFeaturePath,i,176);
    orientation(videoObject15,visualFeaturePath,i,176);

    spatialContrast(videoObject16,visualFeaturePath,i,198);
    timeContrast(videoObject16,visualFeaturePath,i,198);
    orientation(videoObject16,visualFeaturePath,i,198);

    spatialContrast(videoObject17,visualFeaturePath,i,204);
    timeContrast(videoObject17,visualFeaturePath,i,204);
    orientation(videoObject17,visualFeaturePath,i,204);

    spatialContrast(videoObject18,visualFeaturePath,i,220);
    timeContrast(videoObject18,visualFeaturePath,i,220);
    orientation(videoObject18,visualFeaturePath,i,220);

    spatialContrast(videoObject19,visualFeaturePath,i,264);
    timeContrast(videoObject19,visualFeaturePath,i,264);
    orientation(videoObject19,visualFeaturePath,i,264);

    spatialContrast(videoObject20,visualFeaturePath,i,272);
    timeContrast(videoObject20,visualFeaturePath,i,272);
    orientation(videoObject20,visualFeaturePath,i,272);

    spatialContrast(videoObject21,visualFeaturePath,i,348);
    timeContrast(videoObject21,visualFeaturePath,i,348);
    orientation(videoObject21,visualFeaturePath,i,348);
    
    close all;

end


    function spatialContrast(video,path,v,radius)

        vidFrames = read(video);
        
        nFrames = video.NumberOfFrames;
        
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
        
        save(strcat(path,'v',int2str(v),'-radius-',int2str(radius),'-spatialContrast'),'spatialContrast');
        
        f = figure;
        
        plot(1:nFrames,Michelson);
        
        print(f,'-djpeg',strcat(path,'v',int2str(v),'-radius-',int2str(radius),'-spatialContrast-plot')); 
        
        
    end


    function timeContrast(video,path,v,radius)
        
       
        vidFrames = read(video);
        
        nFrames = video.NumberOfFrames;
        
        for k=1:nFrames
           
            frame = vidFrames(:,:,:,k);
            
            if k ~= nFrames
            
                nextFrame = vidFrames(:,:,:,k+1);
                
                S = abs(nextFrame - frame);
                
                meanS(k) = mean(mean(mean(S)));
            
            end
            
        end
        
        timeContrast.pixelWiseChange = meanS;
        
        save(strcat(path,'v',int2str(v),'-radius-',int2str(radius),'-timeContrast'),'timeContrast');
        
        f = figure;
        
        plot(1:nFrames-1,meanS);
        
        print(f,'-djpeg',strcat(path,'v',int2str(v),'-radius-',int2str(radius),'-timeContrast-plot')); 
        
    end

    function orientation(video,path,v,radius)
        
       vidFrames = read(video);
       
       nFrames = video.NumberOfFrames;
       
       for k=1:nFrames
          
           frame = vidFrames(:,:,:,k);
           
           grayframe = rgb2gray(frame);
           
           [gradient or] = canny(grayframe,1);
         
           angle(k) = mean(mean(or));        
           
       end
        
       orientation.angle = angle;
       
       orientation.sinAngle = sin(angle);
       
        save(strcat(path,'v',int2str(v),'-radius-',int2str(radius),'-orientation'),'orientation');
        
        f = figure;
        
        plot(1:nFrames,orientation.sinAngle);
        
        print(f,'-djpeg',strcat(path,'v',int2str(v),'-radius-',int2str(radius),'-orientationSin-plot')); 
        
    end

end