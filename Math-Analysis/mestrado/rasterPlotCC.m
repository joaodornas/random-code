function rasterPlotCC

video1 = mmreader('video_out-v1-r34.avi');
video2 = mmreader('video_out-v2-r34.avi');
video3 = mmreader('video_out-v3-r34.avi');
video4 = mmreader('video_out-v4-r34.avi');

video5 = mmreader('video_out-v5-r22.avi');
video6 = mmreader('video_out-v6-r22.avi');
% video7 = mmreader('video_out-v7-r22.avi');
video8 = mmreader('video_out-v8-r22.avi');
video9 = mmreader('video_out-v9-r22.avi');
video10 = mmreader('video_out-v10-r22.avi');
video11 = mmreader('video_out-v11-r22.avi');

video1x4 = mmreader('video_out-v1-r136.avi');
video2x4 = mmreader('video_out-v2-r136.avi');
video3x4 = mmreader('video_out-v3-r136.avi');
video4x4 = mmreader('video_out-v4-r136.avi');

video5x4 = mmreader('video_out-v5-r88.avi');
video6x4 = mmreader('video_out-v6-r88.avi');
% video7x4 = mmreader('video_out-v7-r88.avi');
video8x4 = mmreader('video_out-v8-r88.avi');
video9x4 = mmreader('video_out-v9-r88.avi');
video10x4 = mmreader('video_out-v10-r88.avi');
video11x4 = mmreader('video_out-v11-r88.avi');

video5x8 = mmreader('video_out-v5-r176.avi');
video6x8 = mmreader('video_out-v6-r176.avi');
% video7x8 = mmreader('video_out-v7-r176.avi');
video8x8 = mmreader('video_out-v8-r176.avi');
video9x8 = mmreader('video_out-v9-r176.avi');
video10x8 = mmreader('video_out-v10-r176.avi');
video11x8 = mmreader('video_out-v11-r176.avi');


    rasterPlot(video1,1,34)
    rasterPlot(video2,2,34)
    rasterPlot(video3,3,34)
    rasterPlot(video4,4,34)
    
    rasterPlot(video5,5,22)
    rasterPlot(video6,6,22)
    %rasterPlot(video7,7,22)
    rasterPlot(video8,8,22)
    rasterPlot(video9,9,22)
    rasterPlot(video10,10,22)
    rasterPlot(video11,11,22)
    
    rasterPlot(video1x4,1,136)
    rasterPlot(video2x4,2,136)
    rasterPlot(video3x4,3,136)
    rasterPlot(video4x4,4,136)
    
    rasterPlot(video5x4,5,88)
    rasterPlot(video6x4,6,88)
    %rasterPlot(video7x4,7,88)
    rasterPlot(video8x4,8,88)
    rasterPlot(video9x4,9,88)
    rasterPlot(video10x4,10,88)
    rasterPlot(video11x4,11,88)
    
    rasterPlot(video5x8,5,176)
    rasterPlot(video6x8,6,176)
    %rasterPlot(video7x8,7,176)
    rasterPlot(video8x8,8,176)
    rasterPlot(video9x8,9,176)
    rasterPlot(video10x8,10,176)
    rasterPlot(video11x8,11,176)
    
    
    
    
    
    
    function rasterPlot(video,vn,radius);

        nFrames = video.NumberOfFrames;
        
        vidFrames = read(video);
        
        for F=1:nFrames
            
            movframe(F).cdata = vidFrames(:,:,:,F);
        
            imwrite(movframe(F).cdata,strcat('v',int2str(vn),'-radius-',int2str(radius),'-frame-',int2str(F),'.jpg'));
        
        end
        
    end


end