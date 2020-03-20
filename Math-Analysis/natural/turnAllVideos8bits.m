function turnAllVideos8bits


for v=13:23
    
    disp(strcat('video:',int2str(v)));
   
    if (v >= 5) && (v <= 11)
        
        ext = '.avi';
        
    elseif (v >= 12) && (v <= 23)
        
        ext = '.mp4';
        
    end
    
    mainPath = strcat('/Volumes/Data/Cenas Naturais/VIDEOS/VIDEO-',int2str(v),'/');
    
    folderPath = strcat('/Volumes/Data/Cenas Naturais/PICTURES/Random/');
    
    
    get8bitsVideos(v,ext,mainPath);
    
    frame = get8bitsFramesRand(v);
    
    imwrite(frame,strcat(folderPath,'random-frame-video-',int2str(v),'.jpg'));
    
    
end


end

