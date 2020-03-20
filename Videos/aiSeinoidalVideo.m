
function aiSeinoidalVideo(imagePath,nFrames,T)

img = imread(imagePath);

imggray = rgb2gray(img);

height = size(imggray,1);
width = size(imggray,2);

             writerObj = VideoWriter(strcat('v_ai-','T-',int2str(T),'.avi'));
             writerObj.FrameRate = 30;
             open(writerObj);

s = 0;             
for k=1:nFrames
    
    s = s + T;
    
    for i=0:255

       [row, column] = find(imggray == i) ;

       pixels = length(row);

       for j=1:pixels
           
           movie(row(j),column(j)) = abs(imggray(row(j),column(j)).*sin((pi/2)*(s/nFrames)*i)) ;

       end

    end
    
    %imwrite(movie,strcat('frame-',int2str(k),'.jpg'));
    
    writeVideo(writerObj,movie);
    
end

close(writerObj);

end

