function stdImg = stdImg(name,f,y,x)


frame(1:f) = struct('cdata', zeros(y, x, 3, 'uint8'), 'colormap', []);
framegray(1:f) = struct('cdata', zeros(y, x), 'colormap', gray);

for i=1:f
 
    fullname = strcat(name, int2str(i), '.png');
    frame(i).cdata = imread(fullname);
    framegray(i).cdata = rgb2gray(frame(i).cdata); 
    
    %imwrite(framegray(i).cdata,['gray-' int2str(i) '.png']);

end

clear frame;

imgM = framegray(1).cdata;

for i=2:f
   
    imgM(:,:,i) = framegray(i).cdata;
    
end

clear framegray;

imgM = double(imgM);
stdM = std(imgM,0,3);

clear imgM;

stdImg = stdM;

end