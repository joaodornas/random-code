function imgC = imgconv(imgpath,kernel,option)


img = imread(imgpath);

imgR = img(:,:,1);
imgG = img(:,:,2);
imgB = img(:,:,3);

if option == 1
    
    imgCR = conv2(double(kernel),double(imgR));
    imgCG = conv2(double(kernel),double(imgG));
    imgCB = conv2(double(kernel),double(imgB));

else if option == 2
        
    imgCR = filtro(double(imgR),double(kernel));
    imgCG = filtro(double(imgG),double(kernel));
    imgCB = filtro(double(imgB),double(kernel));
       
end
    
imgCR = uint8(imgCR);
imgCG = uint8(imgCG);
imgCB = uint8(imgCB);

imgC = cat(3,imgCR,imgCG,imgCB);

imshow(imgC)


end