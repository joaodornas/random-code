function scramble(picture)

    Im = imread(picture);
    
    ImSize = size(Im);

    for layer = 1:ImSize(3)
        
        ImFourier(:,:,layer) = fft2(Im(:,:,layer));       
    %Fast-Fourier transform
       
        Amp(:,:,layer) = abs(ImFourier(:,:,layer)); 
        %Amp(:,:,layer) = abs(fft2(rand(ImSize(1), ImSize(2))));
    %amplitude spectrum
  
        %Phase(:,:,layer) = angle(ImFourier(:,:,layer));
        Phase(:,:,layer) = angle(fft2(rand(ImSize(1), ImSize(2))));
    %generate random phase structure
        
        ImScrambled(:,:,layer) = ifft2(Amp(:,:,layer).*exp(sqrt(-1)*(Phase(:,:,layer))));   
    %combine Amp and Phase then perform inverse Fourier
    
    end

    ImScrambled = real(ImScrambled); 
    %get rid of imaginery part in image (due to rounding error)
    
    imwrite(ImScrambled,'scrambled.jpg','jpg');
    
   
    
end
