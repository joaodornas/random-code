function randFourier(picture)

    Im = imread(picture);
    
    if size(Im,3) == 3
        
        grayIm = rgb2gray(Im);
        
        Im = grayIm;
        
    end
  
    ImFourier = fft2(Im);       

    Amp = abs(ImFourier); 
  
    Phase = angle(ImFourier);
    
    w = size(Im,1);
    h = size(Im,2);
  
    permAmpW = randperm(w);
    permAmpH = randperm(h);
    
    randAmp = Amp(permAmpW,1:h) ;
    randAmp = randAmp(1:w,permAmpH);
    
    permPhaseW = randperm(w);
    permPhaseH = randperm(h);
    
    randPhase = Phase(permPhaseW,1:h) ;
    randPhase = randPhase(1:w,permPhaseH);
        
    ImScrambledAmp = ifft2(randAmp.*exp(sqrt(-1)*(Phase)));   
    
    ImScrambledPhase = ifft2(Amp.*exp(sqrt(-1)*(randPhase))); 
   
    ImScrambledAmpreal = real(ImScrambledAmp); 
    
    ImScrambledPhasereal = real(ImScrambledPhase); 
    
    imwrite(ImScrambledAmpreal,strcat(picture(1:end-4),'-ImScrambledAmpreal.jpg'),'jpg');
    
    imwrite(ImScrambledPhasereal,strcat(picture(1:end-4),'-ImScrambledPhasereal.jpg'),'jpg');
    
end