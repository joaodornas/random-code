function dornasCorrelFFTWOtime(name,frames,kind,mainPath,with_corr)

tic

frames = double(frames);

if with_corr == 1

    F = fftn(frames);
    ZF = conj(F);

    FFTR = F.*ZF;

    R = ifftn(FFTR);
    
    nTime = size(R,3);
    
    RF = zeros(size(R,1),size(R,2),size(R,3));
    
    for n=1:nTime
        
        frame(:,:) = R(:,:,n);
        
        RF_frame = fft2(frame);
        
        RF(:,:,n) = RF_frame;
        
    end
    
elseif with_corr == 0
    
    nFrames = size(frames,3);
    
    for n=1:nFrames
        
       frame(:,:) = frames(:,:,n); 
        
       RF_frame = fft2(frame);
       
       RF(:,:,n) = RF_frame;
              
    end
    
end


if with_corr == 1
    
    data = struct('name',name,'kind',kind,'R',R,'RF',RF);
    
    corr = 'wi-corr';
    
elseif with_corr == 0
    
    data = struct('name',name,'kind',kind,'RF',RF);
    
    corr = 'wo-corr';
    
end

save(strcat(mainPath,name,'-dornasFourier-analysis-WO-time-',corr,'-',kind),'data');

toc

end

