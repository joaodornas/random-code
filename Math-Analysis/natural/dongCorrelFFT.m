function dongCorrelFFT(name,reshapeframes,kind,mainPath)

tic

frames = double(reshapeframes);

F = fft2(frames);
ZF = conj(F);

FFTR = F.*ZF;

R = ifft2(FFTR);

RF = fft2(R);

data = struct('name',name,'kind',kind,'R',R,'RF',RF);

save(strcat(mainPath,name,'-dongFourier-analysis-',kind),'data');

clear F;
clear ZF;
clear FFTR;
clear R;
clear RF;
clear data;
clear reshapeframes;

toc

end

