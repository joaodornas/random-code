function dornasCorrelFFT(name,frames,kind,mainPath)

tic

frames = double(frames);

F = fftn(frames);
ZF = conj(F);

FFTR = F.*ZF;

R = ifftn(FFTR);

RF = fftn(R);

data = struct('name',name,'kind',kind,'R',R,'RF',RF);

save(strcat(mainPath,name,'-dornasFourier-analysis-',kind),'data');

clear F;
clear ZF;
clear FFTR;
clear R;
clear RF;
clear data;
clear frames;

toc

end

