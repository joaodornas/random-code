function dongCorrelFFTmovie(video)

tic

name = video(1:end-4);

videoin = mmreader(video);

nFrames = videoin.NumberOfFrames;

width = videoin.Width;

height = videoin.Height;

movie = zeros(nFrames,height,width);

matlabpool open;

parfor i=1:nFrames
    
    frame3D = read(videoin,i);
    
    frame = frame3D(:,:,1);
    
    movie(i,:,:) = frame(:,:);
    
end

frames = double(movie);

F = fftn(frames);

ZF = conj(F);

FFTR = F.*ZF;

R = ifftn(FFTR);

RF = fftn(R);

data = struct('name',name,'RF',RF);

save(strcat(name,'-dongFourier-analysis'),'data','-v7.3');

matlabpool close;

toc

end

