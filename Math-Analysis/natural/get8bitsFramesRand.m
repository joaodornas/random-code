function frame = get8bitsFramesRand(v)

videoin = VideoReader(['_MATRIZ-v' int2str(v) '-' '8bit' '-' int2str(1024) 'x' int2str(720) '.avi']);

nFrames = videoin.NumberOfFrames;

aleatorio = round(rand(1,nFrames)*nFrames);

frame = read(videoin,aleatorio(1));


end

