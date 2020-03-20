function video2reverseDirect(video_name,nFrames)

videoin = mmreader(video_name);
%videoin = video_name;

%nFrames = videoin.NumberOfFrames;
vidHeight = videoin.Height;
vidWidth = videoin.Width;
FrameRate = videoin.FrameRate;

% Preallocate movie structure.
mov(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
           'colormap', []);

% Read one frame at a time.
for k = 1 : nFrames
    mov(k).cdata = read(videoin, k);
end

mov_reverse = fliplr(mov);

videoout = VideoWriter('filme_direto.mp4');
videoout.FrameRate = FrameRate;
open(videoout);

for k = 1 : nFrames
   
    writeVideo(videoout, mov(k).cdata);
    
end

close(videoout);

videoout = VideoWriter('filme_reverso.mp4');
videoout.FrameRate = FrameRate;
open(videoout);

for k = 1 : nFrames
   
    writeVideo(videoout, mov_reverse(k).cdata);
    
end

close(videoout);

end