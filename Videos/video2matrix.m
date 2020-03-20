function h = video2matrix(video_name)

nFrames = video_name.NumberOfFrames;
vidHeight = video_name.Height;
vidWidth = video_name.Width;

% Preallocate movie structure.
mov(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
           'colormap', []);

% Read one frame at a time.
for k = 1 : nFrames
    mov(k).cdata = read(video_name, k);
end
h = mov;
h;