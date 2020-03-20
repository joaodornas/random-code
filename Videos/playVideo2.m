function h = playVideo2(video_name)

% xyloObj = mmreader('xylophone.mpg');

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

% Size a figure based on the video's width and height.
hf = figure;
set(hf, 'position', [150 150 vidWidth vidHeight])

% Play back the movie once at the video's frame rate.
h = movie(hf, mov, 1, video_name.FrameRate);
h;