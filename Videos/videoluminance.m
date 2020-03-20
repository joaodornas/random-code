function h = videoluminance(video_name)

nFrames = video_name.NumberOfFrames;
vidHeight = video_name.Height;
vidWidth = video_name.Width;

% Preallocate movie structure.
movrgb(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
           'colormap', []);

movgray0(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth),...
           'colormap', gray);
       
movgray64(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth),...
           'colormap', gray);

movgray128(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth),...
           'colormap', gray);

movgray192(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth),...
           'colormap', gray);

       
% Read one frame at a time.
for k = 1 : nFrames
    
    movrgb(k).cdata = read(video_name, k);

    movgray0(k).cdata = rgb2gray(movrgb(k).cdata);
    
    movgray64(k).cdata = movgray0(k).cdata + 64;
    movgray128(k).cdata = movgray0(k).cdata + 128;
    movgray192(k).cdata = movgray0(k).cdata + 192;
    
    movgray64(k).cdata(movgray64(k).cdata>255) = movgray64(k).cdata(movgray64(k).cdata>255) - 256; 
    movgray128(k).cdata(movgray128(k).cdata>255) = movgray128(k).cdata(movgray128(k).cdata>255) - 256;
    movgray192(k).cdata(movgray192(k).cdata>255) = movgray192(k).cdata(movgray192(k).cdata>255) - 256;

end

writerobj0 = VideoWriter('8bit-720x480-0.mp4');
writerobj0.FrameRate = video_name.FrameRate;
open(writerobj0);

writerobj64 = VideoWriter('8bit-720x480-64.mp4');
writerobj64.FrameRate = video_name.FrameRate;
open(writerobj64);

writerobj128 = VideoWriter('8bit-720x480-128.mp4');
writerobj128.FrameRate = video_name.FrameRate;
open(writerobj128);

writerobj192 = VideoWriter('8bit-720x480-192.mp4');
writerobj192.FrameRate = video_name.FrameRate;
open(writerobj192);

for k = 1 : nFrames
    
    writeVideo(writerobj0, movgray0(k).cdata);
    writeVideo(writerobj64, movgray64(k).cdata);
    writeVideo(writerobj128, movgray128(k).cdata);
    writeVideo(writerobj192, movgray192(k).cdata);
    
end

close(writerobj0);
close(writerobj64);
close(writerobj128);
close(writerobj192);

clear movrgb;       
clear movgray0;
clear movgray64;
clear movgray128;
clear movgray192;

end





