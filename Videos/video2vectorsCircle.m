function h = video2vectorsCircle(video_name,radius,g)

nFrames = video_name.NumberOfFrames;
vidHeight = video_name.Height;
vidWidth = video_name.Width;

% Preallocate movie structure.

movrgb(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
           'colormap', []);
       
movgray(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth),...
           'colormap', gray);
       
movgraycrop(1:nFrames) = ...
    struct('cdata', zeros(2*radius, 2*radius),...
           'colormap', gray);
       
mask = imcircle(2*radius);

hStart = (vidHeight/2) - radius + 1;
hEnd = (vidHeight/2) + radius;
wStart = (vidWidth/2) - radius + 1;
wEnd = (vidWidth/2) + radius;

h = fspecial('gaussian',g,g);


%Read one frame at a time.

for k = 1 : nFrames
    
    movrgb(k).cdata = read(video_name, k);

    %imwrite(movrgb(k).cdata,['movrgb' int2str(k) '.bmp']);

    movgray(k).cdata = rgb2gray(movrgb(k).cdata);
    
    %imwrite(movgray(k).cdata,['movgray-' int2str(k) '.bmp']);
    
    movgraycrop(k).cdata = movgray(k).cdata(hStart:hEnd,wStart:wEnd);
    
    %imwrite(movgraycrop(k).cdata,['movgraycrop' int2str(k) '.bmp']);
    
    movgraycrop(k).cdata(~mask) = 0;

    %imwrite(movgraycrop(k).cdata,['movgraycropCircle-' int2str(k) '.bmp']);
    
    movgraycrop(k).cdata = imfilter(movgraycrop(k).cdata,h);
    
    %imwrite(movgraycrop(k).cdata,['movgraycropCircleFilter-' int2str(k) '-' int2str(radius) '.bmp']);
    
end

%movie2avi(movgraycrop,['moviecrop-' int2str(radius) '.avi'],'compression','None');

clear movrgb;       
clear movgray;


% Read one frame at a time.

k = 1;
f = 1;

while (k < 1001)
    
    %X = [1];

    %for i = 1 : 480
        %for j = 1: 720
            %X = [X movgraycrop(k).cdata(i,j)];
        %end
    %end

    %X(1) = [];

    X = reshape(movgraycrop(f).cdata.',[],1);
    
    X = X.';
    
    for i = 1 : 3
        if k < 10
            frame = strcat('000', int2str(k));
        elseif k < 100
            frame = strcat('00', int2str(k));
        elseif k < 1000
            frame = strcat('0', int2str(k));
        else
            frame = int2str(k);
        end
    
        dlmwrite(strcat('marinho', '_', int2str(radius*2), '_', int2str(radius*2), '_', frame, '.txt'), X, 'delimiter', '\t', '-append', 'newline', 'pc');
        
        k = k + 1;
        
        i = i + 1;
    end
    
    f = f + 1;

end

clear movgraycrop;


end





