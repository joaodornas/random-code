function outmovie = redundancy(nFrames,level,width,height,radius)

tic

for v=1:23
    
    readerobj(v) = VideoReader(['_MATRIZ-v' int2str(v) '-8bit-' int2str(width) 'x' int2str(height) '.avi'], 'tag', 'myreader');

end

w = 4;

for v=1:23
    
    for i=1:nFrames

        frame = read(readerobj(v),i); 

        frame = rgb2gray(frame);

        mask = imcircle(2*radius);

        hStart = (height/2) - radius;
        hEnd = (height/2) + radius - 1;
        wStart = (width/2) - radius;
        wEnd = (width/2) + radius - 1;

        h = fspecial('gaussian',w,w);

        frame = frame(hStart:hEnd,wStart:wEnd);

        frame(frame==0) = 1;

        frame(~mask) = 255;

        frame = vertcat(ones(radius,2*radius)*255,frame);

        frame = vertcat(frame,ones(radius,2*radius)*255);

        frame = horzcat(frame,ones(4*radius,radius)*255);

        frame = horzcat(ones(4*radius,radius)*255,frame);

        frame = imfilter(frame,h);

        length = size(frame,2);

        increase = 4;

        start_ = length/2 - radius - increase;
        end_ = length/2 + radius - 1 + increase;

        frame = frame(start_:end_,start_:end_);

        frames(i,:,:) = frame;

        pixels = size(frames,2);

    end
    
    allframes(v,1:nFrames,1:pixels,1:pixels) = frames(1:nFrames,:,:);

end

disp(strcat('Time:',num2str(toc),'-','Frames'));

combFrames1 = randperm(nFrames);

combFrames2 = combFrames1 + 1;

combFrames2(combFrames2==nFrames+1) = nFrames;

gpuFrames1 = gpuArray(frames(:,:,combFrames1));

gpuFrames2 = gpuArray(frames(:,:,combFrames2));

[getmedian getdots] = arrayfun(@dotproduct,gpuFrames1,gpuFrames2,nFrames,pixels);

while  getmedian > level
    
    combFrames1 = randperm(nFrames);
    
    combFrames2 = combFrames1 + 1;
    
    combFrames2(combFrames2==nFrames+1) = nFrames;
    
    gpuFrames1 = gpuArray(frames(:,:,combFrames1));
    
    gpuFrames2 = gpuArray(frames(:,:,combFrames2));

    [getmedian getdots] = arrayfun(@dotproduct,gpuFrames1,gpuFrames2,nFrames,pixels);
    
    new_Frames = frames(:,:,combFrames1);

    disp(strcat('Time:',num2str(toc),'-','Median:',num2string(getmedian)));
    
end

data = struct('getmedian',getmedian,'getdots',getdots);

save('data','data');

videoout = VideoWriter(['v' int2str(v) '-8bit-' int2str(width) 'x' int2str(height)' '-redundancy-' int2str(getmedian)]);
    
videoout.FrameRate = FrameRate;
    
open(videoout);

for i=1:nFrames
    
    writeVideo(videoout, new_Frames(:,:,i));

end

close(videoout);

function [mediandot dots] = dotproduct(frames1,frames2,nFrames,pixels)

       Total_Dot = double(frames1.*frames2);  

       Norma = sqrt(frames1^2.*frames2^2);

       Total = Total_Dot / Norma ;

       dots = Total;

       mediandot = median(median(dots));
    
end


end

