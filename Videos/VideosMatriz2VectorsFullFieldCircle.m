function h = VideosMatriz2VectorsFullFieldCircle

    G = 4;
    
    numFrames = 300;
    
    radius(1).px = 24/2;
    radius(2).px = 32/2;
    radius(3).px = 40/2;
    radius(4).px = 48/2;
    radius(5).px = 56/2;
    
    height = 720;
    width = 1024;
    crop = 100;
    
    FrameRate = 30;
    
    %video_name = mmreader(['_MATRIZ-v' int2str(v) '-8bit-' int2str(width) 'x' int2str(height) '.avi']);
    Videos(1).readerobj = VideoReader(['_MATRIZ-v' int2str(1) '-8bit-' int2str(width) 'x' int2str(height) '.avi'], 'tag', 'myreader');
    Videos(1).v = 1;
%     Videos(2).readerobj = VideoReader(['_MATRIZ-v' int2str(2) '-8bit-' int2str(width) 'x' int2str(height) '.avi'], 'tag', 'myreader');
%     Videos(2).v = 2;
%     Videos(3).readerobj = VideoReader(['_MATRIZ-v' int2str(3) '-8bit-' int2str(width) 'x' int2str(height) '.avi'], 'tag', 'myreader');
%     Videos(3).v = 3;
%     Videos(4).readerobj = VideoReader(['_MATRIZ-v' int2str(4) '-8bit-' int2str(width) 'x' int2str(height) '.avi'], 'tag', 'myreader');
%     Videos(4).v = 4;
%     Videos(5).readerobj = VideoReader(['_MATRIZ-v' int2str(8) '-8bit-' int2str(width) 'x' int2str(height) '.avi'], 'tag', 'myreader');
%     Videos(5).v = 8;
%     Videos(6).readerobj = VideoReader(['_MATRIZ-v' int2str(11) '-8bit-' int2str(width) 'x' int2str(height) '.avi'], 'tag', 'myreader');
%     Videos(6).v = 11;
        
    Togray = 0;
    
    vectorsVideos(Videos,radius,height,width,G,numFrames,Togray,FrameRate,crop);

end

function vectorsVideos(Videos,radius,height,width,g,nFrames,Togray,FrameRate,crop)

                
            nVideos = length(Videos);
            nRadius = length(radius);
            
            for iVideo=1:nVideos
                
                Videos(iVideo).FolderFFFor = strcat('v','-',int2str(Videos(iVideo).v),'-','FullField','-','forward','-',int2str(width-2*crop),'x',int2str(height-2*crop));
                Videos(iVideo).FolderFFBack = strcat('v','-',int2str(Videos(iVideo).v),'-','FullField','-','backward','-',int2str(width-2*crop),'x',int2str(height-2*crop));
            
                mkdir(Videos(iVideo).FolderFFFor);
                mkdir(Videos(iVideo).FolderFFBack); 
            
                for iRadius=1:nRadius
                
                    Videos(iVideo).radius(iRadius).FolderCCFor = strcat('v','-',int2str(Videos(iVideo).v),'-','CenterContour','-','forward','-',int2str(radius(iRadius).px*2),'x',int2str(radius(iRadius).px*2));
                    Videos(iVideo).radius(iRadius).FolderCCBack = strcat('v','-',int2str(Videos(iVideo).v),'-','CenterContour','-','backward','-',int2str(radius(iRadius).px*2),'x',int2str(radius(iRadius).px*2));
                
                    mkdir(Videos(iVideo).radius(iRadius).FolderCCFor);
                    mkdir(Videos(iVideo).radius(iRadius).FolderCCBack);
                
                end
                
            end

            for iVideo=1:nVideos
            
                movgray(1:nFrames) = struct('cdata', zeros(height, width),'colormap', gray);
                
                writerObjFFFor = VideoWriter(strcat('video_out-v-',int2str(Videos(iVideo).v),'-','FullField','-','forward','-',int2str(width-2*crop),'x',int2str(height-2*crop),'.avi'));
                writerObjFFFor.FrameRate = FrameRate;
                open(writerObjFFFor);
                
                writerObjFFBack = VideoWriter(strcat('video_out-v-',int2str(Videos(iVideo).v),'-','FullField','-','backward','-',int2str(width-2*crop),'x',int2str(height-2*crop),'.avi'));
                writerObjFFBack.FrameRate = FrameRate;
                open(writerObjFFBack);
                
                vidFrames = read(Videos(iVideo).readerobj);
                
                for k=1:nFrames

                    frame = vidFrames(:,:,:,k);
                
                    if Togray == 1
                    
                        movgray(k).cdata = rgb2gray(frame);
                        
                    else
                        
                        movgray(k).cdata = frame;
                        
                    end
                    
                    movgray(k).cdata = movgray(k).cdata(crop+1:height-crop,crop+1:width-crop);
                    
                    writeVideo(writerObjFFFor,movgray(k).cdata);
                    
                end
                
                frame50 = vidFrames(crop+1:height-crop,crop+1:width-crop,:,50);
                
                imwrite(movgray(50).cdata,strcat('frame-',int2str(Videos(iVideo).v),'-',int2str(width-2*crop),'-',int2str(height-2*crop),'.jpg'));
                
                saveVectors(movgray,Videos(iVideo).FolderFFFor,nFrames);
                
                movgrayReverse = fliplr(movgray);
                
                for k=1:nFrames

                    writeVideo(writerObjFFBack,movgrayReverse(k).cdata);
                    
                end
                
                saveVectors(movgrayReverse,Videos(iVideo).FolderFFBack,nFrames);
                
                close(writerObjFFFor);
                close(writerObjFFBack);
                
                h = fspecial('gaussian',g,g);
                
                height = height - 2*crop;
                width = width - 2*crop;
                
                for iRadius=1:nRadius

                    movgraycrop(1:nFrames) = struct('cdata', zeros(2*radius(iRadius).px, 2*radius(iRadius).px),'colormap', gray);
                    
                    mask = imcircle(2*radius(iRadius).px);
                    
                    hStart = (height/2) - radius(iRadius).px;
                    hEnd = (height/2) + radius(iRadius).px - 1;
                    wStart = (width/2) - radius(iRadius).px;
                    wEnd = (width/2) + radius(iRadius).px - 1;
                    
                    yellow = uint8([255 255 0]);
                    
                    X = width/2;
                    Y = height/2;
                    
                    yellow = uint8([255 255 0]);

                    %circles = int32([X Y radius(iRadius).px;X Y radius(iRadius).px-1;X Y radius(iRadius).px-2;X Y radius(iRadius).px-3;X Y radius(iRadius).px-4;X Y radius(iRadius).px-5;X Y radius(iRadius).px-6;X Y radius(iRadius).px-7]);
                    circles = int32([X Y radius(iRadius).px;X Y radius(iRadius).px-1;X Y radius(iRadius).px-2;X Y radius(iRadius).px-3]);

                    shapeInserter = vision.ShapeInserter('Shape','Circles','BorderColor','Custom','CustomBorderColor',yellow);

                    frame = step(shapeInserter,frame50,circles);
                    
                    imwrite(frame,strcat('frame-RF-',int2str(Videos(iVideo).v),'-',int2str(radius(iRadius).px*2),'x',int2str(radius(iRadius).px*2),'-',int2str(width/2),'-',int2str(height/2),'.jpg'));
                    
                    writerObjCCFor = VideoWriter(strcat('video_out-v-',int2str(Videos(iVideo).v),'-','CenterContour','-','forward','-',int2str(radius(iRadius).px*2),'x',int2str(radius(iRadius).px*2),'.avi'));
                    writerObjCCFor.FrameRate = FrameRate;
                    open(writerObjCCFor);
                    
                    writerObjCCBack = VideoWriter(strcat('video_out-v-',int2str(Videos(iVideo).v),'-','CenterContour','-','backward','-',int2str(radius(iRadius).px*2),'x',int2str(radius(iRadius).px*2),'.avi'));
                    writerObjCCBack.FrameRate = FrameRate;
                    open(writerObjCCBack);
                    
                    for k=1:nFrames
                        
                        mask_value = 0;
                        
                        %imwrite(movgray(k).cdata,'framegray.jpg');

                        movgraycrop(k).cdata = movgray(k).cdata(hStart:hEnd,wStart:wEnd);

                        movgraycrop(k).cdata(movgraycrop(k).cdata==0) = 1;

                        %imwrite(movgraycrop(k).cdata,strcat('framecrop-',int2str(radius),'.jpg'));

                        movgraycrop(k).cdata(~mask) = mask_value; %%% SEMPRE COLOCAR MASCAR IGUAL A ZERO !!!!             

                        movgraycrop(k).cdata = vertcat(ones(radius(iRadius).px,2*radius(iRadius).px)*mask_value,movgraycrop(k).cdata);

                        movgraycrop(k).cdata = vertcat(movgraycrop(k).cdata,ones(radius(iRadius).px,2*radius(iRadius).px)*mask_value);

                        movgraycrop(k).cdata = horzcat(movgraycrop(k).cdata,ones(4*radius(iRadius).px,radius(iRadius).px)*mask_value);

                        movgraycrop(k).cdata = horzcat(ones(4*radius(iRadius).px,radius(iRadius).px)*mask_value,movgraycrop(k).cdata);

                        movgraycrop(k).cdata = imfilter(movgraycrop(k).cdata,h);

                        lengthcrop = size(movgraycrop(k).cdata,2);

                        increase = 4;

                        start_ = lengthcrop/2 - radius(iRadius).px - increase;
                        end_ = lengthcrop/2 + radius(iRadius).px - 1 + increase;

                        movgraycrop(k).cdata = movgraycrop(k).cdata(start_:end_,start_:end_);

                        writeVideo(writerObjCCFor,movgraycrop(k).cdata);
                        
                    end
                    
                    imwrite(movgraycrop(50).cdata,strcat('framemask-',int2str(Videos(iVideo).v),'-',int2str(radius(iRadius).px*2),'.jpg'));
                    
                    saveVectors(movgraycrop,Videos(iVideo).radius(iRadius).FolderCCFor,nFrames);
                    
                    movgraycropReverse = fliplr(movgraycrop);
                    
                    for k=1:nFrames
                        
                        writeVideo(writerObjCCBack,movgraycropReverse(k).cdata);
                        
                    end
                    
                    saveVectors(movgraycropReverse,Videos(iVideo).radius(iRadius).FolderCCBack,nFrames);
                    
                    close(writerObjCCFor);
                    close(writerObjCCBack);
                    
                end
                
            end

end

function saveVectors(frames,folder,nFrames)

            k = 1;        
            while (k < nFrames + 1)

                if length(size(frames(k).cdata)) == 3
                    
                    frame = frames(k).cdata(:,:,1);
                    
                else
                    
                    frame = frames(k).cdata;
                    
                end
                
                X = reshape(frame.',[],1);

                X = X.';

                if k < 10
                    frameindex = strcat('00', int2str(k));
                elseif k < 100
                    frameindex = strcat('0', int2str(k));
%                 elseif k < 1000
%                     frameindex = strcat('0', int2str(k));
                else
                    frameindex = int2str(k);
                end

                filename = strcat(folder,'_', frameindex);
                
                fullPath = strcat(folder,'/',filename);

                dlmwrite(strcat(fullPath, '.txt'), X, 'delimiter', '\t', '-append', 'newline', 'pc');

                k = k + 1;

            end

    
end
    



