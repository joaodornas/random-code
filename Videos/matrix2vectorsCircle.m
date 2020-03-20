function h = matrix2vectorsCircle(numFrames)

    G = 4;
    
    'radius = 34/2'; '44/2';
    'radius = 68/2';
    'radius = 136/2'; '176/2';
    'radius = 272/2';
    
    %video_name = mmreader(['_MATRIZ-v' int2str(v) '-8bit-' int2str(width) 'x' int2str(height) '.avi']);
    readerobj1 = VideoReader(['_MATRIZ-v' int2str(1) '-8bit-' int2str(1024) 'x' int2str(720) '.avi'], 'tag', 'myreader');
    readerobj2 = VideoReader(['_MATRIZ-v' int2str(2) '-8bit-' int2str(1024) 'x' int2str(720) '.avi'], 'tag', 'myreader');
    readerobj3 = VideoReader(['_MATRIZ-v' int2str(3) '-8bit-' int2str(1024) 'x' int2str(720) '.avi'], 'tag', 'myreader');
    readerobj4 = VideoReader(['_MATRIZ-v' int2str(4) '-8bit-' int2str(1024) 'x' int2str(720) '.avi'], 'tag', 'myreader');
    readerobj5 = VideoReader(['_MATRIZ-CHIMP-v' int2str(5) '-' int2str(1024) 'x' int2str(720) '.avi'], 'tag', 'myreader');
    readerobj6 = VideoReader(['_MATRIZ-CHIMP-v' int2str(6) '-' int2str(1024) 'x' int2str(720) '.avi'], 'tag', 'myreader');
    readerobj7 = VideoReader(['_MATRIZ-v' int2str(7) '-' int2str(1024) 'x' int2str(720) '.avi'], 'tag', 'myreader');
    readerobj8 = VideoReader(['_MATRIZ-v' int2str(8) '-' int2str(1024) 'x' int2str(720) '.avi'], 'tag', 'myreader');
    readerobj9 = VideoReader(['_MATRIZ-v' int2str(9) '-' int2str(1024) 'x' int2str(720) '.avi'], 'tag', 'myreader');
    readerobj10 = VideoReader(['_MATRIZ-v' int2str(10) '-' int2str(1024) 'x' int2str(720) '.avi'], 'tag', 'myreader');
    readerobj11 = VideoReader(['_MATRIZ-v' int2str(11) '-' int2str(1024) 'x' int2str(720) '.avi'], 'tag', 'myreader');
    
    %nFrames = video_name.NumberOfFrames;
    vidHeight = readerobj1.Height;
    vidWidth = readerobj1.Width;
    
%     vectorsVideos(readerobj1,78,vidHeight,vidWidth,G,numFrames,0,1,2);
%       vectorsVideos(readerobj2,78,vidHeight,vidWidth,G,numFrames,0,2,2);
%       vectorsVideos(readerobj3,78,vidHeight,vidWidth,G,numFrames,0,3,2);
%       vectorsVideos(readerobj4,78,vidHeight,vidWidth,G,numFrames,0,4,2);
%       vectorsVideos(readerobj5,78,vidHeight,vidWidth,G,numFrames,1,5,2);
%       vectorsVideos(readerobj6,78,vidHeight,vidWidth,G,numFrames,1,6,2);
%       vectorsVideos(readerobj7,78,vidHeight,vidWidth,G,numFrames,1,7,2);
       vectorsVideos(readerobj8,78,vidHeight,vidWidth,G,numFrames,1,8,2);
%       vectorsVideos(readerobj9,78,vidHeight,vidWidth,G,numFrames,1,9,2);
%       vectorsVideos(readerobj10,78,vidHeight,vidWidth,G,numFrames,1,10,2);
%       vectorsVideos(readerobj11,78,vidHeight,vidWidth,G,numFrames,1,11,2);
    
    
    
    
%         for j=1:1             
%             %vectorsVideos(readerobj1,j*22,vidHeight,vidWidth,G,numFrames,0,1,j);
%             vectorsVideos(readerobj1,j*34,vidHeight,vidWidth,G,numFrames,0,1,j);
%             %vectorsVideos(readerobj1,j*44,vidHeight,vidWidth,G,numFrames,0,1,j);
%             
%             %vectorsVideos(readerobj2,j*22,vidHeight,vidWidth,G,numFrames,0,2,j);
%             vectorsVideos(readerobj2,j*34,vidHeight,vidWidth,G,numFrames,0,2,j);
%             %vectorsVideos(readerobj2,j*44,vidHeight,vidWidth,G,numFrames,0,2,j);
%             
%             %vectorsVideos(readerobj3,j*22,vidHeight,vidWidth,G,numFrames,0,3,j);
%             vectorsVideos(readerobj3,j*34,vidHeight,vidWidth,G,numFrames,0,3,j);
%             %vectorsVideos(readerobj3,j*44,vidHeight,vidWidth,G,numFrames,0,3,j);
%             
%             %vectorsVideos(readerobj4,j*22,vidHeight,vidWidth,G,numFrames,0,4,j);
%             vectorsVideos(readerobj4,j*34,vidHeight,vidWidth,G,numFrames,0,4,j);
%             %vectorsVideos(readerobj4,j*44,vidHeight,vidWidth,G,numFrames,0,4,j);
%             
%             %vectorsVideos(readerobj5,j*22,vidHeight,vidWidth,G,numFrames,1,5,j);
%             vectorsVideos(readerobj5,j*34,vidHeight,vidWidth,G,numFrames,1,5,j);
%             %vectorsVideos(readerobj5,j*44,vidHeight,vidWidth,G,numFrames,1,5,j);
%             
%             %vectorsVideos(readerobj6,j*22,vidHeight,vidWidth,G,numFrames,1,6,j);
%             vectorsVideos(readerobj6,j*34,vidHeight,vidWidth,G,numFrames,1,6,j);
%             vectorsVideos(readerobj6,j*44,vidHeight,vidWidth,G,numFrames,1,6,j);
%             
%             %vectorsVideos(readerobj7,j*22,vidHeight,vidWidth,G,numFrames,1,7,j);
%             vectorsVideos(readerobj7,j*34,vidHeight,vidWidth,G,numFrames,1,7,j);
%             %vectorsVideos(readerobj7,j*44,vidHeight,vidWidth,G,numFrames,1,7,j);
%             
%             %vectorsVideos(readerobj8,j*22,vidHeight,vidWidth,G,numFrames,1,8,j);
%             vectorsVideos(readerobj8,j*34,vidHeight,vidWidth,G,numFrames,1,8,j);
%             %vectorsVideos(readerobj8,j*44,vidHeight,vidWidth,G,numFrames,1,8,j);
%             
%             %vectorsVideos(readerobj9,j*22,vidHeight,vidWidth,G,numFrames,1,9,j);
%             vectorsVideos(readerobj9,j*34,vidHeight,vidWidth,G,numFrames,1,9,j);
%             %vectorsVideos(readerobj9,j*44,vidHeight,vidWidth,G,numFrames,1,9,j);
%             
%             %vectorsVideos(readerobj10,j*22,vidHeight,vidWidth,G,numFrames,1,10,j);
%             vectorsVideos(readerobj10,j*34,vidHeight,vidWidth,G,numFrames,1,10,j);
%             %vectorsVideos(readerobj10,j*44,vidHeight,vidWidth,G,numFrames,1,10,j);
%             
%             vectorsVideos(readerobj11,j*22,vidHeight,vidWidth,G,numFrames,1,11,j);
%             vectorsVideos(readerobj11,j*34,vidHeight,vidWidth,G,numFrames,1,11,j);
%                vectorsVideos(readerobj11,j*44,vidHeight,vidWidth,G,numFrames,1,11,j);
% 
%        end
        
    
        
%      vectorsVideos(readerobj1,22,vidHeight,vidWidth,G,numFrames,0,1,1);
%     vectorsVideos(readerobj1,34,vidHeight,vidWidth,G,numFrames,0,1,1);
%     vectorsVideos(readerobj1,88,vidHeight,vidWidth,G,numFrames,0,1,2);
%     vectorsVideos(readerobj1,136,vidHeight,vidWidth,G,numFrames,0,1,2);
%      vectorsVideos(readerobj1,176,vidHeight,vidWidth,G,numFrames,0,1,3);
%     vectorsVideos(readerobj1,272,vidHeight,vidWidth,G,numFrames,0,1,3);
%    
%      vectorsVideos(readerobj2,22,vidHeight,vidWidth,G,numFrames,0,2,1);
% %     vectorsVideos(readerobj2,34,vidHeight,vidWidth,G,numFrames,0,2,1);
%      vectorsVideos(readerobj2,88,vidHeight,vidWidth,G,numFrames,0,2,2);
% %     vectorsVideos(readerobj2,136,vidHeight,vidWidth,G,numFrames,0,2,2);
%      vectorsVideos(readerobj2,176,vidHeight,vidWidth,G,numFrames,0,2,3);
%     vectorsVideos(readerobj2,272,vidHeight,vidWidth,G,numFrames,0,2,3);
%      
%      vectorsVideos(readerobj3,22,vidHeight,vidWidth,G,numFrames,0,3,1);
% %     vectorsVideos(readerobj3,34,vidHeight,vidWidth,G,numFrames,0,3,1);
%      vectorsVideos(readerobj3,88,vidHeight,vidWidth,G,numFrames,0,3,2);
% %     vectorsVideos(readerobj3,136,vidHeight,vidWidth,G,numFrames,0,3,2);
%      vectorsVideos(readerobj3,176,vidHeight,vidWidth,G,numFrames,0,3,3);
%     vectorsVideos(readerobj3,272,vidHeight,vidWidth,G,numFrames,0,3,3);
%      
%      vectorsVideos(readerobj4,22,vidHeight,vidWidth,G,numFrames,0,4,1);
% %     vectorsVideos(readerobj4,34,vidHeight,vidWidth,G,numFrames,0,4,1);
%      vectorsVideos(readerobj4,88,vidHeight,vidWidth,G,numFrames,0,4,2);
% %     vectorsVideos(readerobj4,136,vidHeight,vidWidth,G,numFrames,0,4,2);
%      vectorsVideos(readerobj4,176,vidHeight,vidWidth,G,numFrames,0,4,3);
%     vectorsVideos(readerobj4,272,vidHeight,vidWidth,G,numFrames,0,4,3);
%    
 %     vectorsVideos(readerobj5,22,vidHeight,vidWidth,G,numFrames,1,5,1);
% %     vectorsVideos(readerobj5,34,vidHeight,vidWidth,G,numFrames,1,5,1);
 %     vectorsVideos(readerobj5,88,vidHeight,vidWidth,G,numFrames,1,5,2);
% %     vectorsVideos(readerobj5,136,vidHeight,vidWidth,G,numFrames,1,5,2);
%      vectorsVideos(readerobj5,176,vidHeight,vidWidth,G,numFrames,1,5,3);
%     vectorsVideos(readerobj5,272,vidHeight,vidWidth,G,numFrames,1,5,3);
% 
%      vectorsVideos(readerobj6,22,vidHeight,vidWidth,G,numFrames,1,6,1);
% %     vectorsVideos(readerobj6,34,vidHeight,vidWidth,G,numFrames,1,6,1);
%      vectorsVideos(readerobj6,88,vidHeight,vidWidth,G,numFrames,1,6,2);
% %     vectorsVideos(readerobj6,136,vidHeight,vidWidth,G,numFrames,1,6,2);
%      vectorsVideos(readerobj6,176,vidHeight,vidWidth,G,numFrames,1,6,3);
%     vectorsVideos(readerobj6,272,vidHeight,vidWidth,G,numFrames,1,6,3);

%      vectorsVideos(readerobj7,22,vidHeight,vidWidth,G,numFrames,1,7,1);
% %     vectorsVideos(readerobj7,34,vidHeight,vidWidth,G,numFrames,1,7,1);
%      vectorsVideos(readerobj7,88,vidHeight,vidWidth,G,numFrames,1,7,2);
% %     vectorsVideos(readerobj7,136,vidHeight,vidWidth,G,numFrames,1,7,2);
%      vectorsVideos(readerobj7,176,vidHeight,vidWidth,G,numFrames,1,7,3);
%     vectorsVideos(readerobj7,272,vidHeight,vidWidth,G,numFrames,1,7,3);

%      vectorsVideos(readerobj8,22,vidHeight,vidWidth,G,numFrames,1,8,1);
% %     vectorsVideos(readerobj8,34,vidHeight,vidWidth,G,numFrames,1,8,1);
%      vectorsVideos(readerobj8,88,vidHeight,vidWidth,G,numFrames,1,8,2);
% %     vectorsVideos(readerobj8,136,vidHeight,vidWidth,G,numFrames,1,8,2);
%      vectorsVideos(readerobj8,176,vidHeight,vidWidth,G,numFrames,1,8,3);
%     vectorsVideos(readerobj8,272,vidHeight,vidWidth,G,numFrames,1,8,3);

%      vectorsVideos(readerobj9,22,vidHeight,vidWidth,G,numFrames,1,9,1);
% %     vectorsVideos(readerobj9,34,vidHeight,vidWidth,G,numFrames,1,9,1);
%      vectorsVideos(readerobj9,88,vidHeight,vidWidth,G,numFrames,1,9,2);
% %     vectorsVideos(readerobj9,136,vidHeight,vidWidth,G,numFrames,1,9,2);
%      vectorsVideos(readerobj9,176,vidHeight,vidWidth,G,numFrames,1,9,3);
%     vectorsVideos(readerobj9,272,vidHeight,vidWidth,G,numFrames,1,9,3);
% 
%      vectorsVideos(readerobj10,22,vidHeight,vidWidth,G,numFrames,1,10,1);
% %     vectorsVideos(readerobj10,34,vidHeight,vidWidth,G,numFrames,1,10,1);
%      vectorsVideos(readerobj10,88,vidHeight,vidWidth,G,numFrames,1,10,2);
% %     vectorsVideos(readerobj10,136,vidHeight,vidWidth,G,numFrames,1,10,2);
%      vectorsVideos(readerobj10,176,vidHeight,vidWidth,G,numFrames,1,10,3);
%     vectorsVideos(readerobj10,272,vidHeight,vidWidth,G,numFrames,1,10,3);
% 
%      vectorsVideos(readerobj11,22,vidHeight,vidWidth,G,numFrames,1,11,1);
% %     vectorsVideos(readerobj11,34,vidHeight,vidWidth,G,numFrames,1,11,1);
%      vectorsVideos(readerobj11,88,vidHeight,vidWidth,G,numFrames,1,11,2);
% %     vectorsVideos(readerobj11,136,vidHeight,vidWidth,G,numFrames,1,11,2);
%      vectorsVideos(readerobj11,176,vidHeight,vidWidth,G,numFrames,1,11,3);
%     vectorsVideos(readerobj11,272,vidHeight,vidWidth,G,numFrames,1,11,3);

end

function vectorsVideos(readerobj,radius,Height,Width,g,nFrames,Togray,v,condicao)


            mainPath = strcat('/Users/joaodornas/Documents/_Research/Cenas Naturais/VIDEOS/VIDEO-',int2str(v),'/');
                
            vectorFolder = strcat(int2str(condicao),'_','v',int2str(v),'_a_forwards_',int2str(radius));
            
            mkdir(mainPath,vectorFolder);
                
            %%% CRC2,17 nCRC7,84

            %CRC = 57*tan(2.17*pi/180); diameter == 68 pixels
            %nCRC = 57*tan(7.84*pi/180); diameter == 300 pixels 

            movgray(1:nFrames) = ...
            struct('cdata', zeros(Height, Width),...
                   'colormap', gray);

            movgraycrop(1:nFrames) = ...
            struct('cdata', zeros(2*radius, 2*radius),...
                   'colormap', gray);

            mask = imcircle(2*radius);
                
            hStart = (Height/2) - radius;
            hEnd = (Height/2) + radius - 1;
            wStart = (Width/2) - radius;
            wEnd = (Width/2) + radius - 1;
            
            h = fspecial('gaussian',g,g);

             writerObj = VideoWriter(strcat(mainPath,'video_out-v',int2str (v),'-r',int2str(radius),'.avi'));
             writerObj.FrameRate = 30;
             open(writerObj);
            
            vidFrames = read(readerobj);

            for k = 1 : nFrames

                movgray(k).cdata = vidFrames(:,:,:,k);
                
                if Togray == 1
                    
                    movgray(k).cdata = rgb2gray(movgray(k).cdata);
                    
                end
                
                %imwrite(movgray(k).cdata,'framegray.jpg');

                movgraycrop(k).cdata = movgray(k).cdata(hStart:hEnd,wStart:wEnd);
                
                movgraycrop(k).cdata(movgraycrop(k).cdata==0) = 1;

                %imwrite(movgraycrop(k).cdata,strcat('framecrop-',int2str(radius),'.jpg'));

                movgraycrop(k).cdata(~mask) = 255;

                %imwrite(movgraycrop(50).cdata,strcat('framemask-',int2str(radius),'.jpg'));

                 movgraycrop(k).cdata = vertcat(ones(radius,2*radius)*255,movgraycrop(k).cdata);
                 
                 movgraycrop(k).cdata = vertcat(movgraycrop(k).cdata,ones(radius,2*radius)*255);
                 
                 movgraycrop(k).cdata = horzcat(movgraycrop(k).cdata,ones(4*radius,radius)*255);
                 
                 movgraycrop(k).cdata = horzcat(ones(4*radius,radius)*255,movgraycrop(k).cdata);
  
                 movgraycrop(k).cdata = imfilter(movgraycrop(k).cdata,h);
                 
                 length = size(movgraycrop(k).cdata,2);
                 
%                  if radius*2 > 500 
%                      
%                      increase = 5;
%                      
%                  else
%                      
%                      increase = 4;
%                      
%                  end

                 increase = 4;
                 
                 start_ = length/2 - radius - increase;
                 end_ = length/2 + radius - 1 + increase;
                 
                 movgraycrop(k).cdata = movgraycrop(k).cdata(start_:end_,start_:end_);
                
                 imwrite(movgraycrop(50).cdata,strcat(mainPath,'framefilter-','video-',int2str(v),'-',int2str(radius),'-',int2str(50),'.jpg'));
                 
                 writeVideo(writerObj,movgraycrop(k).cdata);

            end

            %movie2avi(movgraycrop,strcat('video_out-v',int2str (v),'-r',int2str(radius),'.avi'),'compression','None');

            %close(writerObj);
            
            movgraycropReverse = fliplr(movgraycrop);

            clear movgray;

%             k = 1;        
%             while (k < nFrames + 1)
% 
%                 X = reshape(movgraycrop(k).cdata.',[],1);
% 
%                 X = X.';
% 
% %                 XR = reshape(movgraycropReverse(k).cdata.',[],1);
% % 
% %                 XR = XR.';
% 
%                 if k < 10
%                     frameindex = strcat('00', int2str(k));
%                 elseif k < 100
%                     frameindex = strcat('0', int2str(k));
% %                 elseif k < 1000
% %                     frameindex = strcat('0', int2str(k));
%                 else
%                     frameindex = int2str(k);
%                 end
% 
%                 filename = strcat(vectorFolder,'_', frameindex);
%                 
%                 fullPath = strcat(mainPath,vectorFolder,'/',filename);
% 
%                 dlmwrite(strcat(fullPath, '.txt'), X, 'delimiter', '\t', '-append', 'newline', 'pc');
%                 
% %               pathr = strcat(condicao,'_','v',int2str(v),'_b_backwards_',int2str(radius),'_', frameindex);
% 
% %               dlmwrite(strcat(pathr, '.txt'), XR, 'delimiter', '\t', '-append', 'newline', 'pc');
% 
%                 k = k + 1;
% 
%             end

            clear movgraycrop;
            clear movgraycropRerverse;

    end
    



