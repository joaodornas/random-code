function matrix2vectorsMeanGray(FrameRate,nFrames,width,height,videoN,removeTransparent,framePace)



for v=1:1
    
    movieMatrix(1:nFrames) = struct('cdata',zeros(height,width));
    
    videoin = mmreader(['_MATRIZ-v' int2str(v) '-8bit-' int2str(width) 'x' int2str(height) '.avi']);

    for k = 1 : nFrames

        frame = read(videoin,k);
        
         if removeTransparent == 1
           
            frame(frame==0) = 1;
            
        end

        movieMatrix(k).cdata = frame;
        
        frameMeanGray(k) = mean2(movieMatrix(k).cdata);
        
    end
    
    movieMeanPixel = uint8(mean2(frameMeanGray));
    
    frameGray(1) = struct('cdata', ones(height,width,'uint8'));
    
    frameGray(1).cdata = frameGray(1).cdata.*movieMeanPixel;
    
    frameGray(1).cdata = round(frameGray(1).cdata);

    imwrite(frameGray(1).cdata, strcat('v',int2str(v),'.jpg'));
            
%     k = 1;
%     f = 1;
% 
%     while (f < (nFrames + 1))
% 
%         X = reshape(frameGray(1).cdata.',[],1);
% 
%         X = X.';
% 
%          for i=1:framePace
%              
%             if k < 10
%                 frameindex = strcat('00', int2str(k));
%             elseif k < 100
%                 frameindex = strcat('0', int2str(k));
%             else
%                 frameindex = int2str(k);
%             end
% 
%             path = strcat('v',int2str(v),'_c_gray_',int2str(width),'_', frameindex);
% 
%             dlmwrite(strcat(path, '.txt'), X, 'delimiter', '\t', '-append', 'newline', 'pc');
% 
%             k = k + 1;
%             
%          end
%          
%          f = f + 1;
% 
%     end
%     
%     clear movieMatrix;
%   
%     clear videoin;    
    
end



end

