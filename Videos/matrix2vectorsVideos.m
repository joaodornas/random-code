function matrix2vectorsVideos(FrameRate,nFrames,width,height,videoN,framePace,reverseRand,vectors)



for v=1:videoN
    
    movieMatrix(1:nFrames) = struct('cdata',zeros(height,width));
    
    movieMatrixReverse(1:nFrames) = struct('cdata',zeros(height,width));
   
    videoin = mmreader(['_MATRIZ-v' int2str(v) '-8bit-' int2str(width) 'x' int2str(height) '.avi']);

    videoout = VideoWriter(['v' int2str(v) '-8bit-' int2str(width) 'x' int2str(height) ' - direct']);
    
    videoout.FrameRate = FrameRate;
    
    open(videoout);

    for k = 1 : nFrames

        frame = read(videoin,k);

        writeVideo(videoout, frame);
        
        movieMatrix(k).cdata = frame;
        
    end

   close(videoout);
    
   if reverseRand == 1
       
       movieMatrixReverse = fliplr(movieMatrix);
       
       name = 'reverse';
       
   else
       
       vector = randperm(length(movieMatrix));
       
       movieMatrixReverse = movieMatrix(vector);
       
       name = 'rand';
       
       save(['v' int2str(v) '-8bit-' int2str(width) 'x' int2str(height) ' - ' name '-frameOrder'],'vector');
       
   end
    
   videooutreverse = VideoWriter(['v' int2str(v) '-8bit-' int2str(width) 'x' int2str(height) ' - ' name]);
    
   videooutreverse.FrameRate = FrameRate;
    
   open(videooutreverse);
    
    for k = 1 : nFrames

        frame = movieMatrixReverse(k).cdata;
        
        writeVideo(videooutreverse, frame);

    end

   close(videooutreverse);
    
   if vectors == 1
       
        k = 1;
        f = 1;

        while (f < (nFrames + 1))

            X = reshape(movieMatrix(f).cdata(:,:,1).',[],1);

            X = X.';

            XR = reshape(movieMatrixReverse(f).cdata(:,:,1).',[],1);

            XR = XR.';

    %         if removeTransparent == 1
    %            
    %             X(X==0) = 1;
    %             XR(XR==0) = 1;
    %             
    %         end

                for i=1:framePace

                    if k < 10

                        frameindex = strcat('00', int2str(k));

                    elseif k < 100

                        frameindex = strcat('0', int2str(k));

                    else

                        frameindex = int2str(k);

                    end

                    path = strcat('v',int2str(v),'_a_direto_',int2str(width),'_', frameindex);

                    pathr = strcat('v',int2str(v),'_b_',name,'_',int2str(width),'_', frameindex);

                    dlmwrite(strcat(path, '.txt'), X, 'delimiter', '\t', '-append', 'newline', 'pc');

                    dlmwrite(strcat(pathr, '.txt'), XR, 'delimiter', '\t', '-append', 'newline', 'pc');

                    k = k + 1;

                end

                f = f + 1;

        end
        
   end
    
   clear movieMatrix;
   clear movieMatrixReverse;
    
   clear videoin;    
    
end



end

