function rfmapstimulusBar(side,spotSide,background,nCondicoes)

tic

a = 4*spotSide;
b = spotSide;

nOrientations = 360/22.5;

nFrames = ( (side - 2*(2*a - 1) - mod(side - 2*(2*a - 1),4))/4 )^2 * 2;

for D=1:1
    
    stimulusMap(1:side,1:side,1:nFrames) = background;
        
    F = 1;
    
        for x=2*a:4:(side - 2*a)
            
            for y=2*a:4:(side - 2*a)
                
                [X Y] = calculateEllipse(x, y, a, b, (D-1)*22.5, 360);

                X = round(X);
                Y = round(Y);

                Y = abs(Y);
                X = abs(X);

                nPoints = size(X,1);

                for i=1:nPoints

                    if Y(i) == 0; Y(i) = 1; end

                    if X(i) == 0; X(i) = 1; end

                    stimulusMap(Y(i),X(i),F) = 1;
                    stimulusMap(Y(i),X(i),F+1) = 255;

                end

                coordenadas(D,F,1,:) = X;
                coordenadas(D,F,2,:) = Y;
                coordenadas(D,F,3,:) = x;
                coordenadas(D,F,4,:) = y;

                coordenadas(D,F+1,1,:) = X;
                coordenadas(D,F+1,2,:) = Y;
                coordenadas(D,F+1,3,:) = x;
                coordenadas(D,F+1,4,:) = y;

                
                for i=1:side

                    for j=1:side


                        if  ((i-x)*cos(-(D-1)*22.5*pi/180)-(j-y)*sin(-(D-1)*22.5*pi/180)).^2/a^2 + ((i-x)*sin(-(D-1)*22.5*pi/180)+(j-y)*cos(-(D-1)*22.5*pi/180)).^2/b^2 <= 1

                            stimulusMap(j,i,F) = 1;
                            stimulusMap(j,i,F+1) = 255;

                        end

                    end
              
                end

                F = F + 2;
                
            end
        end
        
        filepathDirection = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Receptive-Field/bar-oriented/direction-',int2str(D),'/');
        
        mkdir(filepathDirection);
        
        save(strcat(filepathDirection,'Coordenadas-ellipse-position-direction-',int2str(D)),'coordenadas');
    
    for C=1:nCondicoes
        
        folder = strcat(int2str(D),'_','condition-',int2str(C),'-stimulusMap-Bar-side-',int2str(side),'-spotSide-',int2str(spotSide));
        
        filepathVideos = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Receptive-Field/bar-oriented/direction-',int2str(D),'/','videos','/');
        
        mkdir(filepathVideos);
        
        writerObj = VideoWriter(strcat(filepathVideos,'stimulusMap-Bar-side-',int2str(side),'-spotSide-',int2str(spotSide),'-direction-',int2str(D),'-background-',int2str(background),'-condition-',int2str(C),'.avi'));
        writerObj.FrameRate = 30;
        open(writerObj);
    
        order = 0;
        
        while size(order,2) < nFrames + 1

            next = randi(nFrames,1);

            exist = 0;

            for O=1:size(order,2)

                if next == order(O);

                    exist = 1;

                end

            end

            if exist == 0

                order = [order next];

            end

        end

        order(1) = [];

        newStimulusMap = stimulusMap(:,:,order);

        filepathFrameOrder = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Receptive-Field/bar-oriented/direction-',int2str(D),'/','frame-order','/');
         
        mkdir(filepathFrameOrder);
        
        save(strcat(filepathFrameOrder,'Frame-order-direction-',int2str(D),'-condition-',int2str(C)),'order');

        filepathVectors = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Receptive-Field/bar-oriented/orientation-',int2str(D),'/','vectors','/');
            
        mkdir(filepathVectors,folder);
        
        path_to_info = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Receptive-Field/bar-oriented/_info.txt');
        new_path_to_info = strcat(filepathVectors,folder);
        copyfile(path_to_info,new_path_to_info);
        
        for F=1:nFrames

            frame = newStimulusMap(:,:,F);

            g = 2;

            h = fspecial('gaussian',g,g);

            newFrame = imfilter(frame,h);

            [G,newFrame] = gaborfilter(newFrame,2,4,16,pi/3);

            newFrame = real(newFrame);
            newFrame = uint8(newFrame);

            writeVideo(writerObj,uint8(newFrame));

            if F < 10
                %frameindex = strcat('0000', int2str(F));
            %elseif F < 100
                frameindex = strcat('000', int2str(F));
            elseif F < 100
                frameindex = strcat('00', int2str(F));
            elseif F < 1000 
                frameindex = strcat('0', int2str(F));
            else
                frameindex = int2str(F);
            end

            X = reshape(newFrame.',[],1);

            X = X.';
            
            framePath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Receptive-Field/bar-oriented/direction-',int2str(D),'/','vectors','/',int2str(D),'_','condition-',int2str(C),'-stimulusMap-Bar-side-',int2str(side),'-spotSide-',int2str(spotSide),'/',int2str(D),'_','condition-',int2str(C),'-stimulusMap-Bar-side-',int2str(side),'-spotSide-',int2str(spotSide));
            
            dlmwrite(strcat(framePath,'_', frameindex, '.txt'), X, 'delimiter', '\t', '-append', 'newline', 'pc');

        end

    end
  
end

close(writerObj);

toc

end