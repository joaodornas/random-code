function h = mapCalibrate

nFrames = 100;

picture = ones(768,1024)*255;

picture(1:4,1:4) = 0;
picture(188:196,1:4) = 0;
picture(380:392,1:4) = 0;
picture(572:580,1:4) = 0;
picture(764:768,1:4) = 0;

picture(1:4,252:260) = 0;
picture(188:196,252:260) = 0;
picture(380:392,252:260) = 0;
picture(572:580,252:260) = 0;
picture(764:768,252:260) = 0;

picture(1:4,508:516) = 0;
picture(188:196,508:516) = 0;
picture(380:392,508:516) = 0;
picture(572:580,508:516) = 0;
picture(764:768,508:516) = 0;

picture(1:4,764:772) = 0;
picture(188:196,764:772) = 0;
picture(380:392,764:772) = 0;
picture(572:580,764:772) = 0;
picture(764:768,764:772) = 0;

picture(1:4,1020:1024) = 0;
picture(188:196,1020:1024) = 0;
picture(380:392,1020:1024) = 0;
picture(572:580,1020:1024) = 0;
picture(764:768,1020:1024) = 0;

colormap gray;

imshow(picture);

imwrite(picture,'picture.jpg');

% 
%             k = 1;        
%             while (k < nFrames + 1)
% 
%                 X = reshape(picture.',[],1);
% 
%                 X = X.';
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
%                     path = strcat('mapRF-1024-720','_', frameindex);
% 
%                     dlmwrite(strcat(path, '.txt'), X, 'delimiter', '\t', '-append', 'newline', 'pc');
% 
%                     k = k + 1;
% 
%             end


end