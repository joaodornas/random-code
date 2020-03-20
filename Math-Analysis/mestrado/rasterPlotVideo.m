%function rasterPlotVideo(radius,nFrames,width,height)

%%%%%%%%%%%%%%%%%%%    FORWARDS-BACKWARDS

%Registro 2 (Arnaldo) - 27/08/12 - 10:00 - 1 protocolo (v?lido) 

%NSP005a01 v3
%centro X: 523 Y: 274

%rasterPlot(523,274,video3,3,'005a01');


%Registro 3 (Arnaldo) - 29/08/12 - 10:00 - 3 protocolos (v?lidos)

%NSP006a01&a02 v3,v4
%centro X: 604, Y: 402 

%rasterPlot(604,402,video3,3,'006a01');
%rasterPlot(604,402,video4,4,'006a02');


%NSP006b01 v1
%centro X: 587, Y: 403

%rasterPlot(587,403,video1,1,'006b01');



%Registro 4 (Arnaldo) - 30/08/12 - 09:30 - 3 protocolos (v?lidos)

%NSP007a01&a02&a03 v4,v1,v3
%centro X: 598 Y: 435
%centro X: 679 Y: 496

%rasterPlot(598,435,video4,4,'007a01');
%rasterPlot(598,435,video1,1,'007a02');
%rasterPlot(598,435,video3,3,'007a03');

% rasterPlot(679,496,video1,1,'007a01');
% rasterPlot(679,496,video4,4,'007a02');


%Registro 5 (Arnaldo) - 31/08/12 - 10:00 - 3 protocolos (v?lidos)

%NSP008a01&a02&a03 v1,v3,v4
%centro E1: X: 602, Y: 362 / E: X: 569, Y: 356

%rasterPlot(602,362,video1,1,'008a01');
%rasterPlot(602,362,video3,3,'008a02');
%rasterPlot(602,362,video4,4,'008a03');

%rasterPlot(569,356,video1,1,'008a01');
%rasterPlot(569,356,video3,3,'008a02');
%rasterPlot(569,356,video4,4,'008a03');


%Registro 6 (Arnaldo) - 03/09/12 - 10:45 - 2 protocolos (v?lidos)

%NSP009a01 v2
%centro X: 164, Y: 378

%rasterPlot(164,378,video2,2,'009a01');


%NSP009b01 v4
%centro X: 587, Y: 370

%rasterPlot(587,370,video4,4,'009b01');

%Registro 7 (Arnaldo) - 04/09/12 - 11:30 - 2 protocolos (v?lidos)

%NSP010a01&a02 v4
%centro X: 724, Y: 443

%rasterPlot(724,443,video4,4,'010a02');


%Registro 8 (Arnaldo) - 05/09/12 - 09:30 - 2 protocolos (v?lidos)

%NSP011a01&a02 v3,v4
%centro X: 689, Y: 442

%rasterPlot(689,442,video3,3,'011a01');
%rasterPlot(689,442,video4,4,'011a02');

%Registro 9 (Arnaldo) - 06/09/12 - 10:30 - 2 protocolos (v?lidos)

%NSP012a01&a02 v3,v4
%centro X: 683, Y: 458

%rasterPlot(683,458,video3,3,'012a01');
%rasterPlot(683,458,video4,4,'012a02');

%Registro 10 (Pantro) - 16/10/12 - 14:00 - 1 protocolo (v?lido)

%NSP013a01 v1
%centro X: 560, Y: 408

%rastertPlot(560,408,video1,1,'013a01');


function t = rasterPlotVideo(X,Y,vn,registro,width,height,nFrames,radius)
    
    g = 4;

    video = VideoReader(['_MATRIZ-v' int2str(vn) '-8bit-' int2str(width) 'x' int2str(height) '.avi'], 'tag', strcat('myreader',int2str(vn)))

    movgray(1:nFrames) = struct('cdata',zeros(height, width),'colormap', gray);
    movgraycrop(1:nFrames) = struct('cdata', zeros(2*radius, 2*radius),'colormap',gray);

    writerObj = VideoWriter(strcat('NSP',registro,'-v',int2str(vn),'-x', int2str(X),'-y', int2str(Y),'.avi'));
    writerObj.FrameRate = 30;
    open(writerObj);
    vidFrames = read(video);

    mask = imcircle(2*radius);

    hStart = Y - radius;
    hEnd = Y + radius - 1;
    wStart = X - radius;
    wEnd = X + radius - 1;

            for k = 1 : nFrames

                movgray(k).cdata = vidFrames(:,:,:,k);

                %imwrite(movgray(k).cdata,'framegray.jpg');

                movgraycrop(k).cdata = movgray(k).cdata(hStart:hEnd,wStart:wEnd);

                %imwrite(movgraycrop(k).cdata,strcat('framecrop-',int2str(radius),'.jpg'));

                movgraycrop(k).cdata(~mask) = 0;

                h = fspecial('gaussian',g,g);

                movgraycrop(k).cdata = imfilter(movgraycrop(k).cdata,h);

                imwrite(movgraycrop(k).cdata,strcat('NSP',registro,'-v',int2str(vn),'-x', int2str(X),'-y', int2str(Y),'-frame-',int2str(k),'.jpg'));

                %imwrite(movgraycrop(50).cdata,strcat('framefilter-',int2str(radius),'-',int2str(50),'.jpg'));

                writeVideo(writerObj,movgraycrop(k).cdata);

            end

    %movie2avi(movgraycrop,strcat('video_out-v',int2str (v),'-r',int2str(radius),'.avi'),'compression','None');

    close(writerObj);

    clear movgraycrop;
    clear movgray;
    clear writerObj;

end