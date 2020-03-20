function [ForwardMovieCRF, BackwardMovieCRF] = ForBackCRFmovie(name,v,nFrames,width,height,radius,save,mainPath)

tic

videoin = VideoReader(['_MATRIZ-v' int2str(v) '-8bit-' int2str(1024) 'x' int2str(720) '.avi']);

w = 4;

% if strcmp(OS,'Mac')
%     
%     mainPath = strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Dong/');
%     
% elseif strcmp(OS,'Win')
%     
%     mainPath = strcat('Z:\DATA\Forward-Backward\Power-Law\Dong\');
%     
% end

if save == 1
    
    writerObj = VideoWriter(strcat(mainPath,name,'-video_out-v',int2str(v),'-r',int2str(radius),'-w',int2str(width),'-h',int2str(height),'-Forward.avi'));
    writerObj.FrameRate = 30;
    open(writerObj);

end

for i=1:nFrames
    
   frame = read(videoin,i);
   
   frame = frame(:,:,1);

   mask = imcircle(2*radius);
    
%    hStart = (height/2) - radius;
%    hEnd = (height/2) + radius - 1;
%    wStart = (width/2) - radius;
%    wEnd = (width/2) + radius - 1;

   hStart = height - radius;
   hEnd = height + radius - 1;
   wStart = width - radius;
   wEnd = width + radius - 1;
   
   frame = frame(hStart:hEnd,wStart:wEnd);
    
   h = fspecial('gaussian',w,w);
   
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
    
   ForwardMovieCRF(:,:,i) = frame(:,:);

   if save == 1
       
        writeVideo(writerObj,frame);  
   
   end
       
end

if save == 1
    
    close(writerObj);
    
end

BackwardMovieCRF = flipdim(ForwardMovieCRF,3);


if save == 1
    
    writerObj = VideoWriter(strcat(mainPath,name,'-video_out-v',int2str (v),'-r',int2str(radius),'-w',int2str(width),'-h',int2str(height),'-Backward.avi'));
    writerObj.FrameRate = 30;
    open(writerObj);

end

for i=1:nFrames
    
    frame(:,:) = BackwardMovieCRF(:,:,i);

    if save == 1
        
        writeVideo(writerObj,frame);
        
    end

end

toc

end

