function dong(v,nFrames,width,height,radius)

tic

videoin = mmreader(['_MATRIZ-v' int2str(v) '-8bit-' int2str(width) 'x' int2str(height) '.avi']);

w = 4;

for i=1:nFrames
    
   frame = read(videoin,i);
   
   frame = frame(:,:,1);

   mask = imcircle(2*radius);
    
   hStart = (height/2) - radius;
   hEnd = (height/2) + radius - 1;
   wStart = (width/2) - radius;
   wEnd = (width/2) + radius - 1;
   
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
   
   frame = reshape(frame.',[],1).';
    
   frames(:,i) = frame; 
    
end

pixels = size(frames,1);

R = zeros(pixels,nFrames);

for y=1:pixels
        
    for t=1:nFrames-1
           
        for yy=1:pixels-y
                    
            for tt=1:nFrames-t
                        
                R(y,t) = R(y,t) + ( ( frames(y + yy,t + tt).*frames(y,t) ) );
            
            end
            
        end
        
    end    
    
end

R = R ./ (pixels.*nFrames) ;

RF = zeros(pixels,nFrames);

for fy=1:pixels
        
    for w=1:nFrames
                
        for j=1:pixels-fy
                    
            for k=1:nFrames-w
                        
                RF(fy,w) = RF(fy,w) + ( R(j,k).*exp(sqrt(-1)*2*pi*(j*fy + k*w)) );
            
             end
                    
         end
                
    end
    
end

RF = Rf ./ (pixels.*nFrames) ;

A = abs(RF);

save('A','A');

plot(log10(A(1:pixels,1,1)),-log10((1:pixels)/pixels));

toc

end

