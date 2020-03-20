
movie = VideoReader('entropies.avi');

idx_frame(1) = kECmax;
str_frame{1} = 'kECmax';
idx_frame(2) = kECmin;
str_frame{2} = 'kECmin';
idx_frame(3) = kEMmax;
str_frame{3} = 'kEMmax';
idx_frame(4) = kEMmin;
str_frame{4} = 'kEMmin';
idx_frame(5) = kEBmax;
str_frame{5} = 'kEBmax';
idx_frame(6) = kEBmin;
str_frame{6} = 'kEBmin';
idx_frame(7) = kEQmax;
str_frame{7} = 'kEQmax';
idx_frame(8) = kEQmin;
str_frame{8} = 'kEQmin';
idx_frame(9) = kVCmax-1;
str_frame{9} = 'kVCmax';
idx_frame(10) = kVCmin;
str_frame{10} = 'kVCmin';
idx_frame(11) = kVMmax;
str_frame{11} = 'kVMmax';
idx_frame(12) = kVMmin;
str_frame{12} = 'kVMmin';
idx_frame(13) = kVBmax;
str_frame{13} = 'kVBmax';
idx_frame(14) = kVBmin;
str_frame{14} = 'kVBmin';
idx_frame(15) = kVQmax;
str_frame{15} = 'kVQmax';
idx_frame(16) = kVQmin;
str_frame{16} = 'kVQmin';

idx_frame = idx_frame.*60;

for k=1:length(idx_frame)
    
    frame = read(movie,idx_frame(k));
    imwrite(frame,strcat(str_frame{k},'.jpg'),'JPEG');
    
    clear frame;
    
end

close(movie);
