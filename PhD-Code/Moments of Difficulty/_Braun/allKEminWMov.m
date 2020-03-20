clear all

kemin = zeros(1,100);

for i=1:100
    
        
      if i > 1; load('allKEminWMov.mat'); end

      load(strcat('MOT-Mov-10-s-',int2str(i)));

      kemin(i) = MOT.kemin;
      
      save('allKEminWMov','kemin');
      
      clear all

end

clear all