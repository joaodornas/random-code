clear all

whoMin = zeros(1,100);
whoMax = zeros(1,100);

minK = 1;
maxK = 1;

for i=1:100
    
        
      if i > 1; load('whoIsMinMaxColor.mat'); end

      load(strcat('MOT-TrajectoriesAnalysis-5-m-',int2str(i)));

      minLcolor = Lcolor(kVCmin);
      maxHcolor = Hcolor(kVCmax);

      if (1 - minLcolor) > maxHcolor

        whoMin(minK) = i;
        minK = minK + 1;
        kVC(i) = kVCmin;
        
      else

        whoMax(maxK) = i;
        maxK = maxK + 1;
        kVC(i) = kVCmax;
        
      end
      
      whoMin(whoMin == 0) = [];
      whoMax(whoMax == 0) = [];
      
      save('whoIsMinMaxColor','whoMin','whoMax','kVC','minK','maxK');
      
      clear all

end

clear all