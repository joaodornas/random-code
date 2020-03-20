clear all

entropyMax = zeros(1,100);
entropyMin = zeros(1,100);

keminMax = zeros(1,100);
keminMin = zeros(1,100);

for i=1:100
    
        
      if i > 1; load('MinMaxEntropyColor.mat'); end

      load(strcat('MOT-TrajectoriesAnalysis-5-m-',int2str(i)));

      entropyMax(i) = Ecolor(kECmax);
      entropyMin(i) = Ecolor(kECmin);

      keminMax(i) = kECmax;
      keminMin(i) = kECmin;
      
      save('MinMaxEntropyColor','entropyMin','entropyMax','keminMax','keminMin');
      
      clear all

end

clear all