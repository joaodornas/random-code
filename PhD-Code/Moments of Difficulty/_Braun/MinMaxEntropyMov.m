clear all

entropyMax = zeros(1,100);
entropyMin = zeros(1,100);

keminMax = zeros(1,100);
keminMin = zeros(1,100);

for i=1:100
    
        
      if i > 1; load('MinMaxEntropyMov.mat'); end

      load(strcat('MOT-TrajectoriesAnalysis-5-m-',int2str(i)));

      entropyMax(i) = Emotion(kEMmax);
      entropyMin(i) = Emotion(kEMmin);

      keminMax(i) = kEMmax;
      keminMin(i) = kEMmin;
      
      save('MinMaxEntropyMov','entropyMin','entropyMax','keminMax','keminMin');
      
      clear all

end

clear all