

prefix = 'MOT-Analyze-V2-10-m-';

uniques = load('uniqueKEmin.mat');

idx_uniques_EMmin = uniques.idx_MOT_unique_EMmin;
idx_uniques_EMmax = uniques.idx_MOT_unique_EMmax;

for i=1:length(idx_uniques_EMmin)
    
    disp(strcat('min:',int2str(i)));

    datafile = load(strcat(prefix,int2str(idx_uniques_EMmin(i)),'.mat'));
   
    kEMmin = datafile.kEMmin;
    
    EmotionMin(i) = datafile.Emotion(kEMmin);
    
    clear datafile
    
end

for i=1:length(idx_uniques_EMmax)
    
    disp(strcat('max:',int2str(i)));

    datafile = load(strcat(prefix,int2str(idx_uniques_EMmax(i)),'.mat'));
   
    kEMmax = datafile.kEMmax;
    
    EmotionMax(i) = datafile.Emotion(kEMmax);
    
    clear datafile
    
end

if length(idx_uniques_EMmin) > length(idx_uniques_EMmax)
    
    A = idx_uniques_EMmin;
    B = idx_uniques_EMmax;
    
else
    
    B = idx_uniques_EMmin;
    A = idx_uniques_EMmax;
    
end

idx_uniques_EM = A(ismember(A,B));

maxEntropy = max(EmotionMax);
minEntropy = min(EmotionMin);

nBins = 3;

[X,Y] = hist(round(minEntropy):round(maxEntropy)-1,nBins);

info.binsCenter(1) = round(minEntropy);
info.binsCenter(2) = Y(2);
info.binsCenter(3) = round(maxEntropy);
info.minEntropy = round(minEntropy);
info.maxEntropy = round(maxEntropy)-1;

nK = 12000;

allEmotion = zeros(length(idx_uniques_EM),nK);

for i=1:length(idx_uniques_EM)
    
    disp(strcat('unique:',int2str(i)));

    datafile = load(strcat(prefix,int2str(idx_uniques_EM(i)),'.mat'));
   
    Emotion = datafile.Emotion;
    
    allEmotion(i,:) = Emotion;
    
    clear datafile
    
end

nCenters = length(info.binsCenter);

nSamples = length(idx_uniques_EM);

Time = 10;
FrameRate = 60;
maxNK = 12000;
sFrames = 3;
nWindow = 4;

limitFramesInKs = floor((Time*FrameRate/2)/sFrames) + nWindow;

for nC=1:nCenters
    
   for i=1:length(idx_uniques_EM)
      
       line = abs(allEmotion(i,:) - info.binsCenter(nC));
       idx_min = find(line == min(line));
       
       ki(i) = idx_min(1);
       ei(i) = allEmotion(i,idx_min(1));
       
   end
    
   nearestValues = ei - info.binsCenter(nC);
   [X,I] = sort(nearestValues);
   
   allKi = ki(I);
   
   idx_to_remove = [];
   remove = false;
   for kall=1:length(allKi)
       
       if (allKi(kall) <= limitFramesInKs) || ((maxNK - allKi(kall)) <= limitFramesInKs)
           
           idx_to_remove = [idx_to_remove, kall];
           
           remove = true;
           
       end
       
   end
   
   if remove; I(idx_to_remove) = []; end
   
   if nSamples > length(I); nSamples = length(I); end
   
   info.centers(nC).binsMOTidx = idx_uniques_EM(I(1:nSamples));
   info.centers(nC).binsMOTk = ki(I(1:nSamples));
   info.centers(nC).entropies = ei(I(1:nSamples)); 
   
end

save('infoBinsCentersMOTMotionEntropy.mat','info'),