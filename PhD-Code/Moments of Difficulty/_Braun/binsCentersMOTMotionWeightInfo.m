

prefix = 'MOT-Analyze-V2-10-m-';

uniques = load('uniqueKEmin.mat');

idx_uniques_VMmin = uniques.idx_MOT_unique_VMmin;
idx_uniques_VMmax = uniques.idx_MOT_unique_VMmax;

for i=1:length(idx_uniques_VMmin)
    
    disp(strcat('min:',int2str(i)));

    datafile = load(strcat(prefix,int2str(idx_uniques_VMmin(i)),'.mat'));
   
    kVMmin = datafile.kVMmin;
    kVMballmin = datafile.kVMballmin;
    
    VmotionMin(i) = datafile.Wmotion(kVMballmin,kVMmin);
    
    clear datafile
    
end

for i=1:length(idx_uniques_VMmax)
    
    disp(strcat('max:',int2str(i)));

    datafile = load(strcat(prefix,int2str(idx_uniques_VMmax(i)),'.mat'));
   
    kVMmax = datafile.kVMmax;
    kVMballmax = datafile.kVMballmax;
    
    VmotionMax(i) = datafile.Wmotion(kVMballmax,kVMmax);
    
    clear datafile
    
end

if length(idx_uniques_VMmin) > length(idx_uniques_VMmax)
    
    A = idx_uniques_VMmin;
    B = idx_uniques_VMmax;
    
else
    
    B = idx_uniques_VMmin;
    A = idx_uniques_VMmax;
    
end

idx_uniques_VM = A(ismember(A,B));

maxWeight = max(VmotionMax);
minWeight = min(VmotionMin);

nBins = 3;

[X,Y] = hist(round(minWeight):round(maxWeight),nBins);

info.binsCenter(1) = round(minWeight);
info.binsCenter(2) = Y(2);
info.binsCenter(3) = round(maxWeight);
info.minWeight = round(minWeight);
info.maxWeight = round(maxWeight);

nK = 12000;
nBalls = 16;

allVcolor = zeros(length(idx_uniques_VM),nBalls,nK);

for i=1:length(idx_uniques_VM)
    
    disp(strcat('unique:',int2str(i)));

    datafile = load(strcat(prefix,int2str(idx_uniques_VM(i)),'.mat'));
   
    Wmotion = datafile.Wmotion;
    
    allVmotion(i,:,:) = Wmotion(:,:);
    
    clear datafile
    
end

nCenters = length(info.binsCenter);

nSamples = length(idx_uniques_VM);

Time = 10;
FrameRate = 60;
maxNK = 12000;
sFrames = 3;

limitFramesInKs = floor((Time*FrameRate/2)/sFrames);

for nC=1:nCenters
    
   for i=1:length(idx_uniques_VM)
      
       matrix = squeeze(abs(allVmotion(i,:,:) - info.binsCenter(nC)));
       [idx_min_ball, idx_min_kf] = find(matrix == min(min(matrix)));
       
       ki(i) = idx_min_kf(1);
       balli(i) = idx_min_ball(1);
       wi(i) = allVmotion(i,idx_min_ball(1),idx_min_kf(1));
       
   end
    
   nearestValues = abs(wi - info.binsCenter(nC));
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
   
   info.centers(nC).binsMOTidx = idx_uniques_VM(I(1:nSamples));
   info.centers(nC).binsMOTk = ki(I(1:nSamples));
   info.centers(nC).weights = wi(I(1:nSamples)); 
   info.centers(nC).balls = balli(I(1:nSamples));
   
end

save('infoBinsCentersMOTMotionWeight.mat','info'),