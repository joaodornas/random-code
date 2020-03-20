

prefix = 'MOT-Analyze-V2-10-m-';

uniques = load('uniqueKEmin.mat');

idx_uniques_VCmin = uniques.idx_MOT_unique_VCmin;
idx_uniques_VCmax = uniques.idx_MOT_unique_VCmax;

for i=1:length(idx_uniques_VCmin)
    
    disp(strcat('min:',int2str(i)));

    datafile = load(strcat(prefix,int2str(idx_uniques_VCmin(i)),'.mat'));
   
    kVCmin = datafile.kVCmin;
    kVCballmin = datafile.kVCballmin;
    
    VcolorMin(i) = datafile.Wcolor(kVCballmin,kVCmin);
    
    clear datafile
    
end

for i=1:length(idx_uniques_VCmax)
    
    disp(strcat('max:',int2str(i)));

    datafile = load(strcat(prefix,int2str(idx_uniques_VCmax(i)),'.mat'));
   
    kVCmax = datafile.kVCmax;
    kVCballmax = datafile.kVCballmax;
    
    VcolorMax(i) = datafile.Wcolor(kVCballmax,kVCmax);
    
    clear datafile
    
end

if length(idx_uniques_VCmin) > length(idx_uniques_VCmax)
    
    A = idx_uniques_VCmin;
    B = idx_uniques_VCmax;
    
else
    
    B = idx_uniques_VCmin;
    A = idx_uniques_VCmax;
    
end

idx_uniques_VC = A(ismember(A,B));

maxWeight = max(VcolorMax);
minWeight = min(VcolorMin);

nBins = 3;

[X,Y] = hist(round(minWeight):round(maxWeight),nBins);

info.binsCenter(1) = round(minWeight);
info.binsCenter(2) = Y(2);
info.binsCenter(3) = round(maxWeight);
info.minWeight = round(minWeight);
info.maxWeight = round(maxWeight);

nK = 12000;
nBalls = 16;

allVcolor = zeros(length(idx_uniques_VC),nBalls,nK);

for i=1:length(idx_uniques_VC)
    
    disp(strcat('unique:',int2str(i)));

    datafile = load(strcat(prefix,int2str(idx_uniques_VC(i)),'.mat'));
   
    Wcolor = datafile.Wcolor;
    
    allVcolor(i,:,:) = Wcolor(:,:);
    
    clear datafile
    
end

nCenters = length(info.binsCenter);

nSamples = length(idx_uniques_VC);

Time = 10;
FrameRate = 60;
maxNK = 12000;
sFrames = 3;

limitFramesInKs = floor((Time*FrameRate/2)/sFrames);

for nC=1:nCenters
    
   for i=1:length(idx_uniques_VC)
      
       matrix = squeeze(abs(allVcolor(i,:,:) - info.binsCenter(nC)));
       [idx_min_ball, idx_min_kf] = find(matrix == min(min(matrix)));
       
       ki(i) = idx_min_kf(1);
       balli(i) = idx_min_ball(1);
       wi(i) = allVcolor(i,idx_min_ball(1),idx_min_kf(1));
       
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
   
   info.centers(nC).binsMOTidx = idx_uniques_VC(I(1:nSamples));
   info.centers(nC).binsMOTk = ki(I(1:nSamples));
   info.centers(nC).weights = wi(I(1:nSamples)); 
   info.centers(nC).balls = balli(I(1:nSamples));
   
end

save('infoBinsCentersMOTColorWeight.mat','info'),