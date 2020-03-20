
clear all

velocities = [20 25 30 35 40 45 50];
nMOT = 10;

prefixAnalyzeArea = 'MOT-Analyze-V2-11-m-Area-';

for iVelocity=velocities
    
    uniques = load(strcat('uniqueKEminArea-',int2str(iVelocity),'.mat'));

    idx_unique_areas = uniques.idx_MOT_unique_Area;

    nBins = 3;
    nKs = 13200;

    analyzedAreas = zeros(length(idx_unique_areas),nKs);

    for iMOTAnalyzeArea=1:length(idx_unique_areas)

        dataset = load(strcat(prefixAnalyzeArea,int2str(iVelocity),'-',int2str(idx_unique_areas(iMOTAnalyzeArea)),'.mat'));

        comb(iMOTAnalyzeArea) = dataset.minArea.comb;

        color(iMOTAnalyzeArea) = dataset.minArea.color;

        minareas = squeeze(dataset.areaComb(color(iMOTAnalyzeArea),comb(iMOTAnalyzeArea),:));

        analyzedAreas(iMOTAnalyzeArea,:) = minareas; 

        targets(iMOTAnalyzeArea,:) = dataset.allColorComb(color(iMOTAnalyzeArea)).comb(comb(iMOTAnalyzeArea),:);

        clear dataset   
        clear minareas

    end

    maxArea = max(max(analyzedAreas));
    minArea = min(min(analyzedAreas));

    minInterval = minArea;
    maxInterval = maxArea*0.5;

    [X,Y] = hist(minInterval:maxInterval,nBins);

    info.binsCenter(1) = minInterval;
    info.binsCenter(2) = Y(2);
    info.binsCenter(3) = maxInterval;
    info.minArea = round(minArea);
    info.maxArea = round(maxArea);

    nCenters = length(info.binsCenter);

    nSamples = length(idx_unique_areas);

    Time = 10;
    FrameRate = 60;
    maxNK = 13200;
    sFrames = 3;

    limitFramesInKs = ceil((Time*FrameRate/2)/sFrames);

    for nC=1:nCenters

       for i=1:length(idx_unique_areas)

           line = abs(analyzedAreas(i,:) - info.binsCenter(nC));
           idx_min = find(line == min(line));

           ki(i) = idx_min(1);
           ai(i) = analyzedAreas(i,idx_min(1));

       end

       nearestValues = abs(ai - info.binsCenter(nC));
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

       info.centers(nC).binsMOTidx = idx_unique_areas(I(1:nSamples));
       info.centers(nC).binsMOTk = ki(I(1:nSamples));
       info.centers(nC).area = ai(I(1:nSamples)); 
       info.centers(nC).color = color(I(1:nSamples));
       info.centers(nC).comb = comb((I(1:nSamples)));

       for iSamples=1:nSamples

        info.centers(nC).targets(iSamples,:) = targets((I(iSamples)),:);

       end

    end

    save(strcat('infoBinsCentersMOTArea-',int2str(iVelocity),'.mat'),'info');

end

clear all