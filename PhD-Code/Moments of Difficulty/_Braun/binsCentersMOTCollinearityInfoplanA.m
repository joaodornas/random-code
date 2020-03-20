
clear all

velocities = [20 25 30 35 40 45 50];
nMOT = 10;

prefixAnalyzeCollinearity = 'MOT-Analyze-V2-11-m-Collinear-';

for iVelocity=velocities
    
    uniques = load(strcat('uniqueKEminCollinearity-',int2str(iVelocity),'.mat'));

    idx_unique_collinearity = uniques.idx_MOT_unique_collinearity;

    nBins = 3;
    nKs = 13200;

    analyzedCollinearity = zeros(length(idx_unique_collinearity),nKs);

    for iMOTAnalyzeCollinearity=1:length(idx_unique_collinearity)

        dataset = load(strcat(prefixAnalyzeCollinearity,int2str(iVelocity),'-',int2str(idx_unique_collinearity(iMOTAnalyzeCollinearity)),'.mat'));

        comb(iMOTAnalyzeCollinearity) = dataset.minCollinearity.comb;

        color(iMOTAnalyzeCollinearity) = dataset.minCollinearity.color;

        mincollinearity = squeeze(dataset.collinearComb(color(iMOTAnalyzeCollinearity),comb(iMOTAnalyzeCollinearity),:));

        analyzedCollinearity(iMOTAnalyzeCollinearity,:) = mincollinearity; 

        targets(iMOTAnalyzeCollinearity,:) = dataset.allColorComb(color(iMOTAnalyzeCollinearity)).comb(comb(iMOTAnalyzeCollinearity),:);

        clear dataset   
        clear minareas

    end

    maxCollinearity = max(max(analyzedCollinearity));
    minCollinearity = min(min(analyzedCollinearity));

    minInterval = minCollinearity;
    %maxInterval = maxCollinearity*0.01;
    maxInterval = 0.1;

    [X,Y] = hist(minInterval:maxInterval,nBins);

    info.binsCenter(1) = minInterval;
    info.binsCenter(2) = Y(2);
    info.binsCenter(3) = maxInterval;
    info.minCollinear = round(minCollinearity);
    info.maxCollinear = round(maxCollinearity);

    nCenters = length(info.binsCenter);

    nSamples = length(idx_unique_collinearity);

    Time = 10;
    FrameRate = 60;
    maxNK = 13200;
    sFrames = 3;

    limitFramesInKs = ceil((Time*FrameRate/2)/sFrames);

    for nC=1:nCenters

       for i=1:length(idx_unique_collinearity)

           line = abs(analyzedCollinearity(i,:) - info.binsCenter(nC));
           idx_min = find(line == min(line));

           ki(i) = idx_min(1);
           ai(i) = analyzedCollinearity(i,idx_min(1));

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

       info.centers(nC).binsMOTidx = idx_unique_collinearity(I(1:nSamples));
       info.centers(nC).binsMOTk = ki(I(1:nSamples));
       info.centers(nC).collinearity = ai(I(1:nSamples)); 
       info.centers(nC).color = color(I(1:nSamples));
       info.centers(nC).comb = comb((I(1:nSamples)));

       for iSamples=1:nSamples

        info.centers(nC).targets(iSamples,:) = targets((I(iSamples)),:);

       end

    end

    save(strcat('infoBinsCentersMOTCollinearity-',int2str(iVelocity),'.mat'),'info');

end

clear all