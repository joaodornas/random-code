

prefix = 'MOT-Analyze-V2-11-m-';

velocities = [20 25 30 35 40 45 50];
nMOT = 10;

for iVelocity=velocities
    
    uniques = load(strcat('uniqueKEmin-',int2str(iVelocity),'.mat'));

    idx_uniques_ECmin = uniques.idx_MOT_unique_ECmin;
    idx_uniques_ECmax = uniques.idx_MOT_unique_ECmax;

    for i=1:length(idx_uniques_ECmin)

        disp(strcat('min:',int2str(i)));

        datafile = load(strcat(prefix,int2str(iVelocity),'-',int2str(idx_uniques_ECmin(i)),'.mat'));

        kECmin = datafile.kECmin;

        EcolorMin(i) = datafile.Ecolor(kECmin);

        clear datafile

    end

    for i=1:length(idx_uniques_ECmax)

        disp(strcat('max:',int2str(i)));

        datafile = load(strcat(prefix,int2str(iVelocity),'-',int2str(idx_uniques_ECmax(i)),'.mat'));

        kECmax = datafile.kECmax;

        EcolorMax(i) = datafile.Ecolor(kECmax);

        clear datafile

    end

    if length(idx_uniques_ECmin) > length(idx_uniques_ECmax)

        A = idx_uniques_ECmin;
        B = idx_uniques_ECmax;

    else

        B = idx_uniques_ECmin;
        A = idx_uniques_ECmax;

    end

    idx_uniques_EC = A(ismember(A,B));

    maxEntropy = max(EcolorMax);
    minEntropy = min(EcolorMin);

    nBins = 3;

    [X,Y] = hist(round(minEntropy):round(maxEntropy),nBins);

    info.binsCenter(1) = round(minEntropy);
    info.binsCenter(2) = Y(2);
    info.binsCenter(3) = round(maxEntropy);
    info.minEntropy = round(minEntropy);
    info.maxEntropy = round(maxEntropy);

    nK = 13200;

    allEcolor = zeros(length(idx_uniques_EC),nK);

    for i=1:length(idx_uniques_EC)

        disp(strcat('unique:',int2str(i)));

        datafile = load(strcat(prefix,int2str(iVelocity),'-',int2str(idx_uniques_EC(i)),'.mat'));

        Ecolor = datafile.Ecolor;

        allEcolor(i,:) = Ecolor;

        clear datafile

    end

    nCenters = length(info.binsCenter);

    nSamples = length(idx_uniques_EC);

    Time = 10;
    FrameRate = 60;
    maxNK = 13200;
    sFrames = 3;

    limitFramesInKs = floor((Time*FrameRate/2)/sFrames);

    for nC=1:nCenters

       for i=1:length(idx_uniques_EC)

           line = abs(allEcolor(i,:) - info.binsCenter(nC));
           idx_min = find(line == min(line));

           ki(i) = idx_min(1);
           ei(i) = allEcolor(i,idx_min(1));

       end

       nearestValues = abs(ei - info.binsCenter(nC));
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

       info.centers(nC).binsMOTidx = idx_uniques_EC(I(1:nSamples));
       info.centers(nC).binsMOTk = ki(I(1:nSamples));
       info.centers(nC).entropies = ei(I(1:nSamples)); 

    end

    save(strcat('infoBinsCentersMOTColorEntropy-',int2str(iVelocity),'.mat'),'info'),

end