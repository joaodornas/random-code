
%%% MOMENTS OF DIFFICULTY

%% AREA - MIN
dataset = load('infoBinsCentersMOTArea.mat');

areaMinMOTidx = dataset.info.centers(1).binsMOTidx;

clear dataset

%% COLOR - MIN ENTROPY
dataset = load('infoBinsCentersMOTColorEntropy.mat');

colorMinEntropyMOTidx = dataset.info.centers(1).binsMOTidx;

clear dataset

%% MOTION - MAX SIMILARITY
dataset = load('infoBinsCentersMOTMotionWeight.mat');

motionMaxSimilarityMOTidx = dataset.info.centers(3).binsMOTidx;

clear dataset

%% MOTION - MIN ENTROPY
dataset = load('infoBinsCentersMOTMotionEntropy.mat');

motionMinEntropyMOTidx = dataset.info.centers(1).binsMOTidx;

clear dataset

%%% TEST UNIQUENESS

allMOTidx = [areaMinMOTidx, colorMinEntropyMOTidx, motionMaxSimilarityMOTidx, motionMinEntropyMOTidx];

[uniqueMOT, nuniqueMOT] = count_unique(allMOTidx);
