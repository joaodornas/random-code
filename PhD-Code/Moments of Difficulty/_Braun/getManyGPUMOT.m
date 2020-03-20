
clear all

ItemsN = 16;

gpuItemsN = gpuArray(ItemsN);

gpuMOT = arrayfun(@GenerateBraunianWalkTrajectories,gpuItemsN);

MOT = gather(gpuMOT);