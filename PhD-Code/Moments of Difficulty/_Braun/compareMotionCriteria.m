

MOTidx = input('MOTidx:');

datasetNormal = load(strcat('MOT-Analyze-10-m-',int2str(MOTidx),'.mat'));

datasetBraun = load(strcat('MOT-Analyze-Braun-10-m-',int2str(MOTidx),'.mat'));

figure

plot(datasetNormal.Emotion,'b');
hold on
plot(datasetBraun.Emotion,'r');

kEMNormalmin = datasetNormal.kEMmin;
kEMNormalmax = datasetNormal.kEMmax;

kEMBraunmin = datasetBraun.kEMmin;
kEMBraunmax = datasetBraun.kEMmax;

clear datasetNormal
clear datasetBraun

