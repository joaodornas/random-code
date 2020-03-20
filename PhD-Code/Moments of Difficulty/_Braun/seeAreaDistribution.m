
clear all

prefix = 'MOT-Analyze-10-m-Area-';

nAreas = 500;

for i=1:nAreas
   
    dataset = load(strcat(prefix,int2str(i),'.mat'));
    
    area(i) = dataset.minArea.area;
    
    clear dataset
    
end

[unique nunique] = count_unique(area);

idx_uniques_ones = find(nunique == 1);
unique_areas = unique(idx_uniques_ones);

idx_unique_areas = find(ismember(area,unique_areas));

save('areaInfo.mat','area','idx_unique_areas');
