
start_ = input('Start:');

end_ = input('End:');

prefix = 'MOT-Analyze-10-m-Area-';

for i=start_:end_
    i
    datafile = strcat(prefix,int2str(i),'.mat');
    
    info = load(datafile);
    
    window.minAreaKf(i) = info.minArea.kf;
    
    clear info
    clear datafile
   
end

[uniKf, nuniKf] = count_unique(window.minAreaKf);

repKf = sum(nuniKf>1);

idx_unique_areas = find(nuniKf == 1);
k_unique_Area = uniKf(idx_unique_areas);
idx_MOT_unique_Area = find(ismember(window.minAreaKf,k_unique_Area));

save('uniqueKEminArea.mat');