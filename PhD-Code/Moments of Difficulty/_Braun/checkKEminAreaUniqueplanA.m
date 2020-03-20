
velocities = [20 25 30 35 40 45 50];
nMOT = 10;

prefixanalize = 'MOT-Analyze-V2-11-m-Area-';

for iVelocity=velocities
    
    i = 0;
    
    for iMOT=1:nMOT
        
        i = i + 1;
    
        datafile = strcat(prefixanalize,int2str(iVelocity),'-',int2str(iMOT),'.mat');

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

    save(strcat('uniqueKEminArea-',int2str(iVelocity),'.mat'));
    
    clear uniKf
    clear nuniKf
    clear repKf
    clear idx_unique_areas
    clear k_unique_Area
    clear idx_MOT_unique_Area
    clear window
   
end

clear all