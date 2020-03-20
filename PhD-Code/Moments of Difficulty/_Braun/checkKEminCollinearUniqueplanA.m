
velocities = [20 25 30 35 40 45 50];
nMOT = 10;

prefixanalize = 'MOT-Analyze-V2-11-m-Collinear-';

for iVelocity=velocities
    
    i = 0;
    
    for iMOT=1:nMOT
        
        i = i + 1;
    
        datafile = strcat(prefixanalize,int2str(iVelocity),'-',int2str(iMOT),'.mat');

        info = load(datafile);

        window.minCollinearKf(i) = info.minCollinearity.kf;

        clear info
        clear datafile
    
    end
    
    [uniKf, nuniKf] = count_unique(window.minCollinearKf);

    repKf = sum(nuniKf>1);

    idx_unique_collinearity = find(nuniKf == 1);
    k_unique_collinearity = uniKf(idx_unique_collinearity);
    idx_MOT_unique_collinearity = find(ismember(window.minCollinearKf,k_unique_collinearity));

    save(strcat('uniqueKEminCollinearity-',int2str(iVelocity),'.mat'));
    
    clear uniKf
    clear nuniKf
    clear repKf
    clear idx_unique_collinearity
    clear k_unique_collinearity
    clear idx_MOT_unique_collinearity
    clear window
   
end

clear all