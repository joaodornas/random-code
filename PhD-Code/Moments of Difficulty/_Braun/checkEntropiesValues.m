

prefix = 'MOT-Analyze-10-m-';

uniques = load('uniqueKEmin.mat');

idx_uniques_ECmin = uniques.idx_MOT_unique_ECmin;
idx_uniques_ECmax = uniques.idx_MOT_unique_ECmax;

for i=1:length(idx_uniques_ECmin)

    datafile = load(strcat(prefix,int2str(idx_uniques_ECmin(i)),'.mat'));
   
    kECmin = datafile.kECmin;
    
    EcolorMin(i) = datafile.Ecolor(kECmin);
    
    clear datafile
    
end

for i=1:length(idx_uniques_ECmax)

    datafile = load(strcat(prefix,int2str(idx_uniques_ECmax(i)),'.mat'));
   
    kECmax = datafile.kECmax;
    
    EcolorMax(i) = datafile.Ecolor(kECmax);
    
    clear datafile
    
end