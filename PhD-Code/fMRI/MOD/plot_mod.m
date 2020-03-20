
load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\MOD\stimuli\MOT-14-1-Analyze.mat');
load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\MOD\stimuli\MOT-14-1_5.5_min.mat');
idx_k_MOD_Color_mv1 = [146, 1100, 2135, 2385, 2896, 3220, 3764, 4511, 5426, 5974];

for kemin=idx_k_MOD_Color_mv1
   
    figure_label = strcat('MOD-Color-mv1-k-',int2str(kemin),'.eps');
    plotWindowColorAuto;
    
end

clear all

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\MOD\stimuli\MOT-14-2-Analyze.mat');
load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\MOD\stimuli\MOT-14-2_5.5_min.mat');
idx_k_MOD_Color_mv2 = [106, 1170, 2279, 2472, 3992, 4988, 5386, 5608, 5940, 6589];

for kemin=idx_k_MOD_Color_mv2
   
    figure_label = strcat('MOD-Color-mv2-k-',int2str(kemin),'.eps');
    plotWindowColorAuto;
    
end

clear all