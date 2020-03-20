function real_FC_voxel_AAL_ROI_densities_stat_voxels


%computeEverything;
plotCircularRelationship_v4;



end

function computeEverything

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

[idx_Att_Neg, idx_Att_Pos, idx_Pass_Neg, idx_Pass_Pos] = getIDXVoxelsDensities;

[AAL_ROIs, AAL_Densities] = getAmountVoxelsInAALROIPerDensity(idx_Att_Neg, idx_Att_Pos, idx_Pass_Neg, idx_Pass_Pos,AAL_img,AAL_ROI);

AAL_StatsResults = statistical_inference(AAL_Densities);

Nets_Densities = getAmountVoxelsInFunctionalNetworksPerDensity(idx_Att_Neg, idx_Att_Pos, idx_Pass_Neg, idx_Pass_Pos,AAL_img,AAL_ROI);

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

for iNet=1:length(all_networks_labels)
    
    Nets_Densities.net(iNet).StatsResults = statistical_inference(Nets_Densities.net(iNet).Densities);
    
end

save('Non-Parametric-Voxel-Wise-Per-Density-AAL-FuncNets-AAL-IDs.mat');

end

function [idx_Att_Neg, idx_Att_Pos, idx_Pass_Neg, idx_Pass_Pos] = getIDXVoxelsDensities

zscore = 1.6;

folder = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\densities\rendering\both\clusters-dlh=1.5-both';

filename_Att_Neg = 'LHR-Correlation-Contrast-Track-Neg-z-clu-both';
filename_Att_Pos = 'LHR-Correlation-Contrast-Track-Pos-z-clu-both';

filename_Pass_Neg = 'LHR-Correlation-Contrast-Passive-Neg-z-clu-both';
filename_Pass_Pos = 'LHR-Correlation-Contrast-Passive-Pos-z-clu-both';

Att_Neg_img = loadNifTiFile(folder,filename_Att_Neg);
Att_Pos_img = loadNifTiFile(folder,filename_Att_Pos);

Pass_Neg_img = loadNifTiFile(folder,filename_Pass_Neg);
Pass_Pos_img = loadNifTiFile(folder,filename_Pass_Pos);

idx_Att_Neg.increase = find(Att_Neg_img > zscore);
idx_Att_Neg.decrease = find(Att_Neg_img < -zscore);

idx_Att_Pos.increase = find(Att_Pos_img > zscore);
idx_Att_Pos.decrease = find(Att_Pos_img < -zscore);

idx_Pass_Neg.increase = find(Pass_Neg_img > zscore);
idx_Pass_Neg.decrease = find(Pass_Neg_img < -zscore);

idx_Pass_Pos.increase = find(Pass_Pos_img > zscore);
idx_Pass_Pos.decrease = find(Pass_Pos_img < -zscore);

end

function img = loadNifTiFile(folder,filename)

path(path,folder);

load_img = nifti(strcat(filename,'.nii'));
load_img.dat.fname = strcat(folder,'\',load_img.dat.fname);

img = load_img.dat(:,:,:);

end

function [AAL_ROIs, Densities] = getAmountVoxelsInAALROIPerDensity(idx_Att_Neg, idx_Att_Pos, idx_Pass_Neg, idx_Pass_Pos,AAL_img,AAL_ROI)

% %%% LOAD AAL
% 
% load_aal = nifti('ROI_MNI_V4.nii');
% load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
% 
% AAL_img = load_aal.dat(:,:,:);
% 
% load_roi = load('ROI_MNI_V4_List.mat');
% 
% AAL_ROI = load_roi.ROI;

nNodes = 90;
nDensities = 8;

AAL_ROIs = cell(nNodes+1,nDensities+1);
Densities = zeros(nNodes,nDensities);

AAL_ROIs{1,1} = 'AAL ROI';

for iROI=1:nNodes
    
    AAL_ROIs{iROI+1,1} = AAL_ROI(iROI).Nom_L;
    
end

AAL_ROIs{1,2} = 'Pass-Neg-Down';
AAL_ROIs{1,3} = 'Pass-Neg-Up';
AAL_ROIs{1,4} = 'Pass-Pos-Down';
AAL_ROIs{1,5} = 'Pass-Pos-Up';
AAL_ROIs{1,6} = 'Att-Neg-Down';
AAL_ROIs{1,7} = 'Att-Neg-Up';
AAL_ROIs{1,8} = 'Att-Pos-Down';
AAL_ROIs{1,9} = 'Att-Pos-Up';

for iROI=1:nNodes
    
    idx_ROI = AAL_ROI(iROI).ID;
    idx_voxels = find(AAL_img==idx_ROI);
    
    AAL_ROIs{iROI+1,2} = length(find(ismember(idx_voxels,idx_Pass_Neg.decrease)));
    Densities(iROI,1) = length(find(ismember(idx_voxels,idx_Pass_Neg.decrease)));
    
    AAL_ROIs{iROI+1,3} = length(find(ismember(idx_voxels,idx_Pass_Neg.increase)));
    Densities(iROI,2) = length(find(ismember(idx_voxels,idx_Pass_Neg.increase)));
    
    AAL_ROIs{iROI+1,4} = length(find(ismember(idx_voxels,idx_Pass_Pos.decrease)));
    Densities(iROI,3) = length(find(ismember(idx_voxels,idx_Pass_Pos.decrease)));
    
    AAL_ROIs{iROI+1,5} = length(find(ismember(idx_voxels,idx_Pass_Pos.increase)));
    Densities(iROI,4) = length(find(ismember(idx_voxels,idx_Pass_Pos.increase)));
    
    AAL_ROIs{iROI+1,6} = length(find(ismember(idx_voxels,idx_Att_Neg.decrease)));
    Densities(iROI,5) = length(find(ismember(idx_voxels,idx_Att_Neg.decrease)));
    
    AAL_ROIs{iROI+1,7} = length(find(ismember(idx_voxels,idx_Att_Neg.increase)));
    Densities(iROI,6) = length(find(ismember(idx_voxels,idx_Att_Neg.increase)));
   
    AAL_ROIs{iROI+1,8} = length(find(ismember(idx_voxels,idx_Att_Pos.decrease)));
    Densities(iROI,7) = length(find(ismember(idx_voxels,idx_Att_Pos.decrease))); 
    
    AAL_ROIs{iROI+1,9} = length(find(ismember(idx_voxels,idx_Att_Pos.increase)));
    Densities(iROI,8) = length(find(ismember(idx_voxels,idx_Att_Pos.increase)));

end

end

function Nets_Densities = getAmountVoxelsInFunctionalNetworksPerDensity(idx_Att_Neg, idx_Att_Pos, idx_Pass_Neg, idx_Pass_Pos,AAL_img,AAL_ROI)

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

nDensities = 8;

for iNet=1:length(all_networks_labels)
    
    network_label = all_networks_labels{iNet};
    
    seeds = getFunctionalSeeds_v5(network_label);
    
    nSeeds = length(seeds.ROI);

    Net_Seeds = cell(nSeeds+1,nDensities+1);
    Densities = zeros(nSeeds,nDensities);
   
    Net_Seeds{1,1} = 'Net Seed';

    for iSeed=1:nSeeds
    
        Net_Seeds{iSeed+1,1} = seeds.ROI(iSeed).label;
    
    end
    
    Net_Seeds{1,2} = 'Pass-Neg-Down';
    Net_Seeds{1,3} = 'Pass-Neg-Down';
    Net_Seeds{1,4} = 'Pass-Neg-Up';
    Net_Seeds{1,5} = 'Pass-Neg-Up';
    Net_Seeds{1,6} = 'Pass-Pos-Down';
    Net_Seeds{1,7} = 'Pass-Pos-Down';
    Net_Seeds{1,8} = 'Pass-Pos-Up';
    Net_Seeds{1,9} = 'Pass-Pos-Up';
    Net_Seeds{1,10} = 'Att-Neg-Down';
    Net_Seeds{1,11} = 'Att-Neg-Down';
    Net_Seeds{1,12} = 'Att-Neg-Up';
    Net_Seeds{1,13} = 'Att-Neg-Up';
    Net_Seeds{1,14} = 'Att-Pos-Down';
    Net_Seeds{1,15} = 'Att-Pos-Down';
    Net_Seeds{1,16} = 'Att-Pos-Up';
    Net_Seeds{1,17} = 'Att-Pos-Up';
    
    folder = 'Z:\Dropbox (Uni Magdeburg)\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually';
    path(path,folder);
    
    load_img = nifti(strcat('LHR-All-Subjects-Functional-Parcels-Pop-Map-',network_label,'-seed.nii'));
    load_img.dat.fname = strcat(folder,'\',load_img.dat.fname);

    Net_img = load_img.dat(:,:,:);
    
    for iSeed=1:nSeeds
    
        idx_seed = seeds.ROI(iSeed).idx;
        idx_voxels = find(Net_img==idx_seed);

        Net_Seeds{iSeed+1,2} = length(find(ismember(idx_voxels,idx_Pass_Neg.decrease)));
        Densities(iSeed,1) = length(find(ismember(idx_voxels,idx_Pass_Neg.decrease)));
        Net_Seeds{iSeed+1,3} = getROIsFromIDXVoxels(idx_voxels(ismember(idx_voxels(:),idx_Pass_Neg.decrease(:))),AAL_img,AAL_ROI);

        Net_Seeds{iSeed+1,4} = length(find(ismember(idx_voxels,idx_Pass_Neg.increase)));
        Densities(iSeed,2) = length(find(ismember(idx_voxels,idx_Pass_Neg.increase)));
        Net_Seeds{iSeed+1,5} = getROIsFromIDXVoxels(idx_voxels(ismember(idx_voxels(:),idx_Pass_Neg.increase(:))),AAL_img,AAL_ROI);

        Net_Seeds{iSeed+1,6} = length(find(ismember(idx_voxels,idx_Pass_Pos.decrease)));
        Densities(iSeed,3) = length(find(ismember(idx_voxels,idx_Pass_Pos.decrease)));
        Net_Seeds{iSeed+1,7} = getROIsFromIDXVoxels(idx_voxels(ismember(idx_voxels(:),idx_Pass_Pos.decrease(:))),AAL_img,AAL_ROI);

        Net_Seeds{iSeed+1,8} = length(find(ismember(idx_voxels,idx_Pass_Pos.increase)));
        Densities(iSeed,4) = length(find(ismember(idx_voxels,idx_Pass_Pos.increase)));
        Net_Seeds{iSeed+1,9} = getROIsFromIDXVoxels(idx_voxels(ismember(idx_voxels(:),idx_Pass_Pos.increase(:))),AAL_img,AAL_ROI);

        Net_Seeds{iSeed+1,10} = length(find(ismember(idx_voxels,idx_Att_Neg.decrease)));
        Densities(iSeed,5) = length(find(ismember(idx_voxels,idx_Att_Neg.decrease)));
        Net_Seeds{iSeed+1,11} = getROIsFromIDXVoxels(idx_voxels(ismember(idx_voxels(:),idx_Att_Neg.decrease(:))),AAL_img,AAL_ROI);

        Net_Seeds{iSeed+1,12} = length(find(ismember(idx_voxels,idx_Att_Neg.increase)));
        Densities(iSeed,6) = length(find(ismember(idx_voxels,idx_Att_Neg.increase)));
        Net_Seeds{iSeed+1,13} = getROIsFromIDXVoxels(idx_voxels(ismember(idx_voxels(:),idx_Att_Neg.increase(:))),AAL_img,AAL_ROI);

        Net_Seeds{iSeed+1,14} = length(find(ismember(idx_voxels,idx_Att_Pos.decrease)));
        Densities(iSeed,7) = length(find(ismember(idx_voxels,idx_Att_Pos.decrease))); 
        Net_Seeds{iSeed+1,15} = getROIsFromIDXVoxels(idx_voxels(ismember(idx_voxels(:),idx_Att_Pos.decrease(:))),AAL_img,AAL_ROI);

        Net_Seeds{iSeed+1,16} = length(find(ismember(idx_voxels,idx_Att_Pos.increase)));
        Densities(iSeed,8) = length(find(ismember(idx_voxels,idx_Att_Pos.increase)));
        Net_Seeds{iSeed+1,17} = getROIsFromIDXVoxels(idx_voxels(ismember(idx_voxels(:),idx_Att_Pos.increase(:))),AAL_img,AAL_ROI);

    end
    
    Nets_Densities.net(iNet).Net_Seeds = Net_Seeds;
    Nets_Densities.net(iNet).Densities = Densities;
    Nets_Densities.net(iNet).network_label = network_label;

end


end

function StatsResults = statistical_inference(Densities)

StatsResults{1,1} = 'DENSITY';
StatsResults{1,2} = 'H';
StatsResults{1,3} = 'P';
StatsResults{1,4} = 'Diff-Mean';
StatsResults{1,5} = 'Mean-Voxels-Attention';
StatsResults{1,6} = 'Mean-Voxels-Passive';

x = Densities(:,1);
y = Densities(:,5);
[P,H] = ranksum(x,y);

StatsResults{2,1} = 'Neg-Down';
StatsResults{2,2} = H;
StatsResults{2,3} = P;
StatsResults{2,4} = mean(y)-mean(x);
StatsResults{2,5} = mean(y);
StatsResults{2,6} = mean(x);

x = Densities(:,2);
y = Densities(:,6);
[P,H] = ranksum(x,y);

StatsResults{3,1} = 'Neg-Up';
StatsResults{3,2} = H;
StatsResults{3,3} = P;
StatsResults{3,4} = mean(y)-mean(x);
StatsResults{3,5} = mean(y);
StatsResults{3,6} = mean(x);

x = Densities(:,3);
y = Densities(:,7);
[P,H] = ranksum(x,y);

StatsResults{4,1} = 'Pos-Down';
StatsResults{4,2} = H;
StatsResults{4,3} = P;
StatsResults{4,4} = mean(y)-mean(x);
StatsResults{4,5} = mean(y);
StatsResults{4,6} = mean(x);

x = Densities(:,4);
y = Densities(:,8);
[P,H] = ranksum(x,y);

StatsResults{5,1} = 'Pos-Up';
StatsResults{5,2} = H;
StatsResults{5,3} = P;
StatsResults{5,4} = mean(y)-mean(x);
StatsResults{5,5} = mean(y);
StatsResults{5,6} = mean(x);

end

function ROI_Names_Amounts = getROIsFromIDXVoxels(idx_voxels,AAL_img,AAL_ROI)

% %%% LOAD AAL
% 
% load_aal = nifti('ROI_MNI_V4.nii');
% load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
% 
% AAL_img = load_aal.dat(:,:,:);
% 
% load_roi = load('ROI_MNI_V4_List.mat');
% 
% AAL_ROI = load_roi.ROI;

nNodes = 90;
iiROI = 0;

ROI_Names_Amounts = cell(2,2);

for iROI=1:nNodes
    
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels_AAL = find(AAL_img==idx_ROI);
    
    nVoxelsInside = length(find(ismember(idx_voxels,idx_voxels_AAL)));
    
    if nVoxelsInside > 0
   
        iiROI = iiROI + 1;
        
        ROI_Names_Amounts{iiROI+1,1} = label_ROI{iROI};
        ROI_Names_Amounts{iiROI+1,2} = nVoxelsInside;
        
    end
    
end
   
end

function plotCircularRelationship

load('Non-Parametric-Voxel-Wise-Per-Density-AAL-FuncNets-AAL-IDs.mat');

%%% colors

color(1,:) = [0,0,128];
color(2,:) = [128,0,230];
color(3,:) = [0,191,255];
color(4,:) = [0,100,0];
color(5,:) = [255,255,0];
color(6,:) = [232,105,43];
color(7,:) = [255,0,0];
color(8,:) = [0,0,0];

color = color./255;

%%% switch colors AUD for DMN

AUD_color = color(7,:);
DMN_color = color(8,:);

color(7,:) = DMN_color(:);
color(8,:) = AUD_color(:);

%%% LOAD NETS

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

nAAL_ROI = 90;
nLobesAndSubCortical = 5;
nNet = 8;

% rgb = round(linspace(0,255,nLobesAndSubCortical+1))./255;
% rgb = rgb(1:end-1);
rgb(1) = 192;
rgb(2) = 128;
rgb(3) = 32;
rgb(4) = 224;
rgb(5) = 160;
rgb = rgb./255;

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

for iROI=1:nAAL_ROI
    
    roi_label{iROI} = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    
    if ismember(iROI,idx_frontal)
        
        aal_color = 1;
        
    elseif ismember(iROI,idx_occipital)
        
        aal_color = 2;
        
    elseif ismember(iROI,idx_parietal)
        
        aal_color = 3;
        
    elseif ismember(iROI,idx_temporal)
        
        aal_color = 4;
        
    elseif ismember(iROI,idx_subcortical)
        
        aal_color = 5;
        
    end
    
    roi_color(iROI,:) = [rgb(aal_color), rgb(aal_color), rgb(aal_color)];
    
end

roi_label = roi_label([idx_frontal,idx_occipital,idx_parietal,idx_temporal,idx_subcortical]);
roi_color = roi_color([idx_frontal,idx_occipital,idx_parietal,idx_temporal,idx_subcortical],:);

iiROI = 0;
for iNet=1:nNet
    
    seeds = getFunctionalSeeds_v5(all_networks_labels{iNet});
    
    nROI = length(seeds.ROI);
    Net(iNet).nSeeds = nROI;
    Net(iNet).network_label = all_networks_labels{iNet};
    Net(iNet).color(1:3) = color(iNet,:);
    
    for iROI=1:nROI
       
        iiROI = iiROI + 1;
        
        seed_label{iiROI} = seeds.ROI(iROI).label;
        net_label{iiROI} = Net(iNet).network_label;
        net_color(iiROI,:) = Net(iNet).color(:);
        
    end
    
end

nodes_labels = cellstr({roi_label{1:end},seed_label{1:end}});
nodes_colors = [roi_color(:,1:end);net_color(:,1:end)];

nSeeds = length(seed_label);

r = 1;
theta=linspace(0,2*pi,nAAL_ROI+nSeeds+1);
theta=theta(1:end-1);
x = r .* cos(theta);
y = r .* sin(theta);

[pnu,pnd,ppu,ppd,anu,and,apu,apd] =  getABinaryMatrix(Nets_Densities,nSeeds,nAAL_ROI,AAL_img,AAL_ROI);
%[neg_up, neg_down, pos_up, pos_down, xy] = getAxy(Nets_Densities,pnu,pnd,ppu,ppd,anu,and,apu,apd,x,y,nAAL_ROI,nSeeds,Net);

for i=1:length(x)
    
      xy(i,:) = [x(i) y(i)];
%     xy(i,1) = x(i);
%     xy(i,2) = y(i);
    
end

pcriterion = 0.01;

%%% REMOVE NON-SIGNIFICANT BASED ON WILCOXON

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

nNet = 8;

neg_up = anu;
neg_down = and;
pos_up = apu;
pos_down = apd;

for iNet=1:nNet
    
    nNet_Seeds = Net(iNet).nSeeds;
    
    last_seeds = 0;
    for iiNet=1:(iNet-1)
        
        last_seeds = size(Nets_Densities.net(iiNet).Net_Seeds,1) - 1 + last_seeds;
        
    end
    
    if Nets_Densities.net(iNet).StatsResults{2,3} > pcriterion
        
        neg_down((nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds),:) = 0;
        neg_down(:,(nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds)) = 0;
        
    end
        
    if Nets_Densities.net(iNet).StatsResults{3,3} > pcriterion
        
        neg_up((nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds),:) = 0;
        neg_up(:,(nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds)) = 0;
        
    end
        
    if Nets_Densities.net(iNet).StatsResults{4,3} > pcriterion
        
        pos_down((nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds),:) = 0;
        pos_down(:,(nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds)) = 0;
        
    end
    
    if Nets_Densities.net(iNet).StatsResults{5,3} > pcriterion
        
        pos_up((nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds),:) = 0;
        pos_up(:,(nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds)) = 0;
        
    end
  
end

%%% show nodes and edges
f = figure;

for i=1:length(x)-1
    line([x(i) x(i+1)], [y(i) y(i+1)], 'LineWidth', 15, 'Color',nodes_colors(i,:));
    hold on
end

%gplot(neg_up,xy,'b-');

for i=1:(length(x)-nAAL_ROI)

    for j=(nAAL_ROI+1):length(x)
       
        if neg_up(i,j) == 1
        
            plot([x(i) x(j)],[y(i) y(j)],'LineWidth',3,'Color',net_color(j-nAAL_ROI,:));
            hold on
            
        end
        
    end
    
end

for i=2:length(x)
    h = text(x(i).*1.05,y(i).*1.05,nodes_labels{i},'FontSize',8);
    set(h,'Rotation',(theta(i-1)*180/pi));
end

title('Negative UP');

g = figure;

for i=1:length(x)-1
    line([x(i) x(i+1)], [y(i) y(i+1)], 'LineWidth', 15, 'Color',nodes_colors(i,:));
    hold on
end

%gplot(pos_down,xy,'b-');

for i=1:(length(x)-nAAL_ROI)

    for j=(nAAL_ROI+1):length(x)
       
        if pos_down(i,j) == 1
        
            plot([x(i) x(j)],[y(i) y(j)],'LineWidth',3,'Color',net_color(j-nAAL_ROI,:));
            
        end
        
    end
    
end

for i=2:length(x)
    h = text(x(i).*1.05,y(i).*1.05,nodes_labels{i},'FontSize',8);
    set(h,'Rotation',(theta(i-1)*180/pi));
end

title('Positive Down');

end

function plotCircularRelationship_v2

load('Non-Parametric-Voxel-Wise-Per-Density-AAL-FuncNets-AAL-IDs-grouped-seeds.mat');

%%% colors

color(1,:) = [0,0,128];
color(2,:) = [128,0,230];
color(3,:) = [0,191,255];
color(4,:) = [0,100,0];
color(5,:) = [255,255,0];
color(6,:) = [232,105,43];
color(7,:) = [255,0,0];
color(8,:) = [0,0,0];

color = color./255;

%%% switch colors AUD for DMN

AUD_color = color(7,:);
DMN_color = color(8,:);

color(7,:) = DMN_color(:);
color(8,:) = AUD_color(:);

%%% LOAD NETS

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

nAAL_ROI = 90;
nLobesAndSubCortical = 5;
nNet = 8;

% rgb = round(linspace(0,255,nLobesAndSubCortical+1))./255;
% rgb = rgb(1:end-1);
rgb(1) = 192;
rgb(2) = 128;
rgb(3) = 32;
rgb(4) = 224;
rgb(5) = 160;
rgb = rgb./255;

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

for iROI=1:nAAL_ROI
    
    roi_label{iROI} = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    
    if ismember(iROI,idx_frontal)
        
        aal_color = 1;
        
    elseif ismember(iROI,idx_occipital)
        
        aal_color = 2;
        
    elseif ismember(iROI,idx_parietal)
        
        aal_color = 3;
        
    elseif ismember(iROI,idx_temporal)
        
        aal_color = 4;
        
    elseif ismember(iROI,idx_subcortical)
        
        aal_color = 5;
        
    end
    
    roi_color(iROI,:) = [rgb(aal_color), rgb(aal_color), rgb(aal_color)];
    
end

roi_label = roi_label([idx_frontal,idx_occipital,idx_parietal,idx_temporal,idx_subcortical]);
roi_color = roi_color([idx_frontal,idx_occipital,idx_parietal,idx_temporal,idx_subcortical],:);

iiROI = 0;
for iNet=1:nNet
    
    seeds = getFunctionalSeeds_v5(all_networks_labels{iNet});
    
    nROI = length(seeds.ROI);
    Net(iNet).nSeeds = nROI;
    Net(iNet).network_label = all_networks_labels{iNet};
    Net(iNet).color(1:3) = color(iNet,:);
    
    for iROI=1:nROI
       
        iiROI = iiROI + 1;
        
        seed_label{iiROI} = seeds.ROI(iROI).label;
        net_label{iiROI} = Net(iNet).network_label;
        net_color(iiROI,:) = Net(iNet).color(:);
        
    end
    
end

nodes_labels = cellstr({roi_label{1:end},seed_label{1:end}});
nodes_colors = [roi_color(:,1:end);net_color(:,1:end)];

nSeeds = length(seed_label);

% r = 1;
% theta=linspace(0,2*pi,nAAL_ROI+nSeeds+1);
% theta=theta(1:end-1);
% x = r .* cos(theta);
% y = r .* sin(theta);

r = 1;
theta=linspace(0,2*pi,nAAL_ROI+1);
theta=theta(1:end-1);
x = r .* cos(theta);
y = r .* sin(theta);

% for i=1:length(x)
%     
%       xy(i,:) = [x(i) y(i)];
% %     xy(i,1) = x(i);
% %     xy(i,2) = y(i);
%     
% end

%%% show nodes and edges
f = figure;

biggest_amount_voxels = 0;

for iNet=1:iNet
    
    nSeeds = size(grouped_seeds_Nets_Densities.net(iNet).Net_Seeds,1) - 1;
    
    iSeedAAL = 0;
    
    for iSeed=1:nSeeds
        
        aal_names = grouped_seeds_Nets_Densities.net(iNet).Net_Seeds{iSeed+1,13};
        
        for iaal=1:size(aal_names,1)
            
            if ~isempty(aal_names(iaal,1))
                
                roi_per_seed_label = strrep(aal_names(iaal,1),'_','-');
                nVoxels_roi_per_seed = aal_names(iaal,2);
                nVoxels = nVoxels_roi_per_seed{1};
                
                if biggest_amount_voxels < nVoxels; biggest_amount_voxels = nVoxels; end
                
                for iiaal=1:length(roi_label)
                    
                    if strcmp(roi_per_seed_label,roi_label(iiaal));
                        
                        iSeedAAL = iSeedAAL + 1;
                        
                        net_seed_aal_x(iSeedAAL) = x(iiaal) .* nVoxels_roi_per_seed{1};
                        net_seed_aal_y(iSeedAAL) = y(iiaal) .* nVoxels_roi_per_seed{1};
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
    h = compass(net_seed_aal_x,net_seed_aal_y); 
    set(h,'color', color(iNet,:))
    
    ylim([-biggest_amount_voxels, biggest_amount_voxels]);
    xlim([-biggest_amount_voxels, biggest_amount_voxels]);
    
    hold on
    
    clear net_seed_aal_x
    clear net_seed_aal_y;
    
end

x = x .* biggest_amount_voxels;
y = y .* biggest_amount_voxels;

for i=1:length(x)-1
    line([x(i) x(i+1)], [y(i) y(i+1)], 'LineWidth', 15, 'Color',roi_color(i,:));
    hold on
end

for i=1:length(x)
    h = text(x(i).*1.05,y(i).*1.05,roi_label{i},'FontSize',8);
    set(h,'Rotation',(theta(i)*180/pi));
end
        

end

function plotCircularRelationship_v3

load('Non-Parametric-Voxel-Wise-Per-Density-AAL-FuncNets-AAL-IDs-grouped-seeds.mat');

%%% colors

color(1,:) = [0,0,128];
color(2,:) = [128,0,230];
color(3,:) = [0,191,255];
color(4,:) = [0,100,0];
color(5,:) = [255,255,0];
color(6,:) = [232,105,43];
color(8,:) = [255,0,0];
color(7,:) = [0,0,0];

color_label{1} = 'b';
color_label{2} = 'm';
color_label{3} = 'c';
color_label{4} = 'g';
color_label{5} = 'y';
color_label{6} = 'g';
color_label{8} = 'b';
color_label{7} = 'r';

color = color./255;

%%% LOAD NETS

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

nAAL_ROI = 90;
nNet = 8;

% rgb = round(linspace(0,255,nLobesAndSubCortical+1))./255;
% rgb = rgb(1:end-1);
rgb(1) = 192;
rgb(2) = 128;
rgb(3) = 32;
rgb(4) = 224;
rgb(5) = 160;
rgb = rgb./255;

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

for iROI=1:nAAL_ROI
    
    roi_label{iROI} = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    
    if ismember(iROI,idx_frontal)
        
        aal_color = 1;
        
    elseif ismember(iROI,idx_occipital)
        
        aal_color = 2;
        
    elseif ismember(iROI,idx_parietal)
        
        aal_color = 3;
        
    elseif ismember(iROI,idx_temporal)
        
        aal_color = 4;
        
    elseif ismember(iROI,idx_subcortical)
        
        aal_color = 5;
        
    end
    
    roi_color(iROI,:) = [rgb(aal_color), rgb(aal_color), rgb(aal_color)];
    
end

%roi_label = roi_label{[idx_frontal,idx_occipital,idx_parietal,idx_temporal,idx_subcortical]};
new_aal_order = [idx_frontal,idx_occipital,idx_parietal,idx_temporal,idx_subcortical];
for inew=1:length(new_aal_order)
    new_roi_label{inew} = roi_label{new_aal_order(inew)};
end
roi_label = new_roi_label;
roi_color = roi_color([idx_frontal,idx_occipital,idx_parietal,idx_temporal,idx_subcortical],:);

r = 1;
theta=linspace(0,2*pi,nAAL_ROI+1);
theta=theta(1:end-1);
x = r .* cos(theta);
y = r .* sin(theta);

biggest_amount_voxels = 0;

for iNet=1:nNet
    
    nSeeds = 1;
    
    net_info(iNet).network_label = Nets_Densities.net(iNet).network_label;
    net_info(iNet).color = color_label{iNet};
   
    aal_names = grouped_seeds_Nets_Densities.net(iNet).Net_Seeds{nSeeds+1,13};

    for iaal=1:size(aal_names,1)
            
        aal_roi_label = strrep(aal_names{iaal,1},'_','-');
        nVoxels = aal_names{iaal,2};
                
        if biggest_amount_voxels < nVoxels; biggest_amount_voxels = nVoxels; end
                
            for iiaal=1:length(roi_label)
                    
                if strcmp(aal_roi_label,roi_label(iiaal));
                        
                    net_info(iNet).net_aal_x(iaal) = x(iiaal) .* nVoxels;
                    net_info(iNet).net_aal_y(iaal) = y(iiaal) .* nVoxels;
                    net_info(iNet).net_aal_voxels(iaal) = nVoxels;
                    net_info(iNet).net_aal_theta(iaal) = theta(iiaal);
                    net_info(iNet).net_color{iaal} = color_label{iNet};
    
                end
                    
            end
                
    end

end

r_max = 1000;

all_theta = [];
all_net = [];
for iNet=1:nNet
    
    nAAL = length(net_info(iNet).net_aal_theta);

    for iAAL=1:nAAL
        
       all_theta = [all_theta, ones(1,net_info(iNet).net_aal_voxels(iAAL)).*net_info(iNet).net_aal_theta(iAAL)];
       all_net = [all_net, ones(1,net_info(iNet).net_aal_voxels(iAAL))*iNet];
        
    end 
    
end

[all_theta_s,I] = sort(all_theta);
all_net = all_net(I);

[uni_theta,count_theta] = count_unique(all_theta_s);

for iTheta=1:length(uni_theta)
    
    this_theta = uni_theta(iTheta);
   
    idx = find(all_theta_s==this_theta);
    
    idx_net = all_net(idx);
    
    [uni_net,count_net] = count_unique(idx_net);
    
    [count_net_s,I_net] = sort(count_net,'descend');
    
    for iNet=1:length(uni_net)
        
        n_voxels_net = round( r_max * sqrt( ( count_net_s(iNet) / biggest_amount_voxels  ) ) );
        
        [t,r] = rose(ones(1,n_voxels_net).*this_theta,nAAL_ROI);
        hold on
        set(gca,'ylim',[-r_max r_max]);
        set(gca,'xlim',[-r_max r_max]);
        
%         x_data = get(h,'Xdata');
%         y_data = get(h,'Ydata');
%         g = patch(x_data,y_data,'b');
%         set(g,'FaceColor',color(I(iNet),:));
%         
%         H = gca;
%         set(H,'ylim',[-r_max r_max]);
%         set(H,'xlim',[-r_max r_max]);

        h = polar(0,r_max);
        delete(h);
        set(gca,'Nextplot','add');
        
        [x_pol,y_pol] = pol2cart(t,r);
        g = patch(reshape(x_pol,4,[]),reshape(y_pol,4,[]),'b');
        set(g,'FaceColor',color(uni_net(I_net(iNet)),:));
        set(gca,'ylim',[-r_max r_max]);
        set(gca,'xlim',[-r_max r_max]);

    end
    
end

g = figure;
z = rose(ones(1,r_max),nAAL_ROI);
set(gca,'ylim',[-r_max r_max]);
set(gca,'xlim',[-r_max r_max]);

% all_theta = [];
% all_net = [];
% for iNet=1:nNet
%     
%     nAAL = length(net_info(iNet).net_aal_theta);
% 
%     for iAAL=1:nAAL
%         
%        all_theta = [all_theta, ones(1,net_info(iNet).net_aal_voxels(iAAL)).*net_info(iNet).net_aal_theta(iAAL)];
%        all_net = [all_net, ones(1,net_info(iNet).net_aal_voxels(iAAL))*iNet];
%         
%     end 
%     
% end
% 
% [all_theta_s,I] = sort(all_theta);
% all_net = all_net(I);
% 
% [uni_theta,count_theta] = count_unique(all_theta_s);
% 
% for iTheta=1:length(uni_theta)
%     
%     this_theta = uni_theta(iTheta);
%    
%     idx = find(all_theta_s==this_theta);
%     
%     idx_net = all_net(idx);
%     
%     [uni_net,count_net] = count_unique(idx_net);
%     
%     [count_net_s,I_net] = sort(count_net,'descend');
%     
%     for iNet=1:length(uni_net)
%         
%         n_voxels_net = round( r_max * sqrt( ( count_net_s(iNet) / biggest_amount_voxels  ) ) );
%         
%         [t,r] = rose(ones(1,n_voxels_net).*this_theta,nAAL_ROI);
%         hold on
%         set(gca,'ylim',[-r_max r_max]);
%         set(gca,'xlim',[-r_max r_max]);
%         
% %         x_data = get(h,'Xdata');
% %         y_data = get(h,'Ydata');
% %         g = patch(x_data,y_data,'b');
% %         set(g,'FaceColor',color(I(iNet),:));
% %         
% %         H = gca;
% %         set(H,'ylim',[-r_max r_max]);
% %         set(H,'xlim',[-r_max r_max]);
% 
%         h = polar(0,r_max);
%         delete(h);
%         set(gca,'Nextplot','add');
%         
%         [x_pol,y_pol] = pol2cart(t,r);
%         g = patch(reshape(x_pol,4,[]),reshape(y_pol,4,[]),'b');
%         set(g,'FaceColor',color(uni_net(I_net(iNet)),:));
%         set(gca,'ylim',[-r_max r_max]);
%         set(gca,'xlim',[-r_max r_max]);
% 
%     end
%     
% end
% 
% g = figure;
% z = rose(ones(1,r_max),nAAL_ROI);
% set(gca,'ylim',[-r_max r_max]);
% set(gca,'xlim',[-r_max r_max]);

% ylim([-biggest_amount_voxels, biggest_amount_voxels]);
% xlim([-biggest_amount_voxels, biggest_amount_voxels]);
%     
% hold on
%     
x_new = x .* r_max;
y_new = y .* r_max;

for i=1:length(x)-1
    line([x_new(i) x_new(i+1)], [y_new(i) y_new(i+1)], 'LineWidth', 15, 'Color',roi_color(i,:));
    hold on
end

for i=1:length(x)
    h = text(x_new(i).*1.05,y_new(i).*1.05,roi_label{i},'FontSize',8);
    set(h,'Rotation',(theta(i)*180/pi));
end
        

end

function plotCircularRelationship_v4

load('Non-Parametric-Voxel-Wise-Per-Density-AAL-FuncNets-AAL-IDs-grouped-seeds.mat');

%%% colors

color(1,:) = [0,0,128];
color(2,:) = [128,0,230];
color(3,:) = [0,191,255];
color(4,:) = [0,100,0];
color(5,:) = [255,255,0];
color(6,:) = [232,105,43];
color(8,:) = [255,0,0];
color(7,:) = [0,0,0];

color_label{1} = 'b';
color_label{2} = 'm';
color_label{3} = 'c';
color_label{4} = 'g';
color_label{5} = 'y';
color_label{6} = 'g';
color_label{8} = 'b';
color_label{7} = 'r';

color = color./255;

%%% LOAD NETS

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

nAAL_ROI = 90;
nNet = 8;

% rgb = round(linspace(0,255,nLobesAndSubCortical+1))./255;
% rgb = rgb(1:end-1);
rgb(1) = 192;
rgb(2) = 128;
rgb(3) = 32;
rgb(4) = 224;
rgb(5) = 160;
rgb = rgb./255;

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

for iROI=1:nAAL_ROI
    
    roi_label{iROI} = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    
    if ismember(iROI,idx_frontal)
        
        aal_color = 1;
        
    elseif ismember(iROI,idx_occipital)
        
        aal_color = 2;
        
    elseif ismember(iROI,idx_parietal)
        
        aal_color = 3;
        
    elseif ismember(iROI,idx_temporal)
        
        aal_color = 4;
        
    elseif ismember(iROI,idx_subcortical)
        
        aal_color = 5;
        
    end
    
    roi_color(iROI,:) = [rgb(aal_color), rgb(aal_color), rgb(aal_color)];
    
end

%roi_label = roi_label{[idx_frontal,idx_occipital,idx_parietal,idx_temporal,idx_subcortical]};
%new_aal_order = [idx_frontal,idx_occipital,idx_parietal,idx_temporal,idx_subcortical];
new_aal_order = [idx_frontal(1:2:end-1),idx_occipital(1:2:end-1),idx_parietal(1:2:end-1),idx_temporal(1:2:end-1),idx_subcortical(1:2:end-1),idx_frontal(2:2:end),idx_occipital(2:2:end),idx_parietal(2:2:end),idx_temporal(2:2:end),idx_subcortical(2:2:end)];
for inew=1:length(new_aal_order)
    new_roi_label{inew} = roi_label{new_aal_order(inew)};
end
roi_label = new_roi_label;
roi_color = roi_color([idx_frontal,idx_occipital,idx_parietal,idx_temporal,idx_subcortical],:);

r = 1;
theta=linspace(0,2*pi,nAAL_ROI+1);
theta=theta(1:end-1);
x = r .* cos(theta);
y = r .* sin(theta);

for iTheta=1:length(theta)
    
    ROI_net(iTheta).net = 0;
    ROI_net(iTheta).nVoxels = 0;

end

biggest_amount_voxels = 0;

for iNet=1:nNet
    
    nSeeds = 1;
    
    net_info(iNet).network_label = Nets_Densities.net(iNet).network_label;
    net_info(iNet).color = color_label{iNet};
   
    aal_names = grouped_seeds_Nets_Densities.net(iNet).Net_Seeds{nSeeds+1,13};

    for iaal=1:size(aal_names,1)
            
        aal_roi_label = strrep(aal_names{iaal,1},'_','-');
        nVoxels = aal_names{iaal,2};
                
        if biggest_amount_voxels < nVoxels; biggest_amount_voxels = nVoxels; end
                
            for iiaal=1:length(roi_label)
                    
                if strcmp(aal_roi_label,roi_label(iiaal));
                        
                    net_info(iNet).net_aal_x(iaal) = x(iiaal) .* nVoxels;
                    net_info(iNet).net_aal_y(iaal) = y(iiaal) .* nVoxels;
                    net_info(iNet).net_aal_voxels(iaal) = nVoxels;
                    net_info(iNet).net_aal_theta(iaal) = theta(iiaal);
                    net_info(iNet).net_color{iaal} = color_label{iNet};
                    
                    ROI_net(iiaal).net(end+1) = iNet;
                    ROI_net(iiaal).nVoxels(end+1) = nVoxels;
    
                end
                    
            end
                
    end

end

r_max = 1000;

for iTheta=1:length(theta)
    
    ROI_net(iTheta).net(ROI_net(iTheta).net == 0) = [];
    ROI_net(iTheta).nVoxels(ROI_net(iTheta).nVoxels == 0) = [];
    
    if ~isempty(ROI_net(iTheta).net)
    
        this_theta_nets = ROI_net(iTheta).net;
        this_theta_nvoxels = ROI_net(iTheta).nVoxels;
    
        [nVoxels_s, I] = sort(this_theta_nvoxels,'descend');
        Nets_s = this_theta_nets(I);

        for iNet=1:length(Nets_s)

            n_voxels_net = round( r_max * sqrt( ( nVoxels_s(iNet) / biggest_amount_voxels  ) ) );

            [t,r] = rose(ones(1,n_voxels_net).*theta(iTheta),nAAL_ROI);
            hold on
            set(gca,'ylim',[-r_max r_max]);
            set(gca,'xlim',[-r_max r_max]);

            h = polar(0,r_max);
            delete(h);
            set(gca,'Nextplot','add');

            [x_pol,y_pol] = pol2cart(t,r);
            g = patch(reshape(x_pol,4,[]),reshape(y_pol,4,[]),'b');
            set(g,'FaceColor',color(Nets_s(iNet),:));
            set(gca,'ylim',[-r_max r_max]);
            set(gca,'xlim',[-r_max r_max]);

        end
        
    end
    
end

g = figure;
z = rose(ones(1,r_max),nAAL_ROI);
set(gca,'ylim',[-r_max r_max]);
set(gca,'xlim',[-r_max r_max]);

x_new = x .* r_max;
y_new = y .* r_max;

for i=1:length(x)-1
    line([x_new(i) x_new(i+1)], [y_new(i) y_new(i+1)], 'LineWidth', 15, 'Color',roi_color(i,:));
    hold on
end

for i=1:length(x)
    h = text(x_new(i).*1.05,y_new(i).*1.05,roi_label{i},'FontSize',8);
    set(h,'Rotation',(theta(i)*180/pi));
end
        

end


