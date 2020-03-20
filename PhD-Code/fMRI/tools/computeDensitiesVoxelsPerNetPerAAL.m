
clear all

% idx_density = 13; %%% Attention Negative Up
% idx_density = 15; %%% Attention Positive Down
idx_density = 17; %%% Attention Positive Up
% idx_density_stat = 3; %%% Attention Negative Up - stat
% idx_density_stat = 4; %%% Attention Positive Down - stat
idx_density_stat = 5; %%% Attention Positive Up - stat

load('Non-Parametric-Voxel-Wise-Per-Density-AAL-FuncNets-AAL-IDs-grouped-seeds.mat');

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

nAAL_ROI = 90;
nNet = 8;

color(1,:) = [0,0,128];
color(2,:) = [128,0,230];
color(3,:) = [0,191,255];
color(4,:) = [0,100,0];
color(5,:) = [255,255,0];
color(6,:) = [232,105,43];
color(7,:) = [0,0,0];
color(8,:) = [255,0,0];

color = color./255;

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
new_aal_order = [idx_frontal(1:2:end-1),idx_occipital(1:2:end-1),idx_parietal(1:2:end-1),idx_temporal(1:2:end-1),idx_subcortical(1:2:end-1),fliplr(idx_subcortical(2:2:end)),fliplr(idx_temporal(2:2:end)),fliplr(idx_parietal(2:2:end)),fliplr(idx_occipital(2:2:end)),fliplr(idx_frontal(2:2:end))];
for inew=1:length(new_aal_order)
    new_roi_label{inew} = roi_label{new_aal_order(inew)};
end
roi_label = new_roi_label;
roi_color = roi_color(new_aal_order,:);

theta = linspace(0,2*pi,nAAL_ROI+1);
theta_net = theta;
theta = theta(1:end-1);

roi_theta = theta(new_aal_order);

densities{1,1} = 'AAL';

for iNet=1:nNet
    
    densities{1,1+iNet} = all_networks_labels{iNet};
    
end

for iROI=1:nAAL_ROI
    
    densities{1+iROI,1} = strrep(AAL_ROI(iROI).Nom_L,'_','-');
    ROI_labels_normal{iROI} = strrep(AAL_ROI(iROI).Nom_L,'_','-');

end

for iROI=1:nAAL_ROI
    
    ROI_net(iROI).net(1) = 0;
    ROI_net(iROI).nVoxels(1) = 0;

end

for iNet=1:nNet
    
    if Nets_Densities.net(iNet).StatsResults{idx_density_stat,2} == 1
        
        nSeeds = 1;

        net_info(iNet).network_label = Nets_Densities.net(iNet).network_label;

        aal_names = grouped_seeds_Nets_Densities.net(iNet).Net_Seeds{nSeeds+1,idx_density};

        for iaal=1:size(aal_names,1)

            aal_roi_label = strrep(aal_names{iaal,1},'_','-');
            nVoxels = aal_names{iaal,2};

                for iiaal=1:length(roi_label)

                    if strcmp(aal_roi_label,ROI_labels_normal(iiaal));
                        
                        ROI_net(iiaal).net(end+1) = iNet;
                        ROI_net(iiaal).nVoxels(end+1) = nVoxels;

                    end

                end
                
        end
    
    end

end

for iROI=1:nAAL_ROI
    
    ROI_net(iROI).net(1) = [];
    ROI_net(iROI).nVoxels(1) = [];

end

for iROI=1:nAAL_ROI
   
    ROI_net(iROI).densities = ROI_net(iROI).nVoxels ./ sum(ROI_net(iROI).nVoxels);
    
end

for iROI=1:nAAL_ROI
   
    n_this_roi_nets = length(ROI_net(iROI).net);

    for iNet=1:n_this_roi_nets
        
        densities{1+iROI,1+ROI_net(iROI).net(iNet)} = ROI_net(iROI).densities(iNet);
        
    end
        
end

count_net = zeros(1,nNet);

for iROI=1:nAAL_ROI

    for iNet=1:nNet
        
        if ~isempty(densities{1+iROI,1+iNet})
            
            count_net(iNet) = count_net(iNet) + 1;
        
        end
        
    end
        
end

[count_net_s,I] = sort(count_net,'ascend');

jump = 0;

is_there_change = zeros(1,nAAL_ROI);

for iNet=1:nNet
    
    r1 = iNet + jump;
    
    x1 = r1*cos(theta_net);
    y1 = r1*sin(theta_net);
    
    plot(x1,y1,'color',color(I(iNet),:),'LineWidth',2);
    hold on
    
    for ijump=0.01:0.01:1
        
        r2 = iNet + ijump + jump;
    
        x2 = r2*cos(theta_net);
        y2 = r2*sin(theta_net);
    
        plot(x2,y2,'color',color(I(iNet),:),'LineWidth',2);
        hold on
 
    end
    
    jump = jump + 1;
    
    r_max = r2;
    
    for iROI=1:nAAL_ROI
    
        if ~isempty(densities{1+iROI,1+I(iNet)})
            
            idx_theta = find(new_aal_order == iROI);
            
            if densities{1+iROI,1+I(iNet)} == 1; marker = '*'; end
            if densities{1+iROI,1+I(iNet)} < 1 && densities{1+iROI,1+I(iNet)} > 0.75; marker = '*'; end
            if densities{1+iROI,1+I(iNet)} == 0.75; marker = 'x'; end
            if densities{1+iROI,1+I(iNet)} < 0.75 && densities{1+iROI,1+I(iNet)} > 0.5; marker = 'x'; end
            if densities{1+iROI,1+I(iNet)} == 0.5; marker = '+'; end
            if densities{1+iROI,1+I(iNet)} < 0.5 && densities{1+iROI,1+I(iNet)} > 0.25;  marker = '+'; end
            if densities{1+iROI,1+I(iNet)} == 0.25; marker = 'o'; end
            if densities{1+iROI,1+I(iNet)} < 0.25; marker = 'o'; end
            
            if I(iNet) == 5; marker_color = [0 0 0]; else marker_color = [255 255 255]./255; end
            
            plot(((x2(idx_theta)-x1(idx_theta))/2) + x1(idx_theta),((y2(idx_theta)-y1(idx_theta))/2) + y1(idx_theta),marker,'color',marker_color,'MarkerSize',13);
        
            is_there_change(idx_theta) = 1;
            
        end
         
    end
    
end

r = 1;
x = r .* cos(theta);
y = r .* sin(theta);

x_new = x .* ( r_max + 1);
y_new = y .* ( r_max + 1);

for i=1:length(x)-1
    line([x_new(i) x_new(i+1)], [y_new(i) y_new(i+1)], 'LineWidth', 15, 'Color',roi_color(i,:));
    hold on
end

for i=1:length(x)
    h = text(x_new(i).*1.05,y_new(i).*1.05,roi_label{i},'FontSize',16);
    set(h,'Rotation',(theta(i)*180/pi));
end

for i=1:length(is_there_change)
   
    if is_there_change(i)
        
        line([0 x_new(i)],[0 y_new(i)],'LineWidth',0.5,'Color',[0 0 0]);
        
    end
    
end



