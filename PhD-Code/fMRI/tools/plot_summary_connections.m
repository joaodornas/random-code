
clear all

load('FC-Voxels-AAL-ROI-corr-KMeans-MeaningFul-FC-Att-Stim-Only.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Functional-Parcels.mat');

nClusters = 758;

nROI = 90;
for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

FC_strength = zeros(1,length(MeaningFul_Attention)-1);
GC_FROM_strength = zeros(1,length(MeaningFul_Attention)-1);
GC_TO_strength = zeros(1,length(MeaningFul_Attention)-1);

for iMean=2:length(MeaningFul_Attention)
    
    idx_cluster(iMean-1) = MeaningFul_Attention{iMean,1};
    
    cluster_label{iMean-1} = MeaningFul_Attention{iMean,2};

    if ~isempty(length(MeaningFul_Attention{iMean,7}))
        FC_strength(iMean-1) = length(MeaningFul_Attention{iMean,7}); % FC - direct
    else
        FC_strength(iMean-1) = 0;
    end
    if ~isempty(length(MeaningFul_Attention{iMean,8}))
        FC_strength(iMean-1) = FC_strength(iMean-1) + length(MeaningFul_Attention{iMean,8}); % FC - indirect
    else
        FC_strength(iMean-1) = FC_strength(iMean-1) + 0;
    end

    if ~isempty(length(MeaningFul_Attention{iMean,11}))
        GC_FROM_strength(iMean-1) = length(MeaningFul_Attention{iMean,11}); % GC - direct - FROM
    else
        GC_FROM_strength(iMean-1) = 0;
    end
    if ~isempty(length(MeaningFul_Attention{iMean,12}))
        GC_FROM_strength(iMean-1) = GC_FROM_strength(iMean-1) + length(MeaningFul_Attention{iMean,12}); % GC - indirect - FROM
    else
        GC_FROM_strength(iMean-1) = GC_FROM_strength(iMean-1) + 0;
    end
    
    if ~isempty(length(MeaningFul_Attention{iMean,15}))
        GC_TO_strength(iMean-1) = length(MeaningFul_Attention{iMean,15}); % GC - direct - TO
    else
        GC_TO_strength(iMean-1) = 0;
    end
    if~isempty(length(MeaningFul_Attention{iMean,16}))
        GC_TO_strength(iMean-1) = GC_TO_strength(iMean-1) + length(MeaningFul_Attention{iMean,16}); % GC - indirect - TO
    else
        GC_TO_strength(iMean-1) = GC_TO_strength(iMean-1) + 0;
    end

end

FC_strength = FC_strength ./ max(FC_strength);

GC_FROM_percent = GC_FROM_strength ./ (GC_FROM_strength + GC_TO_strength);
GC_TO_percent = GC_TO_strength ./ (GC_FROM_strength + GC_TO_strength);

[u_cluster_label, c_cluster_label] = count_unique(cluster_label);

for iUnique=1:length(u_cluster_label)

    for iROI=1:nROI
           
       if strcmp(ROI(iROI).label,u_cluster_label(iUnique))
           
           n_clusters(iUnique) = ROI(iROI).nClusters;
       
       end
       
    end
    
end

f = figure;

r = 1;
angles = linspace(0,2*pi,sum(n_clusters)+1);
angles(end) = [];
x = r.*cos(angles);
y = r.*sin(angles);

% plot(x,y,'k');
% hold on

plot(x,y,'k*','MarkerSize',2);
hold on

for iUnique=1:length(u_cluster_label)
   
    if iUnique > 1
        
        last_clusters = sum(n_clusters(1:iUnique-1));
        
    else
        
        last_clusters = 0;
        
    end
        
    position = round(n_clusters(iUnique)/2) + last_clusters;
   
    x_label = (r+0.20).*cos(angles(position));
    y_label = (r+0.20).*sin(angles(position));
    h = text(x_label,y_label,u_cluster_label(iUnique));
    set(h,'Rotation',angles(position)*180/pi);
    
    position = n_clusters(iUnique) + last_clusters;
     
    step = angles(2) - angles(1);
    
    x_division_up = (r+0.05).*cos(angles(position) + step/2);
    x_division = r.*cos(angles(position) + step/2);
    y_division_up = (r+0.05).*sin(angles(position) + step/2);
    y_division = r.*sin(angles(position) + step/2);
    
    plot([x_division x_division_up],[y_division y_division_up],'k','LineWidth',2);
    hold on
    
end

for iMean=1:length(idx_cluster)
    
    for iROI=1:nROI
        
        if idx_cluster(iMean) >= ROI_info{iROI,4} && idx_cluster(iMean) <= ROI_info{iROI,5}
        
            order_clusters(iMean) = idx_cluster(iMean) - ROI_info{iROI,4} + 1;
            
        end
        
    end
    
end

reduction = 0.2;

for iUnique=1:length(u_cluster_label)
    
   for iCluster=1:c_cluster_label(iUnique)
   
       if iUnique > 1
           
           start = sum(c_cluster_label(1:iUnique-1));
           
       else
           
           start = 0;
       
       end
       
       order = order_clusters(start+iCluster);
       
       position_hubs = start+iCluster;
       
       if iUnique > 1
           
           position = sum(n_clusters(1:iUnique-1)) + order;
       
       else
           
           position = order;
           
       end
       
        x_up = (r+GC_FROM_percent(position_hubs)*reduction).*cos(angles(position));
        y_up = (r+GC_FROM_percent(position_hubs)*reduction).*sin(angles(position));
        
        x_down = (r-GC_TO_percent(position_hubs)*reduction).*cos(angles(position));
        y_down = (r-GC_TO_percent(position_hubs)*reduction).*sin(angles(position));
        
        if FC_strength(position_hubs) ~=0
            plot([x_up x_down],[y_up y_down],'b','LineWidth',FC_strength(position_hubs)*3.5);
            hold on
        end
        
   end
    
end

axis off


    