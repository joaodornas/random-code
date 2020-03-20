
% function check_meaningful_clusters

clear all

% load('FC-Voxels-AAL-ROI-corr-KMeans-MeaningFul-FC-GC-Att-Stim-Only.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-FC-All-Clusters-Contrast.mat');
% load('FC_Voxel_AAL_ROI_kmeans_Granger_Clusters-Mean-Contrast.mat');
load('FC-Voxels-AAL-ROI-corr-KMeans-Info-Functional-Parcels.mat');

nClusters = 758;
nCon = 6;
zcriterion = 2.3;

con_color{1} = 'b';
con_color{2} = 'r';
con_color{3} = 'y';
con_color{4} = 'm';
con_color{5} = 'g';
con_color{6} = 'k';

FC_Contrast_Attention = contrast_attention_z;
GC_Contrast_Attention = Attention_Contrast.Z;
FC_Contrast_Stimulus = contrast_stimulus_z;
GC_Contrast_Stimulus = Stimulus_Contrast.Z;

FC_Contrast_Attention(isnan(FC_Contrast_Attention)) = 0;
GC_Contrast_Attention(isnan(GC_Contrast_Attention)) = 0;
FC_Contrast_Attention(find(FC_Contrast_Attention>(-1)*zcriterion & FC_Contrast_Attention<zcriterion)) = 0; 
GC_Contrast_Attention(find(GC_Contrast_Attention>(-1)*zcriterion & GC_Contrast_Attention<zcriterion)) = 0; 

FC_Contrast_Stimulus(isnan(FC_Contrast_Stimulus)) = 0;
GC_Contrast_Stimulus(isnan(GC_Contrast_Stimulus)) = 0;
FC_Contrast_Stimulus(find(FC_Contrast_Stimulus>(-1)*zcriterion & FC_Contrast_Stimulus<zcriterion)) = 0; 
GC_Contrast_Stimulus(find(GC_Contrast_Stimulus>(-1)*zcriterion & GC_Contrast_Stimulus<zcriterion)) = 0; 

FC_Attention_Only = FC_Contrast_Attention;
FC_Attention_Only(find(FC_Contrast_Stimulus)) = 0;
FC_Stimulus_Only = FC_Contrast_Stimulus;
FC_Stimulus_Only(find(FC_Contrast_Attention)) = 0;
FC_Common = FC_Contrast_Attention & FC_Contrast_Stimulus;

FC_Common_AiSi = FC_Contrast_Attention > 0 & FC_Contrast_Stimulus > 0;
FC_Common_AiSd = FC_Contrast_Attention > 0 & FC_Contrast_Stimulus < 0;
FC_Common_AdSi = FC_Contrast_Attention < 0 & FC_Contrast_Stimulus > 0;
FC_Common_AdSd = FC_Contrast_Attention < 0 & FC_Contrast_Stimulus < 0;

GC_Attention_Only = GC_Contrast_Attention;
GC_Attention_Only(find(GC_Contrast_Stimulus)) = 0;
GC_Stimulus_Only = GC_Contrast_Stimulus;
GC_Stimulus_Only(find(GC_Contrast_Attention)) = 0;
GC_Common = GC_Contrast_Attention & GC_Contrast_Stimulus;

GC_Common_AiSi = GC_Contrast_Attention > 0 & GC_Contrast_Stimulus > 0;
GC_Common_AiSd = GC_Contrast_Attention > 0 & GC_Contrast_Stimulus < 0;
GC_Common_AdSi = GC_Contrast_Attention < 0 & GC_Contrast_Stimulus > 0;
GC_Common_AdSd = GC_Contrast_Attention < 0 & GC_Contrast_Stimulus < 0;

for iCluster=1:nClusters
    
   vector_FC_A = FC_Attention_Only(iCluster,:);
   vector_FC_S = FC_Stimulus_Only(iCluster,:);
   vector_FC_C = FC_Common(iCluster,:);
   
   n_FC_A_p = length(find(vector_FC_A > 0));
   n_FC_A_n = length(find(vector_FC_A < 0));
   n_FC_S_p = length(find(vector_FC_S > 0));
   n_FC_S_n = length(find(vector_FC_S < 0));
   n_FC_C_p = length(find(vector_FC_C > 0));
   
   cluster_FC(1,iCluster) = n_FC_A_p;
   cluster_FC(2,iCluster) = n_FC_A_n;
   cluster_FC(3,iCluster) = n_FC_S_p;
   cluster_FC(4,iCluster) = n_FC_S_n;
   cluster_FC(5,iCluster) = n_FC_C_p;

   vector_GC_A = GC_Attention_Only(iCluster,:);
   vector_GC_S = GC_Stimulus_Only(iCluster,:);
   vector_GC_C = GC_Common(iCluster,:);
   
   n_GC_A_p = length(find(vector_GC_A > 0));
   n_GC_A_n = length(find(vector_GC_A < 0));
   n_GC_S_p = length(find(vector_GC_S > 0));
   n_GC_S_n = length(find(vector_GC_S < 0));
   n_GC_C_p = length(find(vector_GC_C > 0));
   
   cluster_GC_FROM(1,iCluster) = n_GC_A_p;
   cluster_GC_FROM(2,iCluster) = n_GC_A_n;
   cluster_GC_FROM(3,iCluster) = n_GC_S_p;
   cluster_GC_FROM(4,iCluster) = n_GC_S_n;
   cluster_GC_FROM(5,iCluster) = n_GC_C_p;

   vector_GC_A = GC_Attention_Only(:,iCluster);
   vector_GC_S = GC_Stimulus_Only(:,iCluster);
   vector_GC_C = GC_Common(:,iCluster);
   
   n_GC_A_p = length(find(vector_GC_A > 0));
   n_GC_A_n = length(find(vector_GC_A < 0));
   n_GC_S_p = length(find(vector_GC_S > 0));
   n_GC_S_n = length(find(vector_GC_S < 0));
   n_GC_C_p = length(find(vector_GC_C > 0));
  
   cluster_GC_TO(1,iCluster) = n_GC_A_p;
   cluster_GC_TO(2,iCluster) = n_GC_A_n;
   cluster_GC_TO(3,iCluster) = n_GC_S_p;
   cluster_GC_TO(4,iCluster) = n_GC_S_n;
   cluster_GC_TO(5,iCluster) = n_GC_C_p;

   vector_FC_Common_AiSi = FC_Common_AiSi(iCluster,:);
   vector_FC_Common_AiSd = FC_Common_AiSd(iCluster,:);
   vector_FC_Common_AdSi = FC_Common_AdSi(iCluster,:);
   vector_FC_Common_AdSd = FC_Common_AdSd(iCluster,:);
   
   n_FC_Common_AiSi = length(find(vector_FC_Common_AiSi > 0));
   n_FC_Common_AiSd = length(find(vector_FC_Common_AiSd > 0));
   n_FC_Common_AdSi = length(find(vector_FC_Common_AdSi > 0));
   n_FC_Common_AdSd = length(find(vector_FC_Common_AdSd > 0));
   
   cluster_FC_Common(1,iCluster) = n_FC_Common_AiSi;
   cluster_FC_Common(2,iCluster) = n_FC_Common_AiSd;
   cluster_FC_Common(3,iCluster) = n_FC_Common_AdSi;
   cluster_FC_Common(4,iCluster) = n_FC_Common_AdSd;
   
   vector_GC_Common_FROM_AiSi = GC_Common_AiSi(iCluster,:);
   vector_GC_Common_FROM_AiSd = GC_Common_AiSd(iCluster,:);
   vector_GC_Common_FROM_AdSi = GC_Common_AdSi(iCluster,:);
   vector_GC_Common_FROM_AdSd = GC_Common_AdSd(iCluster,:);
   
   n_GC_Common_FROM_AiSi = length(find(vector_GC_Common_FROM_AiSi > 0));
   n_GC_Common_FROM_AiSd = length(find(vector_GC_Common_FROM_AiSd > 0));
   n_GC_Common_FROM_AdSi = length(find(vector_GC_Common_FROM_AdSi > 0));
   n_GC_Common_FROM_AdSd = length(find(vector_GC_Common_FROM_AdSd > 0));
   
   cluster_GC_Common_FROM(1,iCluster) = n_GC_Common_FROM_AiSi;
   cluster_GC_Common_FROM(2,iCluster) = n_GC_Common_FROM_AiSd;
   cluster_GC_Common_FROM(3,iCluster) = n_GC_Common_FROM_AdSi;
   cluster_GC_Common_FROM(4,iCluster) = n_GC_Common_FROM_AdSd;
   
   vector_GC_Common_TO_AiSi = GC_Common_AiSi(:,iCluster);
   vector_GC_Common_TO_AiSd = GC_Common_AiSd(:,iCluster);
   vector_GC_Common_TO_AdSi = GC_Common_AdSi(:,iCluster);
   vector_GC_Common_TO_AdSd = GC_Common_AdSd(:,iCluster);
   
   n_GC_Common_TO_AiSi = length(find(vector_GC_Common_TO_AiSi > 0));
   n_GC_Common_TO_AiSd = length(find(vector_GC_Common_TO_AiSd > 0));
   n_GC_Common_TO_AdSi = length(find(vector_GC_Common_TO_AdSi > 0));
   n_GC_Common_TO_AdSd = length(find(vector_GC_Common_TO_AdSd > 0));
   
   cluster_GC_Common_TO(1,iCluster) = n_GC_Common_TO_AiSi;
   cluster_GC_Common_TO(2,iCluster) = n_GC_Common_TO_AiSd;
   cluster_GC_Common_TO(3,iCluster) = n_GC_Common_TO_AdSi;
   cluster_GC_Common_TO(4,iCluster) = n_GC_Common_TO_AdSd;
   
end

cluster_FC_Common = cluster_FC_Common';
cluster_GC_Common_FROM = cluster_GC_Common_FROM';
cluster_GC_Common_TO = cluster_GC_Common_TO';

con_label{1} = 'Attention (Only) - Increase';
con_label{2} = 'Attention (Only) - Decrease';
con_label{3} = 'Stimulus (Only) - Increase';
con_label{4} = 'Stimulus (Only) - Decrease';
con_label{5} = 'Common';

common_label{1} = 'FC-AiSi';
common_label{2} = 'FC-AiSd';
common_label{3} = 'FC-AdSi';
common_label{4} = 'FC-AiSd';

f = figure;

for icon=1:nCon-1
    con_leg{icon} = strcat('FC:',con_label{icon});
end

step = 0.001;
max_v = 0;
for icon=1:nCon-1
    vector = cluster_FC(icon,:);
    if max(vector) > max_v; max_v = max(vector); end
end
interval = 1:step:max_v;

for icon=1:nCon-1
    
    vector = cluster_FC(icon,:);
    vector(vector==0) = [];

    if ~isempty(vector)

        k = ksdensity(vector(:),interval);
        probability_distribution = k .* step;
        idx = find(probability_distribution==0);
        probability_distribution(idx) = [];
        interval(idx) = [];
        plot(interval,probability_distribution,con_color{icon});

        hold on

    else

        con_leg(icon) = [];

    end
    
end

xlabel('# of connections');
ylabel('# of clusters');
legend(con_leg);
title('FC');

print(f,'-depsc','MeaningfulClusters-Probability-Density-FC-Connections.eps');

f = figure;

for icon=1:nCon-2
    
    subplot(2,2,icon);

    vector = cluster_FC(icon,:);
    vector(vector==0) = [];

    if ~isempty(vector)

        histfit(vector,max(vector),'gamma');
        xlabel('# of connections');
        ylabel('# of clusters');
        title(con_leg{icon});

    end

end

suptitle('FC');

print(f,'-depsc','MeaningfulClusters-Contrasts-Distributions-FC-Connections.eps');

f = figure;
hist(cluster_FC_Common,1:max(cluster_FC_Common(:)));
legend(common_label,'Location','NorthOutside');
xlabel('# of connections');
ylabel('# of clusters (log)');
title('FC');
%set(gca,'YScale','log');

print(f,'-depsc','MeaningfulClusters-Contrasts-Distributions-Common-FC-Connections.eps');

cluster_FC_Common_NoZeros = cluster_FC_Common(:);
cluster_FC_Common_NoZeros(cluster_FC_Common_NoZeros==0) = [];

f = figure;
histfit(cluster_FC_Common_NoZeros,max(cluster_FC_Common_NoZeros(:)),'gamma');
% legend(common_label,'Location','NorthOutside');
xlabel('# of connections');
ylabel('# of clusters (log)');
title('FC');
%set(gca,'YScale','log');

print(f,'-depsc','MeaningfulClusters-Contrasts-Distributions-Common-FC-Connections-noZeros.eps');

con_label{1} = 'Attention (Only) - Increase';
con_label{2} = 'Attention (Only) - Decrease';
con_label{3} = 'Stimulus (Only) - Increase';
con_label{4} = 'Stimulus (Only) - Decrease';
con_label{5} = 'Common';

common_label{1} = 'GC-FROM-AiSi';
common_label{2} = 'GC-FROM-AiSd';
common_label{3} = 'GC-FROM-AdSi';
common_label{4} = 'GC-FROM-AiSd';

f = figure;

for icon=1:nCon-1
    con_leg{icon} = strcat('GC-FROM:',con_label{icon});
end

step = 0.001;
max_v = 0;
for icon=1:nCon-1
    vector = cluster_GC_FROM(icon,:);
    if max(vector) > max_v; max_v = max(vector); end
end
interval = 1:step:max_v;

for icon=1:nCon-1
    
    vector = cluster_GC_FROM(icon,:);
    vector(vector==0) = [];

    if ~isempty(vector)
        
        k = ksdensity(vector(:),interval);
        probability_distribution = k .* step;
        idx = find(probability_distribution==0);
        probability_distribution(idx) = [];
        interval(idx) = [];
        plot(interval,probability_distribution,con_color{icon});

        hold on
    
    else
        
        con_leg(icon) = [];
        
    end
 
end

xlabel('# of connections');
ylabel('# of clusters');
legend(con_leg);
title('GC-FROM');

print(f,'-depsc','MeaningfulClusters-Probability-Density-GC-FROM-Connections.eps');

f = figure;

for icon=1:nCon-2
    
    subplot(2,2,icon);
    
    vector = cluster_GC_FROM(icon,:);
    vector(vector==0) = [];

    if ~isempty(vector)

        histfit(vector,max(vector),'gamma');
        xlabel('# of connections');
        ylabel('# of clusters');
        title(con_leg{icon});

    end
    
end
    
suptitle('GC-FROM');

print(f,'-depsc','MeaningfulClusters-Contrasts-Distributions-GC-FROM-Connections.eps');

f = figure;

hist(cluster_GC_Common_FROM,1:max(cluster_GC_Common_FROM(:)));
legend(common_label,'Location','NorthOutside');
xlabel('# of connections');
ylabel('# of clusters - log');
title('GC-FROM');
%set(gca,'YScale','log');

print(f,'-depsc','MeaningfulClusters-Contrasts-Distributions-Common-GC-FROM-Connections.eps');

cluster_GC_Common_FROM_NoZeros = cluster_GC_Common_FROM(:);
cluster_GC_Common_FROM_NoZeros(cluster_GC_Common_FROM_NoZeros==0) = [];

f = figure;
histfit(cluster_GC_Common_FROM_NoZeros,max(cluster_GC_Common_FROM_NoZeros(:)),'gamma');
% legend(common_label,'Location','NorthOutside');
xlabel('# of connections');
ylabel('# of clusters (log)');
title('GC-FROM');
%set(gca,'YScale','log');

print(f,'-depsc','MeaningfulClusters-Contrasts-Distributions-Common-GC-FROM-Connections-noZeros.eps');

con_label{1} = 'Attention (Only) - Increase';
con_label{2} = 'Attention (Only) - Decrease';
con_label{3} = 'Stimulus (Only) - Increase';
con_label{4} = 'Stimulus (Only) - Decrease';
con_label{5} = 'Common';

common_label{1} = 'GC-TO-AiSi';
common_label{2} = 'GC-TO-AiSd';
common_label{3} = 'GC-TO-AdSi';
common_label{4} = 'GC-TO-AiSd';

f = figure;

for icon=1:nCon-1
    con_leg{icon} = strcat('GC-TO:',con_label{icon});
end

step = 0.001;
max_v = 0;
for icon=1:nCon-1
    vector = cluster_GC_TO(icon,:);
    if max(vector) > max_v; max_v = max(vector); end
end
interval = 1:step:max_v;

for icon=1:nCon-1
    
    vector = cluster_GC_TO(icon,:);
    vector(vector==0) = [];
    
    if ~isempty(vector)
        
        k = ksdensity(vector(:),interval);
        probability_distribution = k .* step;
        idx = find(probability_distribution==0);
        probability_distribution(idx) = [];
        interval(idx) = [];
        plot(interval,probability_distribution,con_color{icon});

        hold on
    
    else
        
        con_leg(icon) = [];
        
    end
    
end

xlabel('# of connections');
ylabel('# of clusters');
legend(con_leg);
title('GC-TO');

print(f,'-depsc','MeaningfulClusters-Probability-Distribution-GC-TO-Connections.eps');

f = figure;

for icon=1:nCon-2
    
    subplot(2,2,icon);
    
    vector = cluster_GC_TO(icon,:);
    vector(vector==0) = [];

    if ~isempty(vector)

        histfit(vector,max(vector),'gamma');
        xlabel('# of connections');
        ylabel('# of clusters');
        title(con_leg{icon});

    end
    
end

suptitle('GC-TO');

print(f,'-depsc','MeaningfulClusters-Contrasts-Distributions-GC-TO-Connections.eps');

f = figure;
        
hist(cluster_GC_Common_TO,1:max(cluster_GC_Common_TO(:)));
legend(common_label,'Location','NorthOutside');
xlabel('# of connections');
ylabel('# of clusters - log');
title('GC-TO');
%set(gca,'YScale','log');
        
suptitle('GC-TO');

print(f,'-depsc','MeaningfulClusters-Contrasts-Distributions-Common-GC-TO-Connections.eps');

cluster_GC_Common_TO_NoZeros = cluster_GC_Common_TO(:);
cluster_GC_Common_TO_NoZeros(cluster_GC_Common_TO_NoZeros==0) = [];

f = figure;
histfit(cluster_GC_Common_TO_NoZeros,max(cluster_GC_Common_TO_NoZeros(:)),'gamma');
% legend(common_label,'Location','NorthOutside');
xlabel('# of connections');
ylabel('# of clusters (log)');
title('GC-TO');
%set(gca,'YScale','log');

print(f,'-depsc','MeaningfulClusters-Contrasts-Distributions-Common-GC-TO-Connections-noZeros.eps');

con_label{1} = 'FC - direct';
con_label{2} = 'FC - indirect';
con_label{3} = 'GC - FROM - direct';
con_label{4} = 'GC - FROM - indirect';
con_label{5} = 'GC - TO - direct';
con_label{6} = 'GC - TO - indirect';

con_legend = {'FC - direct','FC - indirect','GC - FROM - direct','GC - FROM - indirect','GC - TO - direct','GC - TO - indirect'};

connections_Attention = zeros(1,nClusters);
connections_Stimulus = zeros(1,nClusters);
   
for iFC=2:size(MeaningFul_Attention,1)

    c(1).con = MeaningFul_Attention{iFC,7};
    c(2).con = MeaningFul_Attention{iFC,8};
    c(3).con = MeaningFul_Attention{iFC,11};
    c(4).con = MeaningFul_Attention{iFC,12};
    c(5).con = MeaningFul_Attention{iFC,15};
    c(6).con = MeaningFul_Attention{iFC,16};
    
    for icon=1:nCon
    
        connections_Attention(icon,MeaningFul_Attention{iFC,1}) = length(c(icon).con);

    end
    
end

for iFC=2:size(MeaningFul_Stimulus,1)
        
    c(1).con = MeaningFul_Stimulus{iFC,7};
    c(2).con = MeaningFul_Stimulus{iFC,8};
    c(3).con = MeaningFul_Stimulus{iFC,11};
    c(4).con = MeaningFul_Stimulus{iFC,12};
    c(5).con = MeaningFul_Stimulus{iFC,15};
    c(6).con = MeaningFul_Stimulus{iFC,16};
    
    for icon=1:nCon
    
        connections_Stimulus(icon,MeaningFul_Stimulus{iFC,1}) = length(c(icon).con);

    end
        
end

connections_Stimulus_all = sum(connections_Stimulus,1);
connections_Attention_all = sum(connections_Attention,1);

connections_Stimulus_all(connections_Stimulus_all==0) = [];
connections_Attention_all(connections_Attention_all==0) = [];

connections_Common_all = sum(FC_Common + GC_Common,1);
connections_Common_all(connections_Common_all==0) = [];

f = figure;

for icon=1:nCon
    
    subplot(3,2,icon);
    vector = connections_Attention(icon,:);
    [s,i] = sort(vector,'descend');
    vector(vector==0) = [];
    h = histfit(vector,max(vector),'gamma');
    xlabel('# of connections');
    ylabel('# of clusters');
    title(strcat(con_label(icon),':HUB:',int2str(i(1)),':#conns:',int2str(s(1))));

    hold on
    
    disp(strcat('HUB:Attention-Only:',con_label(icon),':',int2str(i(1)),':',int2str(s(1))));

end

suptitle('Increase - Attention (Only)');

print(f,'-depsc','MeaningfulClusters-Connections-Type-Distributions-Attention-Only.eps');

f = figure;

for icon=1:nCon
    
    subplot(3,2,icon);
    vector = connections_Stimulus(icon,:);
    [s,i] = sort(vector,'descend');
    vector(vector==0) = [];
    h = histfit(vector,max(vector),'gamma');
    xlabel('# of connections');
    ylabel('# of clusters');
    title(strcat(con_label(icon),':HUB:',int2str(i(1)),':#conns:',int2str(s(1))));
    hold on
    
    disp(strcat('HUB:Stimulus-Only:',con_label(icon),':',int2str(i(1)),':',int2str(s(1))));

end

suptitle('Increase - Stimulus (Only)');

print(f,'-depsc','MeaningfulClusters-Connections-Type-Distributions-Stimulus-Only.eps');

f = figure;

subplot(3,1,1);
histfit(connections_Attention_all,max(connections_Attention_all),'gamma');
xlabel('# of connections');
ylabel('# of clusters');
title('Increase - Attention (Only)');

hold on

subplot(3,1,2);
histfit(connections_Stimulus_all,max(connections_Stimulus_all),'gamma');
xlabel('# of connections');
ylabel('# of clusters');
title('Increase - Stimulus (Only)');

hold on

subplot(3,1,3);
histfit(connections_Common_all,max(connections_Common_all),'gamma');
xlabel('# of connections');
ylabel('# of clusters');
title('Increase - Common');

suptitle('All Connections');

print(f,'-depsc','MeaningfulClusters-Connections-Type-Distributions-Attention-Stimulus-Common.eps');

% figure;
% 
% subplot(1,2,1);
%     
% step = 0.001;
% max_v = 0;
% for icon=1:nCon
%     vector = connections_Attention(icon,:);
%     if max(vector) > max_v; max_v = max(vector); end
% end
% interval = 1:step:max_v;
% 
% for icon=1:nCon
%     
%     vector = connections_Attention(icon,:);
%     vector(vector==0) = [];
%     k = ksdensity(vector(:),interval);
%     probability_distribution = k .* step;
%     idx = find(probability_distribution==0);
%     probability_distribution(idx) = [];
%     interval(idx) = [];
%     plot(interval,probability_distribution,con_color{icon});
%     xlabel('# of connections');
%     ylabel('# of clusters');
%     
%     hold on
% 
% end
% 
% title('Increase - Attention (Only)');
% 
% legend(con_legend);
% 
% subplot(1,2,2);
% 
% step = 0.001;
% max_v = 0;
% for icon=1:nCon
%     vector = connections_Stimulus(icon,:);
%     if max(vector) > max_v; max_v = max(vector); end
% end
% interval = 1:step:max_v;
% 
% for icon=1:nCon
%     
%     vector = connections_Stimulus(icon,:);
%     vector(vector==0) = [];
%     k = ksdensity(vector(:),interval);
%     probability_distribution = k .* step;
%     idx = find(probability_distribution==0);
%     probability_distribution(idx) = [];
%     interval(idx) = [];
%     plot(interval,probability_distribution,con_color{icon});
%     xlabel('# of connections');
%     ylabel('# of clusters');
%     
%     hold on
%   
% end
% 
% title('Increase - Stimulus (Only)');
% 
% suptitle('Connections Distributions');

% end

