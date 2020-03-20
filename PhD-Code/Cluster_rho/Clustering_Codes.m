%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hierarchical Clustering of Correlation Maps

EPI_Ind = rho;

% Correlation matrix
CorMx = corrcoef(EPI_Ind);
% Distances between correlations profiles
temp = corrcoef(CorMx);
temp(temp==1) = 0;
Dist = 1-squareform(temp);
% built the hierarhical cluster tree
Link = linkage(Dist,'average');

%% Breaking down the tree according to a distance threshold
% the plot of # of clusters vs. the threshold
clstnum=0;
for dist_thr = 0.05:.05:2;
    T_bd = cluster(Link,'cutoff', dist_thr,'criterion','distance');
    clear num_avg;
    for i=1:max(T_bd)
       num_avg(i) = length(find(T_bd==i));
    end
    clstnum(end+1,1) = max(T_bd);
    clstnum(end,2) = sum(num_avg>=8);
end
clstnum(1,:)=[];
figure,plot(0.05:.05:2,(clstnum));ylim([0 100]);xlim([0 2]);
% final choice
dist_thr = 0.4;
T_fc = cluster(Link,'cutoff', dist_thr,'criterion','distance');

%% Breaking down the tree according to a threshold for inconsistent coefficients
% inc_dep = 2;
% Y = inconsistent(Link,inc_dep);
% % figure('color','w'); hist(Y(:,4),500);
% % h = findobj(gca,'Type','patch');set(h,'facecolor','k','edgecolor','k'); xlim([0.8 1.2]);
% % set(gca,'fontname','Arial','fontsize',14,'linewidth',1.5);
% inconsist_thr = prctile(Y(:,4),99); 
%  
% % 
% T = cluster(Link,'cutoff',inconsist_thr,'depth',inc_dep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% K-means Clustering of Correlation Maps
num_Clust = 8;
[T_km,C,sumd,D] = kmeans(CorMx, num_Clust , 'distance', 'correlation', 'display', 'final','replicate',20);


rho_cluster = zeros(size(rho));
k = 0;
for icluster=1:num_Clust
    
    idx = find(T_km == icluster);
    
    for iidx=1:length(idx)
    
        k = k + 1;
        
        column(k) = idx(iidx);
        new_column(k) = k;
        
    end
    
end

rho_cluster(:,new_column) = rho(:,column);
rho_cluster(new_column,:) = rho(column,:);
    

















