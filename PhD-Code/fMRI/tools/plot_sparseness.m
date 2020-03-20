
clear all

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS.mat')

iiCluster = 0;

for iROI=1:90

    nClusters = length(ROI(iROI).clusters);

    for iClusters=1:nClusters

        iiCluster = iiCluster + 1; m_track(iiCluster,:) = squeeze(mean(ROI(iROI).clusters(iClusters).track_mean,1));m_passive(iiCluster,:) = squeeze(mean(ROI(iROI).clusters(iClusters).passive_mean,1));m_rest(iiCluster,:)=squeeze(mean(ROI(iROI).clusters(iClusters).rest_mean,1));

    end
    
end

% figure;
% 
% max_y = max(m_track(1,:));
% min_y = min(m_track(1,:));
% 
% interval = linspace(min_y,max_y,100);
% 
% [N, X] = hist(m_track(1,:),interval);
% 
% idx_zeros = find(N==0);
% N(idx_zeros) = [];
% X(idx_zeros) = [];
% 
% plot(X,N,'b-');
% ylim([0.5 6]);

% figure;
% 
% max_y = max(m_passive(1,:));
% min_y = min(m_passive(1,:));
% 
% interval = linspace(min_y,max_y,1000);
% 
% hist(m_passive(1,:),interval);
% 
% figure;
% 
% max_y = max(m_rest(1,:));
% min_y = min(m_rest(1,:));
% 
% interval = linspace(min_y,max_y,1000);
% 
% hist(m_rest(1,:),interval);

for i=1:3
    
    figure;
    
    max_y = max(m_track(i,:));
    min_y = min(m_track(i,:));

    interval = linspace(min_y,max_y,100);

    [N_track, X] = hist(m_track(i,:),interval);

    idx_zeros = find(N_track==0);
    N_track(idx_zeros) = [];
    X(idx_zeros) = [];

    plot(X,N_track,'b-');
    hold on
    
    max_y = max(m_passive(i,:));
    min_y = min(m_passive(i,:));

    interval = linspace(min_y,max_y,100);

    [N_passive, X] = hist(m_passive(i,:),interval);

    idx_zeros = find(N_passive==0);
    N_passive(idx_zeros) = [];
    X(idx_zeros) = [];

    plot(X,N_passive,'k-');
    hold on
    
    max_y = max(m_rest(i,:));
    min_y = min(m_rest(i,:));

    interval = linspace(min_y,max_y,100);

    [N_rest, X] = hist(m_rest(i,:),interval);

    idx_zeros = find(N_rest==0);
    N_rest(idx_zeros) = [];
    X(idx_zeros) = [];

    plot(X,N_rest,'r-');
    hold on
    
    ylim([0.5 max([N_track, N_passive, N_rest])+0.5]);

end

