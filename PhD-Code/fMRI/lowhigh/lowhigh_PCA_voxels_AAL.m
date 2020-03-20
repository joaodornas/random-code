
function lowhigh_PCA_voxels_AAL

%% SETTINGS

settings_jan_0805;

doTheMath(settings);

clear settings
settings_elena_2905;

doTheMath(settings);

return

function doTheMath(settings)
%% LOAD DATA

lowhigh_load_all_data_FSL;

%% AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

%% AREAS

disp('Frontal');
idx_nodes = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];
doPCA(MOT4Run1,MOT2Run1,RestingStateRun1,AAL_img,AAL_ROI,settings,idx_nodes,1);
doPCA(MOT4Run2,MOT2Run2,RestingStateRun2,AAL_img,AAL_ROI,settings,idx_nodes,2);

disp('Occipital');
idx_nodes = [49 50 51 52 53 54];
doPCA(MOT4Run1,MOT2Run1,RestingStateRun1,AAL_img,AAL_ROI,settings,idx_nodes,1);
doPCA(MOT4Run2,MOT2Run2,RestingStateRun2,AAL_img,AAL_ROI,settings,idx_nodes,2);

disp('Parietal');
idx_nodes = [59 60 61 62];
doPCA(MOT4Run1,MOT2Run1,RestingStateRun1,AAL_img,AAL_ROI,settings,idx_nodes,1);
doPCA(MOT4Run2,MOT2Run2,RestingStateRun2,AAL_img,AAL_ROI,settings,idx_nodes,2);

disp('Temporal');
idx_nodes = [81 82 83 84 85 86 87 88 89 90];
doPCA(MOT4Run1,MOT2Run1,RestingStateRun1,AAL_img,AAL_ROI,settings,idx_nodes,1);
doPCA(MOT4Run2,MOT2Run2,RestingStateRun2,AAL_img,AAL_ROI,settings,idx_nodes,2);

return

function doPCA(MOT4,MOT2,RestingState,AAL_img,AAL_ROI,settings,idx_nodes,run)

% MOT4 = MOT4Run1;
% MOT2 = MOT2Run1;
% RestingState = RestingStateRun1;

nTR = size(MOT4,4);

disp('process voxels');
for iNode=1:length(idx_nodes)
   
    idx = idx_nodes(iNode);
    
    idx_AAL = AAL_ROI(idx).ID;
    
    label_AAL = AAL_ROI(idx).Nom_L;
    
    idx_voxels = find(AAL_img==idx_AAL);
    
    nVoxels = length(idx_voxels);
    
    disp(strcat(label_AAL,':',num2str(nVoxels),'(voxels)'));
    
    for iVoxel=1:nVoxels
        
        idx = idx_voxels(iVoxel);
        
        [idxx(iVoxel), idxy(iVoxel), idxz(iVoxel)] = ind2sub(size(AAL_img),idx);
        
        area_MOT4(iVoxel,:) = MOT4(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        area_MOT2(iVoxel,:) = MOT2(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);
        area_RestingState(iVoxel,:) = RestingState(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),:);

    end
    
    MOT4_mean = mean(area_MOT4,2);
    MOT2_mean = mean(area_MOT2,2);
    RestingState_mean = mean(area_RestingState,2);
        
    [m, I] = sort(MOT4_mean);
    MOT4_idx_voxel = I(end);
    
    [m, I] = sort(MOT2_mean);
    MOT2_idx_voxel = I(end);
    
    [m, I] = sort(RestingState_mean);
    RestingState_idx_voxel = I(end);
    
    disp(strcat('idx_voxels:',int2str(MOT4_idx_voxel),':',int2str(MOT2_idx_voxel),':',int2str(RestingState_idx_voxel)));
    
    MOT4_mean_time_series = area_MOT4(MOT4_idx_voxel,:);
    MOT2_mean_time_series = area_MOT2(MOT2_idx_voxel,:);
    RestingState_mean_time_series = area_RestingState(RestingState_idx_voxel,:);
    
    mn_MOT4 = mean( area_MOT4, 2);                % compute row mean
    X0_area_MOT4 = area_MOT4 - repmat(mn_MOT4,1,nTR);         % subtract row mean
    
    mn_MOT2 = mean( area_MOT2, 2);                % compute row mean
    X0_area_MOT2 = area_MOT2 - repmat(mn_MOT2,1,nTR);         % subtract row mean
    
    mn_RestingState = mean( area_RestingState, 2);                % compute row mean
    X0_area_RestingState = area_RestingState - repmat(mn_RestingState,1,nTR);         % subtract row mean
    
    [U_area_MOT4, S_area_MOT4, V_area_MOT4] = svd( X0_area_MOT4' ); 
    [U_area_MOT2, S_area_MOT2, V_area_MOT2] = svd( X0_area_MOT2' ); 
    [U_area_RestingState, S_area_RestingState, V_area_RestingState] = svd( X0_area_RestingState' ); 
    
    V_area_MOT4 = V_area_MOT4';
    V_area_MOT2 = V_area_MOT2';
    V_area_RestingState =  V_area_RestingState';
    
    Y0_area_MOT4 = V_area_MOT4*X0_area_MOT4;
    Y0_area_MOT2 = V_area_MOT2*X0_area_MOT2;
    Y0_area_RestingState = V_area_RestingState*X0_area_RestingState;
    
    Y0dn_area_MOT4 = Y0_area_MOT4;
    Y0dn_area_MOT4(1,1:nTR) = zeros(1,nTR);                 
    X0dn_area_MOT4 = V_area_MOT4' * Y0dn_area_MOT4;                                         
    Xdn_area_MOT4  = X0dn_area_MOT4 + repmat(mn_MOT4,1,nTR);
    
    Y0dn_area_MOT2 = Y0_area_MOT2;
    Y0dn_area_MOT2(1,1:nTR) = zeros(1,nTR);                 
    X0dn_area_MOT2 = V_area_MOT2' * Y0dn_area_MOT2;                                         
    Xdn_area_MOT2  = X0dn_area_MOT2 + repmat(mn_MOT2,1,nTR);
    
    Y0dn_area_RestingState = Y0_area_RestingState;
    Y0dn_area_RestingState(1,1:nTR) = zeros(1,nTR);                 
    X0dn_area_RestingState = V_area_RestingState' * Y0dn_area_RestingState;                                         
    Xdn_area_RestingState  = X0dn_area_RestingState + repmat(mn_RestingState,1,nTR);
    
    MOT4dn_mean_time_series = Xdn_area_MOT4(MOT4_idx_voxel,:);
    MOT2dn_mean_time_series = Xdn_area_MOT2(MOT2_idx_voxel,:);
    RestingStatedn_mean_time_series = Xdn_area_RestingState(RestingState_idx_voxel,:);
    
    diag_S_area_MOT4 = diag(S_area_MOT4);
    diag_S_area_MOT2 = diag(S_area_MOT2);
    diag_S_area_RestingState = diag(S_area_RestingState);
    
    variance1st_area_MOT4 = 100 * ( diag_S_area_MOT4(1) / sum(diag_S_area_MOT4));
    variance1st_area_MOT2 = 100 * ( diag_S_area_MOT2(1) / sum(diag_S_area_MOT2));
    variance1st_area_RestingState = 100 * ( diag_S_area_RestingState(1) / sum(diag_S_area_RestingState));
    
    n95th_area_MOT4 = get95thvariance(diag_S_area_MOT4);
    n95th_area_MOT2 = get95thvariance(diag_S_area_MOT2);
    n95th_area_RestingState = get95thvariance(diag_S_area_RestingState);
    
    %plotComponents(nVoxels,S_area_MOT4,Y0_area_MOT4,variance1st_area_MOT4,n95th_area_MOT4,settings,'MOT4',label_AAL,run);
    %plotComponents(nVoxels,S_area_MOT2,Y0_area_MOT2,variance1st_area_MOT2,n95th_area_MOT2,settings,'MOT2',label_AAL,run);
    %plotComponents(nVoxels,S_area_RestingState,Y0_area_RestingState,variance1st_area_RestingState,n95th_area_RestingState,settings,'RestingState',label_AAL,run);
    
    plotAllVariances(nVoxels,S_area_MOT4,S_area_MOT2,S_area_RestingState,variance1st_area_MOT4,variance1st_area_MOT2,variance1st_area_RestingState,n95th_area_MOT4,n95th_area_MOT2,n95th_area_RestingState,settings,label_AAL,run);
    
    plotRemovedComponents(Y0_area_MOT4,Y0_area_MOT2,Y0_area_RestingState,MOT4_mean_time_series,MOT2_mean_time_series,RestingState_mean_time_series,MOT4dn_mean_time_series,MOT2dn_mean_time_series,RestingStatedn_mean_time_series,settings,label_AAL,run);
    
    clear area_MOT4
    clear area_MOT2
    clear area_RestingState
    
end

return

function plotRemovedComponents(Y0_area_MOT4,Y0_area_MOT2,Y0_area_RestingState,MOT4_mean_time_series,MOT2_mean_time_series,RestingState_mean_time_series,MOT4dn_mean_time_series,MOT2dn_mean_time_series,RestingStatedn_mean_time_series,settings,label_AAL,run)
    

    f = figure;
    
    subplot(1,3,1);
    
    plot(Y0_area_MOT4(1,:),'b');
    hold on
    plot(Y0_area_MOT2(1,:),'r');
    plot(Y0_area_RestingState(1,:),'y');
    
    legend({'High Attention', 'Low Attention', 'Resting State'});
    
    xlabel('TR');
    ylabel('BOLD Activity');
    
    title('1st PCs');
    
    subplot(1,3,2);
    
    plot(MOT4_mean_time_series,'b');
    hold on
    plot(MOT2_mean_time_series,'r');
    plot(RestingState_mean_time_series,'y');
    
    legend({'High Attention', 'Low Attention', 'Resting State'});
    
    xlabel('TR');
    ylabel('BOLD Activity');
    
    title('HMVoxel');
    
    ylim([min(min([MOT4_mean_time_series, MOT2_mean_time_series, RestingState_mean_time_series])) max(max([MOT4_mean_time_series, MOT2_mean_time_series, RestingState_mean_time_series]))]);
    
    subplot(1,3,3);
    
    plot(MOT4dn_mean_time_series,'b');
    hold on
    plot(MOT2dn_mean_time_series,'r');
    plot(RestingStatedn_mean_time_series,'y');
    
    legend({'High Attention', 'Low Attention', 'Resting State'});
    
    xlabel('TR');
    ylabel('BOLD Activity');
    
    title('HMVoxel-1stPC');
    
    ylim([min(min([MOT4dn_mean_time_series, MOT2dn_mean_time_series, RestingStatedn_mean_time_series])) max(max([MOT4dn_mean_time_series, MOT2dn_mean_time_series, RestingStatedn_mean_time_series]))]);
    
    
    ax=axes('Units','Normal','Position',[.075 .075 .90 .90],'Visible','off');
    set(get(ax,'Title'),'Visible','on');
    
    label_AAL = strrep(label_AAL,'_','-');
    title(strcat(label_AAL,'-Run-',int2str(run)));
    
    print(f,'-djpeg',strcat('Low-High-',settings.folders.subject,'-','All','-','PCA-PCs','-',label_AAL,'-Run-',int2str(run),'.jpeg'));
    print(f,'-depsc',strcat('Low-High-',settings.folders.subject,'-','All','-','PCA-PCs','-',label_AAL,'-Run-',int2str(run),'.eps'));
    
    close all


return

function plotAllVariances(nVoxels,S_area_MOT4,S_area_MOT2,S_area_RestingState,varMOT4, varMOT2, varRestingState, n95thMOT4, n95thMOT2, n95thRestingState, settings,label_AAL,run)

    fs = 12;

    g = figure;
 
    for i=1:3
        
        subplot(1,3,i);

        if i==1, label = 'High Attention'; S = S_area_MOT4; variance1st = varMOT4; n95th = n95thMOT4; end
        if i==2, label = 'Low Attention'; S = S_area_MOT2; variance1st = varMOT2; n95th = n95thMOT2; end
        if i==3, label = 'Resting State'; S = S_area_RestingState; variance1st = varRestingState; n95th = n95thRestingState; end
        
        hold on;
        bar( log( diag(S) ) );              % plot variance of principal components
        hold off;
        title( label );
        
        x = length(log(diag(S)))/4;
        y = min(log(diag(S))) + (max(log(diag(S))) - min(log(diag(S))))/2;
        text(x,y,strcat('#voxels:',num2str(nVoxels)),'FontSize',fs);
        
        y = min(log(diag(S))) + (max(log(diag(S))) - min(log(diag(S))))/2 + 1;
        x = length(log(diag(S)))/4;
        text(x,y,strcat('1st PC:',num2str(variance1st),'%'),'FontSize',fs);
    
        y = min(log(diag(S))) + (max(log(diag(S))) - min(log(diag(S))))/2 + 2;
        x = length(log(diag(S)))/4;
        text(x,y,strcat('#PCs (95%):',num2str(n95th)),'FontSize',fs);
        
        xlabel( 'PCs');
        ylabel( 'log variance' );
    
        xlim([0 length(log(diag(S)))]);
        ylim([min(log(diag(S))) max(log(diag(S)))]);
        
    end
    
    ax=axes('Units','Normal','Position',[.075 .075 .90 .90],'Visible','off');
    set(get(ax,'Title'),'Visible','on');
    
    label_AAL = strrep(label_AAL,'_','-');
    title(strcat(label_AAL,'-Run-',int2str(run)));
    
    print(g,'-djpeg',strcat('Low-High-',settings.folders.subject,'-','All','-','PCA-variance','-',label_AAL,'-Run-',int2str(run),'.jpeg'));
    print(g,'-depsc',strcat('Low-High-',settings.folders.subject,'-','All','-','PCA-variance','-',label_AAL,'-Run-',int2str(run),'.eps'));

return

function plotComponents(nVoxels,S,components,variance1st,n95th,settings,kind,label_AAL,run)

    fs = 12;

    g = figure;
    hold on;
    bar( log( diag(S) ) );              % plot variance of principal components
    hold off;
    title( 'log variance explained by PCs' );
    x = length(log(diag(S)))/2;
    y = min(log(diag(S))) + (max(log(diag(S))) - min(log(diag(S))))/2;
    text(x,y,strcat('# voxels:',num2str(nVoxels)),'FontSize',fs);
    xlabel( 'PCs');
    ylabel( 'log variance' );
    
    print(g,'-djpeg',strcat('Low-High-',settings.folders.subject,'-',kind,'-','PCA-variance','-',label_AAL,'-Run-',int2str(run),'.jpeg'));
    print(g,'-depsc',strcat('Low-High-',settings.folders.subject,'-',kind,'-','PCA-variance','-',label_AAL,'-Run-',int2str(run),'.eps'));
    
    f = figure;
   
    max_y = max([max(components(1,:)), max(components(2,:)), max(components(3,:)), max(components(4,:))]);
    min_y = min([min(components(1,:)), min(components(2,:)), min(components(3,:)), min(components(4,:))]);
    
    plot(components(1,:),'b');
    hold on
    plot(components(2,:),'r');
    plot(components(3,:),'g');
    plot(components(4,:),'k');
    
    legend({'1st','2nd','3rd','4th'});
    
    label_AAL = strrep(label_AAL,'_','-');
    
    title(label_AAL);
    xlabel('TR');
    ylabel('BOLD Activity');
    
    y = min_y + (max_y - min_y)/2;
    x = length(components(1,:))/2;
    
    text(x,y,strcat('1st component variance:',num2str(variance1st),'%'),'FontSize',fs);
    
    y = min_y + (max_y - min_y)/4;
    x = length(components(1,:))/2;
    
    text(x,y,strcat('# components 95% variance:',num2str(n95th)),'FontSize',fs);

    print(f,'-djpeg',strcat('Low-High-',settings.folders.subject,'-',kind,'-','PCA-components','-',label_AAL,'-Run-',int2str(run),'.jpeg'));
    print(f,'-depsc',strcat('Low-High-',settings.folders.subject,'-',kind,'-','PCA-components','-',label_AAL,'-Run-',int2str(run),'.eps'));
    
    close all

return

function iVar = get95thvariance(diagvar)

reached95 = false;
variance = 0;
iVar = 0;

while ~reached95
   
    iVar = iVar + 1;
    
    variance = variance + 100 * ( diagvar(iVar) / sum(diagvar) );
    
    if variance > 95
    
        reached95 = true;
        
    end
    
end

return
