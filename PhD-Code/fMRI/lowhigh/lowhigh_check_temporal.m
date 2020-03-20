
%% LOAD SETTINGS

settings_jan_0805;
%settings_elena_2905;

%%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.standard;
mask = settings.FSL.files.mask.custom;

lowhigh_load_all_data_FSL;

%% PARAMETERS

k = 4;

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;
nTR = size(MOT4Run1,4);

colors{1} = 'b';
colors{2} = 'r';
colors{3} = 'y';

idx_Temporal_Poles = [83 84 87 88];

for iNode=1:length(idx_Temporal_Poles)
   %for iNode=3

    idx_region = AAL_ROI(idx_Temporal_Poles(iNode)).ID;
    label = AAL_ROI(idx_Temporal_Poles(iNode)).Nom_L;
    label = strrep(label,'_','-');

    idx_voxels = find(AAL_img == idx_region);
    nVoxels = length(idx_voxels);
    
    area_MOT4Run1 = zeros(nVoxels,nTR);
    area_MOT2Run1 = zeros(nVoxels,nTR);
    area_RestingStateRun1 = zeros(nVoxels,nTR);
    
    for iVoxel=1:nVoxels
        
        [idxx, idxy, idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
        
        area_MOT4Run1(iVoxel,:) = MOT4Run1(idxx,idxy,idxz,:);
        area_MOT2Run1(iVoxel,:) = MOT2Run1(idxx,idxy,idxz,:);
        area_RestingStateRun1(iVoxel,:) = RestingStateRun1(idxx,idxy,idxz,:);
        
    end
    
    %% PLOT MEAN AMONG VOXELS
    
    mean_voxels_MOT4Run1 = mean(area_MOT4Run1,1);
    mean_voxels_MOT2Run1 = mean(area_MOT2Run1,1);
    mean_voxels_RestingStateRun1 = mean(area_RestingStateRun1,1);
    
    f = figure;
    
    plot(mean_voxels_MOT4Run1,colors{1});
    hold on
    plot(mean_voxels_MOT2Run1,colors{2});
    plot(mean_voxels_RestingStateRun1,colors{3});
    
    print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',label,'.jpeg'));
    print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',label,'.eps'));
    print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',label,'.pdf'));


    %% PLOT HIGHEST MEAN VOXEL

%     mean_time_MOT4Run1 = mean(area_MOT4Run1,2);
%     mean_time_MOT2Run1 = mean(area_MOT2Run1,2);
%     mean_time_RestingStateRun1 = mean(area_RestingStateRun1,2);
%     
%     [x, I] = sort(mean_time_MOT4Run1);
%     idx_hm_voxel_MOT4Run1 = I(end);
%     
%     [x, I] = sort(mean_time_MOT2Run1);
%     idx_hm_voxel_MOT2Run1 = I(end);
%     
%     [x, I] = sort(mean_time_RestingStateRun1);
%     idx_hm_voxel_RestingStateRun1 = I(end);
%     
%     hm_voxel_MOT4Run1 = squeeze(area_MOT4Run1(idx_hm_voxel_MOT4Run1,:));
%     hm_voxel_MOT2Run1 = squeeze(area_MOT2Run1(idx_hm_voxel_MOT2Run1,:));
%     hm_voxel_RestingStateRun1 = squeeze(area_RestingStateRun1(idx_hm_voxel_RestingStateRun1,:));
%     
%     plot(hm_voxel_MOT4Run1,colors{1});
%     hold on
%     plot(hm_voxel_MOT2Run1,colors{2});
%     plot(hm_voxel_RestingStateRun1,colors{3});    


    %% PLOT ALL VOXELS
    
%     for iVoxel=1:nVoxels
%         
%        %plot(area_MOT4Run1(iVoxel,:),colors{1});
%         hold on
%         plot(area_MOT2Run1(iVoxel,:),colors{2});
%         plot(area_RestingStateRun1(iVoxel,:),colors{3});
%     
%     end

    %% PLOT CLUSTERS

%     k = 4;
%     fs = 6;
% 
%     S_MOT4Run1 = mdwtcluster(area_MOT4Run1,'maxclust',k);
%     S_MOT2Run1 = mdwtcluster(area_MOT2Run1,'maxclust',k);
%     S_RestingStateRun1 = mdwtcluster(area_RestingStateRun1,'maxclust',k);
%     
%     MOT4Run1IdxCLU = S_MOT4Run1.IdxCLU;
%     MOT2Run1IdxCLU = S_MOT2Run1.IdxCLU;
%     RestingStateRun1IdxCLU = S_RestingStateRun1.IdxCLU;
%     
%     f = figure;
%     nClusters = k;   
%     j = 0;
%     for i=1:nClusters
%    
%         j = j + 1;
%         subplot(k,3,j);
%         plot(area_MOT4Run1(MOT4Run1IdxCLU(:,1)==i,:)',colors{1});
%         
%         title(strcat(label,'- Cluster:',int2str(i)),'FontSize',fs);
%         ylabel('BOLD activity','FontSize',fs);
%         xlabel('TRs','FontSize',fs);
%         
%         set(gca,'FontSize',fs);
%         
%         hold on
%         
%         j = j + 1;
%         subplot(k,3,j);
%         plot(area_MOT2Run1(MOT2Run1IdxCLU(:,1)==i,:)',colors{2});
%         
%         title(strcat(label,'- Cluster:',int2str(i)),'FontSize',fs);
%         ylabel('BOLD activity','FontSize',fs);
%         xlabel('TRs','FontSize',fs);
%         
%         set(gca,'FontSize',fs);
%         
%         hold on
%         
%         j = j + 1;
%         subplot(k,3,j);
%         plot(area_RestingStateRun1(RestingStateRun1IdxCLU(:,1)==i,:)',colors{3});
%         
%         title(strcat(label,'- Cluster:',int2str(i)),'FontSize',fs);
%         ylabel('BOLD activity','FontSize',fs);
%         xlabel('TRs','FontSize',fs);
%         
%         set(gca,'FontSize',fs);
%         
%         hold on
%         
%     end
%     
%     print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',label,'-wavelet-cluster','.jpeg'));
%     print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',label,'-wavelet-cluster','.eps'));
%     print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',label,'-wavelet-cluster','.pdf'));
%         
    
end
    