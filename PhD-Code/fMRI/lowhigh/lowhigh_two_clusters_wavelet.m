
settings_jan_0805;
%settings_elena_2905;

file = settings.FSL.files.functional.warped;

lowhigh_load_all_data_FSL;


load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

nTR = size(MOT4Run1,4);

for iNode=1:nNodes
    
    iNode
    
    idx_ROI = AAL_ROI(iNode).ID;
    
    idx_voxels = find(AAL_img == idx_ROI);
    
    nVoxels = length(idx_voxels);
    
    area_MOT4Run1 = zeros(nVoxels,nTR);
    
    for iVoxel=1:nVoxels
        
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

        ROI(iNode).area_MOT4Run1(iVoxel,:) = MOT4Run1(idxx,idxy,idxz,:);
        ROI(iNode).nVoxels = nVoxels;
        
    end
        
    k = 2;
    lst2clu = {'s','ca1','ca3','ca6'};

    ROI(iNode).area_MOT4Run1_S = mdwtcluster(ROI(iNode).area_MOT4Run1,'maxclust',k,'lst2clu',lst2clu);

end

f = figure;

l = 0;

fs = 5;

for iNode=1:20

    label = AAL_ROI(iNode).Nom_L;
    label = strrep(label,'_','-');
    
    l = l + 1;
    
    subplot(10,4,l);
    
    plot(mean(ROI(iNode).area_MOT4Run1(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==1,:),1));
    
    nVoxels = length(find(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==1));
    
    title(strcat(label,':',num2str(nVoxels),'-','cluster-1'),'FontSize',fs);
    xlabel('TRs','FontSize',fs);
    ylabel('BOLD activity','FontSize',fs);
    xlim([0 nTR]);
    
    set(gca,'FontSize',fs);
    
    hold on
    
    l = l + 1;
    
    subplot(10,4,l);
    
    plot(mean(ROI(iNode).area_MOT4Run1(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==2,:),1));
    
    nVoxels = length(find(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==2));
    
    title(strcat(label,':',num2str(nVoxels),'-','cluster-2'),'FontSize',fs);
    xlabel('TRs','FontSize',fs);
    ylabel('BOLD activity','FontSize',fs);
    xlim([0 nTR]);
    
    set(gca,'FontSize',fs);
    
end

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-1-20','.jpeg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-1-20','.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-1-20','.pdf'));

close all

f = figure;

l = 0;

for iNode=21:40

    label = AAL_ROI(iNode).Nom_L;
    label = strrep(label,'_','-');
    
    l = l + 1;
    
    subplot(10,4,l);
    
    plot(mean(ROI(iNode).area_MOT4Run1(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==1,:),1));
    
    nVoxels = length(find(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==1));
    
    title(strcat(label,':',num2str(nVoxels),'-','cluster-1'),'FontSize',fs);
    xlabel('TRs','FontSize',fs);
    ylabel('BOLD activity','FontSize',fs);
    xlim([0 nTR]);
    
    set(gca,'FontSize',fs);
    
    hold on
    
    l = l + 1;
    
    subplot(10,4,l);
    
    plot(mean(ROI(iNode).area_MOT4Run1(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==2,:),1));
    
    nVoxels = length(find(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==2));
    
    title(strcat(label,':',num2str(nVoxels),'-','cluster-2'),'FontSize',fs);
    xlabel('TRs','FontSize',fs);
    ylabel('BOLD activity','FontSize',fs);
    xlim([0 nTR]);
    
    set(gca,'FontSize',fs);
    
end

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-21-40','.jpeg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-21-40','.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-21-40','.pdf'));

close all

f = figure;

l = 0;

for iNode=41:60

    label = AAL_ROI(iNode).Nom_L;
    label = strrep(label,'_','-');
    
    l = l + 1;
    
    subplot(10,4,l);
    
    plot(mean(ROI(iNode).area_MOT4Run1(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==1,:),1));
    
    nVoxels = length(find(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==1));
    
    title(strcat(label,':',num2str(nVoxels),'-','cluster-1'),'FontSize',fs);
    xlabel('TRs','FontSize',fs);
    ylabel('BOLD activity','FontSize',fs);
    xlim([0 nTR]);
    
    set(gca,'FontSize',fs);
    
    hold on
    
    l = l + 1;
    
    subplot(10,4,l);
    
    plot(mean(ROI(iNode).area_MOT4Run1(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==2,:),1));
    
    nVoxels = length(find(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==2));
    
    title(strcat(label,':',num2str(nVoxels),'-','cluster-2'),'FontSize',fs);
    xlabel('TRs','FontSize',fs);
    ylabel('BOLD activity','FontSize',fs);
    xlim([0 nTR]);
    
    set(gca,'FontSize',fs);
    
end

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-41-60','.jpeg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-41-60','.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-41-60','.pdf'));

f = figure;

l = 0;

for iNode=61:80

    label = AAL_ROI(iNode).Nom_L;
    label = strrep(label,'_','-');
    
    l = l + 1;
    
    subplot(10,4,l);
    
    plot(mean(ROI(iNode).area_MOT4Run1(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==1,:),1));
    
    nVoxels = length(find(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==1));
    
    title(strcat(label,':',num2str(nVoxels),'-','cluster-1'),'FontSize',fs);
    xlabel('TRs','FontSize',fs);
    ylabel('BOLD activity','FontSize',fs);
    xlim([0 nTR]);
    
    set(gca,'FontSize',fs);
    
    hold on
    
    l = l + 1;
    
    subplot(10,4,l);
    
    plot(mean(ROI(iNode).area_MOT4Run1(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==2,:),1));
    
    nVoxels = length(find(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==2));
    
    title(strcat(label,':',num2str(nVoxels),'-','cluster-2'),'FontSize',fs);
    xlabel('TRs','FontSize',fs);
    ylabel('BOLD activity','FontSize',fs);
    xlim([0 nTR]);
    
    set(gca,'FontSize',fs);
    
end

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-61-80','.jpeg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-61-80','.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-61-80','.pdf'));

f = figure;

l = 0;

for iNode=81:90

    label = AAL_ROI(iNode).Nom_L;
    label = strrep(label,'_','-');
    
    l = l + 1;
    
    subplot(10,4,l);
    
    plot(mean(ROI(iNode).area_MOT4Run1(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==1,:),1));
    
    nVoxels = length(find(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==1));
    
    title(strcat(label,':',num2str(nVoxels),'-','cluster-1'),'FontSize',fs);
    xlabel('TRs','FontSize',fs);
    ylabel('BOLD activity','FontSize',fs);
    xlim([0 nTR]);
    
    set(gca,'FontSize',fs);
    
    hold on
    
    l = l + 1;
    
    subplot(10,4,l);
    
    plot(mean(ROI(iNode).area_MOT4Run1(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==2,:),1));
    
    nVoxels = length(find(ROI(iNode).area_MOT4Run1_S.IdxCLU(:,1)==2));
    
    title(strcat(label,':',num2str(nVoxels),'-','cluster-2'),'FontSize',fs);
    xlabel('TRs','FontSize',fs);
    ylabel('BOLD activity','FontSize',fs);
    xlim([0 nTR]);
    
    set(gca,'FontSize',fs);
    
end

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-81-90','.jpeg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-81-90','.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','MOT4Run1-Wavelet-2-Clusters-81-90','.pdf'));
