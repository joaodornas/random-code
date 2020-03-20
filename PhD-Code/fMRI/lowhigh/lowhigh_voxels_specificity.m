

function lowhigh_voxels_specificity

%% SETTINGS
disp('settings');

%settings_elena_2905;
settings_jan_0805;

%% LOAD DATA
disp('Load Data');

lowhigh_load_all_data_FSL;

%% DO THE MATH

run = 1;
doTheMath(settings,MOT4Run1,MOT2Run1,RestingStateRun1,run);
% run = 2;
% doTheMath(settings,MOT4Run2,MOT2Run2,RestingStateRun2,run);

return

function doTheMath(settings,MOT4,MOT2,RestingState,run)

%% WAVELET CLUSTERING

% disp('Frontal');
% idx_nodes_frontal = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];
% 
% for iidx=1:length(idx_nodes_frontal)
% 
%     applyWavelet(idx_nodes_frontal(iidx),MOT4,'High Attention',run,settings);
%     applyWavelet(idx_nodes_frontal(iidx),MOT2,'Low Attention',run,settings);
%     applyWavelet(idx_nodes_frontal(iidx),RestingState,'RestingState',run,settings);
% 
% end
% 
% disp('Occipital');
% idx_nodes_occipital = [43 44 49 50 51 52 53 54];
% 
% for iidx=1:length(idx_nodes_occipital)
% 
%     applyWavelet(idx_nodes_occipital(iidx),MOT4,'High Attention',run,settings);
%     applyWavelet(idx_nodes_occipital(iidx),MOT2,'Low Attention',run,settings);
%     applyWavelet(idx_nodes_occipital(iidx),RestingState,'RestingState',run,settings);
% 
% end
% 
% disp('Parietal');
% idx_nodes_parietal = [59 60 61 62];
% 
% for iidx=1:length(idx_nodes_parietal)
% 
%     applyWavelet(idx_nodes_parietal(iidx),MOT4,'High Attention',run,settings);
%     applyWavelet(idx_nodes_parietal(iidx),MOT2,'Low Attention',run,settings);
%     applyWavelet(idx_nodes_parietal(iidx),RestingState,'RestingState',run,settings);
% 
% end

% disp('Temporal');
% idx_nodes_temporal = [81 82 83 84 85 86 87 88 89 90];
% 
% for iidx=1:length(idx_nodes_temporal)
% 
%     applyWavelet(idx_nodes_temporal(iidx),MOT4,'High Attention',run,settings);
%     applyWavelet(idx_nodes_temporal(iidx),MOT2,'Low Attention',run,settings);
%     applyWavelet(idx_nodes_temporal(iidx),RestingState,'RestingState',run,settings);
% 
% end

disp('Temporal Pole');
idx_nodes_temporal_pole = [83 84 87 88];

for iidx=1:length(idx_nodes_temporal_pole)

    applyWavelet(idx_nodes_temporal_pole(iidx),MOT4,'High Attention',run,settings);
%     applyWavelet(idx_nodes_temporal_pole(iidx),MOT2,'Low Attention',run,settings);
%     applyWavelet(idx_nodes_temporal_pole(iidx),RestingState,'RestingState',run,settings);

end
return

function applyWavelet(idx_nodes,Condition,Condition_Label,run,settings)

timenow = clock;

disp(strcat(num2str(timenow(1)),'-',num2str(timenow(2)),'-',num2str(timenow(3)),'-',num2str(timenow(4)),'-',num2str(timenow(5)),'-',num2str(timenow(6))));

%% AAL
disp('AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = length(idx_nodes);
iVoxel = 0;

%% VOXELS
disp('get voxels');

nVoxels = 0;
for iNode=1:nNodes
    
    idx_ROI = AAL_ROI(idx_nodes(iNode)).ID;
    
    arealabel = AAL_ROI(idx_nodes(iNode)).Nom_L;
    arealabel = strrep(arealabel,'_','-');
    
    idx_voxels = find(AAL_img == idx_ROI);
    
    nVoxels = nVoxels + length(idx_voxels);
    
end

nTR = size(Condition,4);
area = zeros(nVoxels,nTR);

for iNode=1:nNodes
    
    idx_ROI = AAL_ROI(idx_nodes(iNode)).ID;
    
    idx_voxels = find(AAL_img == idx_ROI);
    
    for iiVoxel=1:length(idx_voxels)
        
        iVoxel = iVoxel + 1;
    
        [idxx, idxy, idxz] = ind2sub(size(AAL_img),idx_voxels(iiVoxel));
        
        area(iVoxel,:) = Condition(idxx,idxy,idxz,:);
        
    end
    
end

%% WAVELET CLUSTER
disp('Wavelet Cluster');

%k = 12;
k = 2;

lst2clu = {'s','ca1','ca3','ca6'};
S = mdwtcluster(area,'maxclust',k,'lst2clu',lst2clu);

S

IdxCLU = S.IdxCLU;
Corr = S.Corr;

maxK = max(IdxCLU(:,1));

disp(strcat('maxK:',int2str(maxK)));

%% PLOT EVERYTHING
disp('plot');

f = figure;
fs = 6;

% color(1).RGB = [0 0 255];
% color(2).RGB = [0 255 0];
% color(3).RGB = [255 0 0];
% color(4).RGB = [250 250 0];
% color(5).RGB = [250 0 250];
% color(6).RGB = [0 250 250];

col = jet;
mi = 0;
ma = maxK;
colors = 1:maxK;
RGB = squeeze(ind2rgb(floor(((colors(:)-mi)/(ma-mi))*size(col,1)),col));

for ik=1:maxK
    
    subplot(ceil(maxK/2),2,ik);

    all_clusters_voxels = area(IdxCLU(:,1)==ik,:);
    mean_all_clusters_voxels = mean(all_clusters_voxels,1);
    
%     plot(area(IdxCLU(:,1)==ik,:)','Color',RGB(ik,:));
%     
%     hold on
    
    plot(mean_all_clusters_voxels,'Color',RGB(ik,:));
    
    hold on
    
    xlabel('TRs','FontSize',fs);
    ylabel('BOLD Activity','FontSize',fs);
    set(gca,'FontSize',fs);
    
    %x = 0.7*length(mean_all_clusters_voxels);
    %y = 0.9*max(max(all_clusters_voxels));
    %text(x,y,num2str(Corr(ik)));
    
    title(strcat(Condition_Label,'-',int2str(run),'-',arealabel,'-','cluster:',int2str(ik)));

end

 print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Condition_Label,'-',int2str(run),'-',arealabel,'-','wavelet-cluster-',int2str(k),'.jpg'));
 print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Condition_Label,'-',int2str(run),'-',arealabel,'-','wavelet-cluster-',int2str(k),'.eps'));
 print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Condition_Label,'-',int2str(run),'-',arealabel,'-','wavelet-cluster-',int2str(k),'.pdf'));
 
 close all

 volume = zeros(size(AAL_img));
 
 iVoxel = 0;
 for iNode=1:nNodes
    
    idx_ROI = AAL_ROI(idx_nodes(iNode)).ID;
    
    idx_voxels = find(AAL_img == idx_ROI);
    
    for iiVoxel=1:length(idx_voxels)
        
        iVoxel = iVoxel + 1;
    
        [idxx, idxy, idxz] = ind2sub(size(AAL_img),idx_voxels(iiVoxel));
        
        ik = IdxCLU(iVoxel,1);
        
        volume(idxx,idxy,idxz) = ik;
        
    end
    
 end

file_label = strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Condition_Label,'-',int2str(run),'-',arealabel,'-','wavelet-cluster-',int2str(k));

fname = strcat(file_label,'.nii');
scl_slope = 1;
scl_inter = 0; 
dim = [size(AAL_img,1),size(AAL_img,2),size(AAL_img,3)];
dtype = 'FLOAT32';
offset = 0;
descrip = 'WAVELET CLUSTER';
nifti_file.mat = load_aal.mat;
nifti_file.mat_intent = load_aal.mat_intent;
nifti_file.mat0 = load_aal.mat0;
nifti_file.mat0_intent = load_aal.mat0_intent;

input_data = volume; 

lowhigh_save_image;
 
rendfile = 'K:\Dropbox (Uni Magdeburg)\_TOOLBOX\spm8\canonical\cortex_20484.surf.gii';
img_dat = fromImageToSPMdat(fname);

brt = 1;
min_value = min(img_dat.t);
max_value = max(img_dat.t);

setGlobalMinMax(min_value,max_value);

snapshot_label = file_label;
setRenderingSnapshot(true,snapshot_label);

spm_render_my(img_dat,brt,rendfile,fname);
 
return