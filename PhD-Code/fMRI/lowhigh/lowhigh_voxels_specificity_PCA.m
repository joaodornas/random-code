

function lowhigh_voxels_specificity_PCA

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
run = 2;
doTheMath(settings,MOT4Run2,MOT2Run2,RestingStateRun2,run);

return

function doTheMath(settings,MOT4,MOT2,RestingState,run)

%% PCA on each area and condition

disp('Frontal');
idx_nodes_frontal = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 24 25 26];
 for iidx=1:length(idx_nodes_frontal)
 
     applyPCA(idx_nodes_frontal(iidx),MOT4,'High Attention',run,settings);
     applyPCA(idx_nodes_frontal(iidx),MOT2,'Low Attention',run,settings);
     applyPCA(idx_nodes_frontal(iidx),RestingState,'RestingState',run,settings);
 
 end

disp('Occipital');
idx_nodes_occipital = [43 44 49 50 51 52 53 54];

for iidx=1:length(idx_nodes_occipital)

    applyPCA(idx_nodes_occipital(iidx),MOT4,'High Attention',run,settings);
    applyPCA(idx_nodes_occipital(iidx),MOT2,'Low Attention',run,settings);
    applyPCA(idx_nodes_occipital(iidx),RestingState,'RestingState',run,settings);

end

disp('Parietal');
idx_nodes_parietal = [59 60 61 62];

for iidx=1:length(idx_nodes_parietal)

    applyPCA(idx_nodes_parietal(iidx),MOT4,'High Attention',run,settings);
    applyPCA(idx_nodes_parietal(iidx),MOT2,'Low Attention',run,settings);
    applyPCA(idx_nodes_parietal(iidx),RestingState,'RestingState',run,settings);

end

disp('Temporal');
idx_nodes_temporal = [81 82 83 84 85 86 87 88 89 90];

for iidx=1:length(idx_nodes_temporal)

    applyPCA(idx_nodes_temporal(iidx),MOT4,'High Attention',run,settings);
    applyPCA(idx_nodes_temporal(iidx),MOT2,'Low Attention',run,settings);
    applyPCA(idx_nodes_temporal(iidx),RestingState,'RestingState',run,settings);

end

return

function applyPCA(idx_nodes,Condition,Condition_Label,run,settings)

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

%% PCA
disp('PCA');

k = 12;

mean_area = mean(area,2);

X0_area = area - repmat(mean_area,1,nTR);

[U_area, S_area, V_area] = svd( X0_area' ); 

V_area = V_area';

Y0_area = V_area*X0_area;

diag_S_area = diag(S_area);

PC_voxels = zeros(1,nVoxels);

for iVoxel=1:nVoxels
    
    time_series = area(iVoxel,:);
    
    for iPC=1:k
        
        PC = Y0_area(iPC,:);
        
        difference(iPC,:) = time_series - PC;
        
    end
    
    diffsqr = difference.^2;
    sse = sqrt(sum(diffsqr,2));
    
    [x, I] = sort(sse);
    
    min_sse_PC = I(1);
    
    PC_voxels(iVoxel) = min_sse_PC;
    
end

maxK = k;

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
mi = 1;
ma = maxK;
colors = 1:maxK;
RGB = squeeze(ind2rgb(floor(((colors(:)-mi)/(ma-mi))*size(col,1)),col));

for ik=1:maxK
    
    subplot(ceil(maxK/2),2,ik);

    plot(Y0_area(ik,:),'Color',RGB(ik,:));
    
    hold on
    
    variance = 100 * ( diag_S_area(ik) / sum(diag_S_area));
    
    xlabel('TRs','FontSize',fs);
    ylabel('BOLD Activity','FontSize',fs);
    set(gca,'FontSize',fs);
    
    title(strcat('PC:',int2str(ik),'-','var:',num2str(variance),'%'));

end

 print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Condition_Label,'-',int2str(run),'-',arealabel,'-','PCA-PCs-',num2str(nVoxels),'.jpg'));
 print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Condition_Label,'-',int2str(run),'-',arealabel,'-','PCA-PCs-',num2str(nVoxels),'.eps'));
 print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Condition_Label,'-',int2str(run),'-',arealabel,'-','PCA-PCs-',num2str(nVoxels),'.pdf'));
 
 close all

 volume = zeros(size(AAL_img));
 
 iVoxel = 0;
 for iNode=1:nNodes
    
    idx_ROI = AAL_ROI(idx_nodes(iNode)).ID;
    
    idx_voxels = find(AAL_img == idx_ROI);
    
    for iiVoxel=1:length(idx_voxels)
        
        iVoxel = iVoxel + 1;
    
        [idxx, idxy, idxz] = ind2sub(size(AAL_img),idx_voxels(iiVoxel));
        
        ik = PC_voxels(iVoxel);
        
        volume(idxx,idxy,idxz) = colors(ik);
        
    end
    
 end

file_label = strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Condition_Label,'-',int2str(run),'-',arealabel,'-','PCA-',num2str(nVoxels));

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

f = figure;

[x, I] = sort(mean_area);

idx_max_six = I(end-5:end);
idx_min_six = I(1:6);

j = 0;
for i=1:length(idx_max_six)
    
    j = j + 1;
    
    subplot(6,2,j);
    
    plot(area(idx_max_six(i),:),'b');
    
    xlabel('TR');
    ylabel('BOLD Activity');
    title(strcat('HM Voxel:',int2str(idx_max_six(i))));
 
    hold on
    
end

for i=1:length(idx_min_six)
    
    j = j + 1;
    
    subplot(6,2,j);
    
    plot(area(idx_min_six(i),:),'r');
    
    xlabel('TR');
    ylabel('BOLD Activity');
    title(strcat('HM Voxel:',int2str(idx_min_six(i))));
 
    hold on
    
end

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Condition_Label,'-',int2str(run),'-',arealabel,'-','PCA-Voxels-',num2str(nVoxels),'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Condition_Label,'-',int2str(run),'-',arealabel,'-','PCA-Voxels-',num2str(nVoxels),'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Condition_Label,'-',int2str(run),'-',arealabel,'-','PCA-Voxels-',num2str(nVoxels),'.pdf'));
 
close all

f = figure;
 
S = S_area;

bar( log( diag(S) ) );              % plot variance of principal components
hold off;
title( arealabel );

x = length(log(diag(S)))/4;
y = min(log(diag(S))) + (max(log(diag(S))) - min(log(diag(S))))/2;
text(x,y,strcat('#voxels:',num2str(nVoxels)),'FontSize',fs);

n95th = get95thvariance(diag_S_area);

y = min(log(diag(S))) + (max(log(diag(S))) - min(log(diag(S))))/2 + 2;
x = length(log(diag(S)))/4;
text(x,y,strcat('#PCs (95%):',num2str(n95th)),'FontSize',fs);

xlabel( 'PCs');
ylabel( 'log variance' );

xlim([0 length(log(diag(S)))]);
ylim([min(log(diag(S))) max(log(diag(S)))]);

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Condition_Label,'-',int2str(run),'-',arealabel,'-','PCA-Variance-',num2str(nVoxels),'.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Condition_Label,'-',int2str(run),'-',arealabel,'-','PCA-Variance-',num2str(nVoxels),'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-',Condition_Label,'-',int2str(run),'-',arealabel,'-','PCA-Variance-',num2str(nVoxels),'.pdf'));
 
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