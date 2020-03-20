
function lowhigh_spotlight_Precuneus_Angular_FC_lobe


settings_jan_0805;
%settings_elena_2905;

%calculateCentroids;

seeCentroids;

%FC_mean_Spotlight_Lobe_IC(settings);

%plotFC(settings);


end

function calculateCentroids

    nCentroids = 40;

    load_aal = nifti('angular_l.nii');
    load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_area_masks\',load_aal.dat.fname);
    Angular_L_img = load_aal.dat(:,:,:);

    load_aal = nifti('angular_r.nii');
    load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_area_masks\',load_aal.dat.fname);
    Angular_R_img = load_aal.dat(:,:,:);

    load_aal = nifti('precuneus_l.nii');
    load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_area_masks\',load_aal.dat.fname);
    Precuneus_L_img = load_aal.dat(:,:,:);

    load_aal = nifti('precuneus_r.nii');
    load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_area_masks\',load_aal.dat.fname);
    Precuneus_R_img = load_aal.dat(:,:,:);

    disp('Angular_L');
    idx_voxels = find(Angular_L_img);
    
    for iVoxel=1:length(idx_voxels)
        
        nInside = round( length(idx_voxels) / nCentroids );
        parcel = ceil(iVoxel/nInside);
        
        parcel
        
        [idxx,idxy,idxz] = ind2sub(size(Angular_L_img),idx_voxels(iVoxel));
    
        Angular_L_img(idxx,idxy,idxz) = parcel;
        
    end
        
    disp('Angular_R');
    idx_voxels = find(Angular_R_img);
    
    for iVoxel=1:length(idx_voxels)
        
        nInside = round( length(idx_voxels) / nCentroids );
        parcel = ceil(iVoxel/nInside);
        
        parcel
        
        [idxx,idxy,idxz] = ind2sub(size(Angular_R_img),idx_voxels(iVoxel));
    
        Angular_R_img(idxx,idxy,idxz) = parcel;
        
    end
    
    disp('Precuneus_L'),
    idx_voxels = find(Precuneus_L_img);
    
    for iVoxel=1:length(idx_voxels)
        
        nInside = round( length(idx_voxels) / nCentroids );
        parcel = ceil(iVoxel/nInside);
        
        parcel
        
        [idxx,idxy,idxz] = ind2sub(size(Precuneus_L_img),idx_voxels(iVoxel));
    
        Precuneus_L_img(idxx,idxy,idxz) = parcel;
        
    end
    
    disp('Precuneus_R');
    idx_voxels = find(Precuneus_R_img);

    for iVoxel=1:length(idx_voxels)
        
        nInside = round( length(idx_voxels) / nCentroids );
        parcel = ceil(iVoxel/nInside);
        
        parcel
        
        [idxx,idxy,idxz] = ind2sub(size(Precuneus_R_img),idx_voxels(iVoxel));
    
        Precuneus_R_img(idxx,idxy,idxz) = parcel;
        
    end
    
    dim = size(load_aal.dat);
    nifti_file = load_aal;
    offset = load_aal.dat.offset;
    scl_slope = load_aal.dat.scl_slope;
    scl_inter = load_aal.dat.scl_inter;
    dtype = 'FLOAT32';

    descrip = 'centroids';

    fname = strcat('angular_l.nii');
    input_data = Angular_L_img; 
    lowhigh_save_image;     
    
    fname = strcat('angular_r.nii');
    input_data = Angular_R_img; 
    lowhigh_save_image;  
    
    fname = strcat('precuneus_l.nii');
    input_data = Precuneus_L_img; 
    lowhigh_save_image;   
    
    fname = strcat('precuneus_r.nii');
    input_data = Precuneus_R_img; 
    lowhigh_save_image;     

end

function seeCentroids

    rendfile = 'K:\Dropbox (Uni Magdeburg)\_TOOLBOX\spm8\canonical\cortex_20484.surf.gii';
    brt = 1;

    filename{1} = 'angular_l.nii';
    filename{2} = 'angular_r.nii';
    filename{3} = 'precuneus_l.nii';
    filename{4} = 'precuneus_r.nii';

    prefix = 'K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_area_centroids\40\';

    for ifile=1:4

        filename{ifile} = strcat(prefix,filename{ifile});

    end

    for ifile=1:4

        img(ifile).dat = fromImageToSPMdat(filename(ifile));

    end

    min_value = 0;
    max_value = 40;

    setGlobalMinMax(min_value,max_value);
    snapshot_label = '';
    setRenderingSnapshot(false,snapshot_label);

    for ifile=1:4

        spm_render_my(img(ifile).dat,brt,rendfile,filename{ifile});

        pause;

    end

end

function FC_mean_Spotlight_Lobe_IC(settings)

    nCentroids = 40;
    
    load_aal = nifti('angular_l.nii');
    load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_area_masks\',load_aal.dat.fname);
    Angular_L_img = load_aal.dat(:,:,:);

    load_aal = nifti('angular_r.nii');
    load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_area_masks\',load_aal.dat.fname);
    Angular_R_img = load_aal.dat(:,:,:);

    load_aal = nifti('precuneus_l.nii');
    load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_area_masks\',load_aal.dat.fname);
    Precuneus_L_img = load_aal.dat(:,:,:);

    load_aal = nifti('precuneus_r.nii');
    load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_area_masks\',load_aal.dat.fname);
    Precuneus_R_img = load_aal.dat(:,:,:);

    disp('Angular_L');
    idx_voxels = find(Angular_L_img);
    
    for iVoxel=1:length(idx_voxels)
        
        nInside = round( length(idx_voxels) / nCentroids );
        parcel = ceil(iVoxel/nInside);
        
        [idxx,idxy,idxz] = ind2sub(size(Angular_L_img),idx_voxels(iVoxel));
    
        Angular_L_img(idxx,idxy,idxz) = parcel;
        
    end
        
    disp('Angular_R');
    idx_voxels = find(Angular_R_img);
    
    for iVoxel=1:length(idx_voxels)
        
        nInside = round( length(idx_voxels) / nCentroids );
        parcel = ceil(iVoxel/nInside);
        
        [idxx,idxy,idxz] = ind2sub(size(Angular_R_img),idx_voxels(iVoxel));
    
        Angular_R_img(idxx,idxy,idxz) = parcel;
        
    end
    
    disp('Precuneus_L'),
    idx_voxels = find(Precuneus_L_img);
    
    for iVoxel=1:length(idx_voxels)
        
        nInside = round( length(idx_voxels) / nCentroids );
        parcel = ceil(iVoxel/nInside);
        
        [idxx,idxy,idxz] = ind2sub(size(Precuneus_L_img),idx_voxels(iVoxel));
    
        Precuneus_L_img(idxx,idxy,idxz) = parcel;
        
    end
    
    disp('Precuneus_R');
    idx_voxels = find(Precuneus_R_img);

    for iVoxel=1:length(idx_voxels)
        
        nInside = round( length(idx_voxels) / nCentroids );
        parcel = ceil(iVoxel/nInside);
        
        [idxx,idxy,idxz] = ind2sub(size(Precuneus_R_img),idx_voxels(iVoxel));
    
        Precuneus_R_img(idxx,idxy,idxz) = parcel;
        
    end
    
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;
get_at_this_preprocessed_step = settings.FSL.folders.custom;

lowhigh_load_all_data_FSL;
    
ICA_preprocessed = 'groupICA';

nTR = 300;

area_lobe_preprocessed_folder = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',ICA_preprocessed,'\','area_lobe');

run = 1;
occipital_ic_run1 = load(strcat(area_lobe_preprocessed_folder,'\','high-rest-run-',int2str(run),'-','occipital_lobe','\','melodic_mix'));
high_occipital_ic_run1 = occipital_ic_run1(1:nTR,1:nCentroids);
rest_occipital_ic_run1 = occipital_ic_run1(nTR+1:end,1:nCentroids);

run = 2;
occipital_ic_run2 = load(strcat(area_lobe_preprocessed_folder,'\','high-rest-run-',int2str(run),'-','occipital_lobe','\','melodic_mix'));
high_occipital_ic_run2 = occipital_ic_run2(1:nTR,1:nCentroids);
rest_occipital_ic_run2 = occipital_ic_run2(nTR+1:end,1:nCentroids);

[centroid_MOT4Run1_angular_l, centroid_MOT4Run2_angular_l, centroid_RestingStateRun1_angular_l, centroid_RestingStateRun2_angular_l] = getCentroids(nTR,nCentroids,MOT4Run1,MOT4Run2,RestingStateRun1,RestingStateRun2,Angular_L_img);

[centroid_MOT4Run1_angular_r, centroid_MOT4Run2_angular_r, centroid_RestingStateRun1_angular_r, centroid_RestingStateRun2_angular_r] = getCentroids(nTR,nCentroids,MOT4Run1,MOT4Run2,RestingStateRun1,RestingStateRun2,Angular_R_img);

[centroid_MOT4Run1_precuneus_l, centroid_MOT4Run2_precuneus_l, centroid_RestingStateRun1_precuneus_l, centroid_RestingStateRun2_precuneus_l] = getCentroids(nTR,nCentroids,MOT4Run1,MOT4Run2,RestingStateRun1,RestingStateRun2,Precuneus_L_img);

[centroid_MOT4Run1_precuneus_r, centroid_MOT4Run2_precuneus_r, centroid_RestingStateRun1_precuneus_r, centroid_RestingStateRun2_precuneus_r] = getCentroids(nTR,nCentroids,MOT4Run1,MOT4Run2,RestingStateRun1,RestingStateRun2,Precuneus_R_img);

%%% ANGULAR L

[FC_spotlight(1).rho_high_angular_l_occipital, FC_spotlight(1).pval_high_angular_l_occipital] = corr([high_occipital_ic_run1,centroid_MOT4Run1_angular_l]);

[FC_spotlight(1).rho_rest_angular_l_occipital, FC_spotlight(1).pval_rest_angular_l_occipital] = corr([rest_occipital_ic_run1,centroid_RestingStateRun1_angular_l]);

[FC_spotlight(2).rho_high_angular_l_occipital, FC_spotlight(2).pval_high_angular_l_occipital] = corr([high_occipital_ic_run2,centroid_MOT4Run2_angular_l]);

[FC_spotlight(2).rho_rest_angular_l_occipital, FC_spotlight(2).pval_rest_angular_l_occipital] = corr([rest_occipital_ic_run2,centroid_RestingStateRun2_angular_l]);

%%% ANGULAR R

[FC_spotlight(1).rho_high_angular_r_occipital, FC_spotlight(1).pval_high_angular_r_occipital] = corr([high_occipital_ic_run1,centroid_MOT4Run1_angular_r]);

[FC_spotlight(1).rho_rest_angular_r_occipital, FC_spotlight(1).pval_rest_angular_r_occipital] = corr([rest_occipital_ic_run1,centroid_RestingStateRun1_angular_r]);

[FC_spotlight(2).rho_high_angular_r_occipital, FC_spotlight(2).pval_high_angular_r_occipital] = corr([high_occipital_ic_run2,centroid_MOT4Run2_angular_r]);

[FC_spotlight(2).rho_rest_angular_r_occipital, FC_spotlight(2).pval_rest_angular_r_occipital] = corr([rest_occipital_ic_run2,centroid_RestingStateRun2_angular_r]);

%%% PRECUNEUS L

[FC_spotlight(1).rho_high_precuneus_l_occipital, FC_spotlight(1).pval_high_precuneus_l_occipital] = corr([high_occipital_ic_run1,centroid_MOT4Run1_precuneus_l]);

[FC_spotlight(1).rho_rest_precuneus_l_occipital, FC_spotlight(1).pval_rest_precuneus_l_occipital] = corr([rest_occipital_ic_run1,centroid_RestingStateRun1_precuneus_l]);

[FC_spotlight(2).rho_high_precuneus_l_occipital, FC_spotlight(2).pval_high_precuneus_l_occipital] = corr([high_occipital_ic_run2,centroid_MOT4Run2_precuneus_l]);

[FC_spotlight(2).rho_rest_precuneus_l_occipital, FC_spotlight(2).pval_rest_precuneus_l_occipital] = corr([rest_occipital_ic_run2,centroid_RestingStateRun2_precuneus_l]);

%%% PRECUNEUS R

[FC_spotlight(1).rho_high_precuneus_r_occipital, FC_spotlight(1).pval_high_precuneus_r_occipital] = corr([high_occipital_ic_run1,centroid_MOT4Run1_precuneus_r]);

[FC_spotlight(1).rho_rest_precuneus_r_occipital, FC_spotlight(1).pval_rest_precuneus_r_occipital] = corr([rest_occipital_ic_run1,centroid_RestingStateRun1_precuneus_r]);

[FC_spotlight(2).rho_high_precuneus_r_occipital, FC_spotlight(2).pval_high_precuneus_r_occipital] = corr([high_occipital_ic_run2,centroid_MOT4Run2_precuneus_r]);

[FC_spotlight(2).rho_rest_precuneus_r_occipital, FC_spotlight(2).pval_rest_precuneus_r_occipital] = corr([rest_occipital_ic_run2,centroid_RestingStateRun2_precuneus_r]);

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Lobe-Spotlight-Precuneus-Angular-All-Runs','.mat'),'FC_spotlight');

end

function [mean_centroid_MOT4Run1, mean_centroid_MOT4Run2, mean_centroid_RestingStateRun1, mean_centroid_RestingStateRun2] = getCentroids(nTR,nCentroids,MOT4Run1,MOT4Run2,RestingStateRun1,RestingStateRun2,area_img)

mean_centroid_MOT4Run1 = zeros(nTR,nCentroids);
mean_centroid_MOT4Run2 = zeros(nTR,nCentroids);
mean_centroid_RestingStateRun1 = zeros(nTR,nCentroids);
mean_centroid_RestingStateRun2 = zeros(nTR,nCentroids);

for iCentroid=1:nCentroids
    
    idx_voxels = find(area_img==iCentroid);
    
    centroid_MOT4Run1 = zeros(nTR,length(idx_voxels));
    centroid_MOT4Run2 = zeros(nTR,length(idx_voxels));
    
    centroid_RestingStateRun1 = zeros(nTR,length(idx_voxels));
    centroid_RestingStateRun2 = zeros(nTR,length(idx_voxels));
    
    for iVoxel=1:length(idx_voxels)
       
        [idxx,idxy,idxz] = ind2sub(size(area_img),idx_voxels(iVoxel));
        
        centroid_MOT4Run1(1:nTR,iVoxel) = MOT4Run1(idxx,idxy,idxz,:);
        centroid_MOT4Run2(1:nTR,iVoxel) = MOT4Run2(idxx,idxy,idxz,:);
        
        centroid_RestingStateRun1(1:nTR,iVoxel) = RestingStateRun1(idxx,idxy,idxz,:);
        centroid_RestingStateRun2(1:nTR,iVoxel) = RestingStateRun2(idxx,idxy,idxz,:);
        
    end
    
     mean_centroid_MOT4Run1(:,iCentroid) = mean(centroid_MOT4Run1,2);
     mean_centroid_MOT4Run2(:,iCentroid) = mean(centroid_MOT4Run2,2);
        
     mean_centroid_RestingStateRun1(:,iCentroid) = mean(centroid_RestingStateRun1,2);
     mean_centroid_RestingStateRun2(:,iCentroid) = mean(centroid_RestingStateRun2,2);

end

end

function plotFC(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Lobe-Spotlight-Precuneus-Angular-All-Runs','.mat'),'FC_spotlight');

Ncomponent = 40;

ROIlabel = 'Angular-L-Run-1';
plotRho(settings,ROIlabel,FC_spotlight(1).rho_high_angular_l_occipital,FC_spotlight(1).pval_high_angular_l_occipital,FC_spotlight(1).rho_rest_angular_l_occipital,FC_spotlight(1).pval_rest_angular_l_occipital,Ncomponent);
ROIlabel = 'Angular-L-Run-2';
plotRho(settings,ROIlabel,FC_spotlight(2).rho_high_angular_l_occipital,FC_spotlight(2).pval_high_angular_l_occipital,FC_spotlight(2).rho_rest_angular_l_occipital,FC_spotlight(2).pval_rest_angular_l_occipital,Ncomponent);

ROIlabel = 'Angular-R-Run-1';
plotRho(settings,ROIlabel,FC_spotlight(1).rho_high_angular_r_occipital,FC_spotlight(1).pval_high_angular_r_occipital,FC_spotlight(1).rho_rest_angular_r_occipital,FC_spotlight(1).pval_rest_angular_r_occipital,Ncomponent);
ROIlabel = 'Angular-R-Run-2';
plotRho(settings,ROIlabel,FC_spotlight(2).rho_high_angular_r_occipital,FC_spotlight(2).pval_high_angular_r_occipital,FC_spotlight(2).rho_rest_angular_r_occipital,FC_spotlight(2).pval_rest_angular_r_occipital,Ncomponent);

ROIlabel = 'Precuneus-L-Run-1';
plotRho(settings,ROIlabel,FC_spotlight(1).rho_high_precuneus_l_occipital,FC_spotlight(1).pval_high_precuneus_l_occipital,FC_spotlight(1).rho_rest_precuneus_l_occipital,FC_spotlight(1).pval_rest_precuneus_l_occipital,Ncomponent);
ROIlabel = 'Precuneus-L-Run-2';
plotRho(settings,ROIlabel,FC_spotlight(2).rho_high_precuneus_l_occipital,FC_spotlight(2).pval_high_precuneus_l_occipital,FC_spotlight(2).rho_rest_precuneus_l_occipital,FC_spotlight(2).pval_rest_precuneus_l_occipital,Ncomponent);

ROIlabel = 'Precuneus-R-Run-1';
plotRho(settings,ROIlabel,FC_spotlight(1).rho_high_precuneus_r_occipital,FC_spotlight(1).pval_high_precuneus_r_occipital,FC_spotlight(1).rho_rest_precuneus_r_occipital,FC_spotlight(1).pval_rest_precuneus_r_occipital,Ncomponent);
ROIlabel = 'Precuneus-R-Run-2';
plotRho(settings,ROIlabel,FC_spotlight(2).rho_high_precuneus_r_occipital,FC_spotlight(2).pval_high_precuneus_r_occipital,FC_spotlight(2).rho_rest_precuneus_r_occipital,FC_spotlight(2).pval_rest_precuneus_r_occipital,Ncomponent);

end

function plotRho(settings,ROIlabel,rho_att,pval_att,rho_rest,pval_rest,Ncomponent)
      
fs = 14;

pcriterion = 0.001;
rholimit   = 0.5;

rhomin = -rholimit + rholimit/32;
rhomax =  rholimit;

rhorange = linspace( rhomin, rhomax, 64 );

kk = find( pval_att > pcriterion );
rho_att(kk) = zeros(size(kk));

kk = find( pval_rest > pcriterion );
rho_rest(kk) = zeros(size(kk));


clrmp = colormap('jet');
clrmp(32,:) = [1 1 1];

f = figure;

subplot(1,2,1);
hold on;
caxis([rhomin rhomax]);
h = pcolor( rho_att );
set( h, 'EdgeColor', 'none');
colormap(clrmp);
plot( 1+ Ncomponent*[1 1], Ncomponent*[0 2], 'k', 'LineWidth', 2 );
plot( Ncomponent*[0 2], 1+ Ncomponent*[1 1], 'k', 'LineWidth', 2 );
hold off;
axis 'equal'; 
Tick = [Ncomponent/2 3*Ncomponent/2];
TickLabel = { 'Lobe' ; 'ROI' };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );
title('attention', 'FontSize', fs );


h = colorbar;
set(h,'YLim',[-0.4 0.4], 'YTick',[-0.4 -0.2 0 0.2 0.4], 'PlotBoxAspectRatio', [1 20 1]);


subplot(1,2,2);
hold on;
caxis([rhomin rhomax]);
h = pcolor( rho_rest );
set( h, 'EdgeColor', 'none');
colormap(clrmp);
plot( 1+ Ncomponent*[1 1], Ncomponent*[0 2], 'k', 'LineWidth', 2 );
plot( Ncomponent*[0 2], 1+ Ncomponent*[1 1], 'k', 'LineWidth', 2 );
hold off;
axis 'equal'; 
Tick = [Ncomponent/2 3*Ncomponent/2];
TickLabel = { 'Lobe' ; 'ROI' };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );
title('resting', 'FontSize', fs );


h = colorbar;
set(h,'YLim',[-0.4 0.4], 'YTick',[-0.4 -0.2 0 0.2 0.4], 'PlotBoxAspectRatio', [1 20 1]);

suptitle( ROIlabel );

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Lobe-Spotlight','-',ROIlabel,'.jpeg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Lobe-Spotlight','-',ROIlabel,'.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-Lobe-Spotlight','-',ROIlabel,'.pdf'));


end