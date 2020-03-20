function lowhigh_granger_functional

settings_jan_0805;
%settings_elena_2905;

%computeGrangerForFunctionalNetworksOnlyOnSeeds(settings);
%computeGrangerForFunctionalNetworksOnlyOnSeedsMatrix(settings);

% computeGrangerForFunctionalNetworksOnlyOnSeedsMean(settings);

%plotGrangerForFunctionalNetworksOnlyOnSeeds(settings);
plotGrangerForFunctionalNetworksOnlyOnSeedsOnlyNetworks(settings);

end


function computeGrangerForFunctionalNetworksOnlyOnSeeds(settings)

analysis_label = 'Granger-Func-Seeds';

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

lowhigh_load_all_data_FSL;

all_networks_labels = {'DAN','VAN','VIS','LAN','DMN','FPC','SMN','AUD'};
%all_networks_labels = {'DAN','VAN','VIS','LAN'};
%all_networks_labels = {'DAN','VIS'};

MNI_dim = [91 109 91];

MNI_x_center = 45;
MNI_y_center = 63;
MNI_z_center = 36;

size_voxels_mm = 2;
ROI_radius_mm = 6 - size_voxels_mm;
ROI_radius_voxels = ROI_radius_mm/size_voxels_mm;

iiROI = 0;
for iNet=1:length(all_networks_labels)
    
    network_label = all_networks_labels{iNet};
    seeds = getFunctionalSeeds(network_label);
    
    nROI = length(seeds.ROI);
    
    for iROI=1:nROI
        
        iiROI = iiROI + 1;
        
        ROI(iiROI).netLabel = seeds.network_label;
        ROI(iiROI).seedLabel = seeds.ROI(iROI).label;
        
        ROI(iiROI).x = seeds.ROI(iROI).x;
        ROI(iiROI).y = seeds.ROI(iROI).y;
        ROI(iiROI).z = seeds.ROI(iROI).z;

        x_center = MNI_x_center + round(ROI(iiROI).x/size_voxels_mm);  
        y_center = MNI_y_center + round(ROI(iiROI).y/size_voxels_mm);
        z_center = MNI_z_center + round(ROI(iiROI).z/size_voxels_mm);

        xgv = (x_center-ROI_radius_voxels):(x_center+ROI_radius_voxels);
        ygv = (y_center-ROI_radius_voxels):(y_center+ROI_radius_voxels);
        zgv = (z_center-ROI_radius_voxels):(z_center+ROI_radius_voxels);

        [ROI(iiROI).X,ROI(iiROI).Y,ROI(iiROI).Z] = meshgrid(xgv,ygv,zgv);
        
        nVoxels = length(find(ROI(iiROI).X));
        
        ROI(iiROI).idx_voxels = zeros(1,nVoxels);
        
        %nVoxels = 25;
        for iVoxel=1:nVoxels
            
            ROI(iiROI).idx_voxels(iVoxel) = sub2ind(MNI_dim,ROI(iiROI).X(iVoxel),ROI(iiROI).Y(iVoxel),ROI(iiROI).Z(iVoxel));
            
        end
    
    end
    
end

all_idx_voxels = zeros(1,length(ROI)*nVoxels);
iiVoxel = 0;
for iiROI=1:length(ROI)
    
    for iVoxel=1:nVoxels
        
        iiVoxel = iiVoxel + 1;
        
        all_idx_voxels(iiVoxel) = ROI(iiROI).idx_voxels(iVoxel);
        
    end
    
end

MOT4Run1_GF = zeros(length(all_idx_voxels),length(all_idx_voxels));
MOT4Run1_Gpval = zeros(length(all_idx_voxels),length(all_idx_voxels));
MOT4Run2_GF = zeros(length(all_idx_voxels),length(all_idx_voxels));
MOT4Run2_Gpval = zeros(length(all_idx_voxels),length(all_idx_voxels));

RestingStateRun1_GF = zeros(length(all_idx_voxels),length(all_idx_voxels));
RestingStateRun1_Gpval = zeros(length(all_idx_voxels),length(all_idx_voxels));
RestingStateRun2_GF = zeros(length(all_idx_voxels),length(all_idx_voxels));
RestingStateRun2_Gpval = zeros(length(all_idx_voxels),length(all_idx_voxels));

start = now;

for iiVoxel=1:length(all_idx_voxels)
    
    tic
    disp(strcat('doing iiVoxel:',int2str(iiVoxel)));
    
    [idxx,idxy,idxz] = ind2sub(MNI_dim,all_idx_voxels(iiVoxel));
    
    MOT4Run1_voxel1 = squeeze(MOT4Run1(idxx,idxy,idxz,:));
    MOT4Run2_voxel1 = squeeze(MOT4Run2(idxx,idxy,idxz,:));
    
    RestingStateRun1_voxel1 = squeeze(RestingStateRun1(idxx,idxy,idxz,:));
    RestingStateRun2_voxel1 = squeeze(RestingStateRun2(idxx,idxy,idxz,:));
    
    for iiiVoxel=(iiVoxel+1):length(all_idx_voxels)
        
        [iidxx,iidxy,iidxz] = ind2sub(MNI_dim,all_idx_voxels(iiiVoxel));
        
        MOT4Run1_voxel2 = squeeze(MOT4Run1(iidxx,iidxy,iidxz,:));
        MOT4Run2_voxel2 = squeeze(MOT4Run2(iidxx,iidxy,iidxz,:));
    
        RestingStateRun1_voxel2 = squeeze(RestingStateRun1(iidxx,iidxy,iidxz,:));
        RestingStateRun2_voxel2 = squeeze(RestingStateRun2(iidxx,iidxy,iidxz,:));
    
        [GF,Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(MOT4Run1_voxel1,MOT4Run1_voxel2); 
        
        MOT4Run1_GF(iiiVoxel,iiVoxel) = GF(2,1);
        MOT4Run1_GF(iiVoxel,iiiVoxel) = GF(1,2);
        MOT4Run1_Gpval(iiiVoxel,iiVoxel) = Gpval(2,1);
        MOT4Run1_Gpval(iiVoxel,iiiVoxel) = Gpval(1,2);
 
%         [GF,Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(MOT4Run2_voxel1,MOT4Run2_voxel2); 
%         
%         MOT4Run2_GF(iiiVoxel,iiVoxel) = GF(2,1);
%         MOT4Run2_GF(iiVoxel,iiiVoxel) = GF(1,2);
%         MOT4Run2_Gpval(iiiVoxel,iiVoxel) = Gpval(2,1);
%         MOT4Run2_Gpval(iiVoxel,iiiVoxel) = Gpval(1,2);
 
%         [GF,Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(RestingStateRun1_voxel1,RestingStateRun1_voxel2); 
%         
%         RestingStateRun1_GF(iiiVoxel,iiVoxel) = GF(2,1);
%         RestingStateRun1_GF(iiVoxel,iiiVoxel) = GF(1,2);
%         RestingStateRun1_Gpval(iiiVoxel,iiVoxel) = Gpval(2,1);
%         RestingStateRun1_Gpval(iiVoxel,iiiVoxel) = Gpval(1,2);
 
%         [GF,Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(RestingStateRun2_voxel1,RestingStateRun2_voxel2); 
%         
%         RestingStateRun2_GF(iiiVoxel,iiVoxel) = GF(2,1);
%         RestingStateRun2_GF(iiVoxel,iiiVoxel) = GF(1,2);
%         RestingStateRun2_Gpval(iiiVoxel,iiVoxel) = Gpval(2,1);
%         RestingStateRun2_Gpval(iiVoxel,iiiVoxel) = Gpval(1,2);
 
    end
    toc
    
end

finished = now;

lasted = finished - start;

save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'.mat'),'all_networks_labels','ROI','all_idx_voxels','MOT4Run1_GF','MOT4Run2_GF','RestingStateRun1_GF','RestingStateRun2_GF','MOT4Run1_Gpval','MOT4Run2_Gpval','RestingStateRun1_Gpval','RestingStateRun2_Gpval','start','finished','lasted');

end

function computeGrangerForFunctionalNetworksOnlyOnSeedsMean(settings)

analysis_label = 'Granger-Func-Seeds-Mean';

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

lowhigh_load_all_data_FSL;

%all_networks_labels = {'DAN','VAN','VIS','LAN','DMN','FPC','SMN','AUD'};
all_networks_labels = {'DAN','VAN','VIS','LAN'};

MNI_dim = [91 109 91];

MNI_x_center = 45;
MNI_y_center = 63;
MNI_z_center = 36;

size_voxels_mm = 2;
ROI_radius_mm = 6 - size_voxels_mm;
ROI_radius_voxels = ROI_radius_mm/size_voxels_mm;

iiROI = 0;
for iNet=1:length(all_networks_labels)
    
    network_label = all_networks_labels{iNet};
    seeds = getFunctionalSeeds(network_label);
    
    nROI = length(seeds.ROI);
    
    for iROI=1:nROI
        
        iiROI = iiROI + 1;
        
        ROI(iiROI).netLabel = seeds.network_label;
        ROI(iiROI).seedLabel = seeds.ROI(iROI).label;
        
        ROI(iiROI).x = seeds.ROI(iROI).x;
        ROI(iiROI).y = seeds.ROI(iROI).y;
        ROI(iiROI).z = seeds.ROI(iROI).z;

        x_center = MNI_x_center + round(ROI(iiROI).x/size_voxels_mm);  
        y_center = MNI_y_center + round(ROI(iiROI).y/size_voxels_mm);
        z_center = MNI_z_center + round(ROI(iiROI).z/size_voxels_mm);

        xgv = (x_center-ROI_radius_voxels):(x_center+ROI_radius_voxels);
        ygv = (y_center-ROI_radius_voxels):(y_center+ROI_radius_voxels);
        zgv = (z_center-ROI_radius_voxels):(z_center+ROI_radius_voxels);

        [ROI(iiROI).X,ROI(iiROI).Y,ROI(iiROI).Z] = meshgrid(xgv,ygv,zgv);
        
        nVoxels = length(find(ROI(iiROI).X));
        
        ROI(iiROI).idx_voxels = zeros(1,nVoxels);
        
        nVoxels = 25;
        for iVoxel=1:nVoxels
            
            ROI(iiROI).idx_voxels(iVoxel) = sub2ind(MNI_dim,ROI(iiROI).X(iVoxel),ROI(iiROI).Y(iVoxel),ROI(iiROI).Z(iVoxel));
            
        end
    
    end
    
end

nTR = size(MOT4Run1,4);
MOT4Run1_all_ROI_mean = zeros(length(ROI),nTR);
MOT4Run2_all_ROI_mean = zeros(length(ROI),nTR);
RestingStateRun1_all_ROI_mean = zeros(length(ROI),nTR);
RestingStateRun2_all_ROI_mean = zeros(length(ROI),nTR);
for iiROI=1:length(ROI)
    
    MOT4Run1_voxels_ts = zeros(nVoxels,nTR);
    MOT4Run2_voxels_ts = zeros(nVoxels,nTR);
    RestingStateRun1_voxels_ts = zeros(nVoxels,nTR);
    RestingStateRun2_voxels_ts = zeros(nVoxels,nTR);
    
    for iVoxel=1:nVoxels
        
        [idxx,idxy,idxz] = ind2sub(MNI_dim,ROI(iiROI).idx_voxels(iVoxel));
        
        MOT4Run1_voxels_ts(iVoxel,:) = squeeze(MOT4Run1(idxx,idxy,idxz,:));
        MOT4Run2_voxels_ts(iVoxel,:) = squeeze(MOT4Run2(idxx,idxy,idxz,:));

        RestingStateRun1_voxels_ts(iVoxel,:) = squeeze(RestingStateRun1(idxx,idxy,idxz,:));
        RestingStateRun2_voxels_ts(iVoxel,:) = squeeze(RestingStateRun2(idxx,idxy,idxz,:));

    end
    
    MOT4Run1_all_ROI_mean(iiROI,:) = mean(MOT4Run1_voxels_ts,1);
    MOT4Run2_all_ROI_mean(iiROI,:) = mean(MOT4Run2_voxels_ts,1);
    
    RestingStateRun1_all_ROI_mean(iiROI,:) = mean(RestingStateRun1_voxels_ts,1);
    RestingStateRun2_all_ROI_mean(iiROI,:) = mean(RestingStateRun2_voxels_ts,1);
        
end

nROI = length(ROI);

MOT4Run1_GF = zeros(nROI,nROI);
MOT4Run1_Gpval = zeros(nROI,nROI);
MOT4Run2_GF = zeros(nROI,nROI);
MOT4Run2_Gpval = zeros(nROI,nROI);

RestingStateRun1_GF = zeros(nROI,nROI);
RestingStateRun1_Gpval = zeros(nROI,nROI);
RestingStateRun2_GF = zeros(nROI,nROI);
RestingStateRun2_Gpval = zeros(nROI,nROI);

for iiROI=1:nROI
    
    disp(strcat('doing iiROI:',int2str(iiROI)));
    
    for iiiROI=(iiROI+1):nROI
        
        [GF,Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(MOT4Run1_all_ROI_mean(iiROI,:),MOT4Run1_all_ROI_mean(iiiROI,:)); 
        
        MOT4Run1_GF(iiiROI,iiROI) = GF(2,1);
        MOT4Run1_GF(iiROI,iiiROI) = GF(1,2);
        MOT4Run1_Gpval(iiiROI,iiROI) = Gpval(2,1);
        MOT4Run1_Gpval(iiROI,iiiROI) = Gpval(1,2);
 
        [GF,Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(MOT4Run2_all_ROI_mean(iiROI,:),MOT4Run2_all_ROI_mean(iiiROI,:)); 
        
        MOT4Run2_GF(iiiROI,iiROI) = GF(2,1);
        MOT4Run2_GF(iiROI,iiiROI) = GF(1,2);
        MOT4Run2_Gpval(iiiROI,iiROI) = Gpval(2,1);
        MOT4Run2_Gpval(iiROI,iiiROI) = Gpval(1,2);
 
        [GF,Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(RestingStateRun1_all_ROI_mean(iiROI,:),RestingStateRun1_all_ROI_mean(iiiROI,:)); 
        
        RestingStateRun1_GF(iiiROI,iiROI) = GF(2,1);
        RestingStateRun1_GF(iiROI,iiiROI) = GF(1,2);
        RestingStateRun1_Gpval(iiiROI,iiROI) = Gpval(2,1);
        RestingStateRun1_Gpval(iiROI,iiiROI) = Gpval(1,2);
 
        [GF,Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(RestingStateRun2_all_ROI_mean(iiROI,:),RestingStateRun2_all_ROI_mean(iiiROI,:)); 
        
        RestingStateRun2_GF(iiiROI,iiROI) = GF(2,1);
        RestingStateRun2_GF(iiROI,iiiROI) = GF(1,2);
        RestingStateRun2_Gpval(iiiROI,iiROI) = Gpval(2,1);
        RestingStateRun2_Gpval(iiROI,iiiROI) = Gpval(1,2);
 
    end
    
end

save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'.mat'),'all_networks_labels','ROI','MOT4Run1_GF','MOT4Run2_GF','RestingStateRun1_GF','RestingStateRun2_GF','MOT4Run1_Gpval','MOT4Run2_Gpval','RestingStateRun1_Gpval','RestingStateRun2_Gpval');

end

function computeGrangerForFunctionalNetworksOnlyOnSeedsMatrix(settings)

analysis_label = 'Granger-Func-Seeds-Matrix';

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

lowhigh_load_all_data_FSL;

%all_networks_labels = {'DAN','VAN','VIS','LAN','DMN','FPC','SMN','AUD'};
all_networks_labels = {'DAN','VAN','VIS','LAN'};

MNI_dim = [91 109 91];
nTR = size(MOT4Run1,4);

MNI_x_center = 45;
MNI_y_center = 63;
MNI_z_center = 36;

size_voxels_mm = 2;
ROI_radius_mm = 6 - size_voxels_mm;
ROI_radius_voxels = ROI_radius_mm/size_voxels_mm;

iiROI = 0;
for iNet=1:length(all_networks_labels)
    
    network_label = all_networks_labels{iNet};
    seeds = getFunctionalSeeds(network_label);
    
    nROI = length(seeds.ROI);
    
    for iROI=1:nROI
        
        iiROI = iiROI + 1;
        
        ROI(iiROI).netLabel = seeds.network_label;
        ROI(iiROI).seedLabel = seeds.ROI(iROI).label;
        
        ROI(iiROI).x = seeds.ROI(iROI).x;
        ROI(iiROI).y = seeds.ROI(iROI).y;
        ROI(iiROI).z = seeds.ROI(iROI).z;

        x_center = MNI_x_center + round(ROI(iiROI).x/size_voxels_mm);  
        y_center = MNI_y_center + round(ROI(iiROI).y/size_voxels_mm);
        z_center = MNI_z_center + round(ROI(iiROI).z/size_voxels_mm);

        xgv = (x_center-ROI_radius_voxels):(x_center+ROI_radius_voxels);
        ygv = (y_center-ROI_radius_voxels):(y_center+ROI_radius_voxels);
        zgv = (z_center-ROI_radius_voxels):(z_center+ROI_radius_voxels);

        [ROI(iiROI).X,ROI(iiROI).Y,ROI(iiROI).Z] = meshgrid(xgv,ygv,zgv);
        
        nVoxels = length(find(ROI(iiROI).X));
        
        ROI(iiROI).idx_voxels = zeros(1,nVoxels);
        
        %nVoxels = 25;
        for iVoxel=1:nVoxels
            
            ROI(iiROI).idx_voxels(iVoxel) = sub2ind(MNI_dim,ROI(iiROI).X(iVoxel),ROI(iiROI).Y(iVoxel),ROI(iiROI).Z(iVoxel));
            
        end
    
    end
    
end

all_idx_voxels = zeros(1,length(ROI)*nVoxels);
iiVoxel = 0;
for iiROI=1:length(ROI)
    
    for iVoxel=1:nVoxels
        
        iiVoxel = iiVoxel + 1;
        
        all_idx_voxels(iiVoxel) = ROI(iiROI).idx_voxels(iVoxel);
        
    end
    
end

MOT4Run1_voxels = zeros(length(all_idx_voxels),nTR);
MOT4Run2_voxels = zeros(length(all_idx_voxels),nTR);
RestingStateRun1_voxels = zeros(length(all_idx_voxels),nTR);
RestingStateRun2_voxels = zeros(length(all_idx_voxels),nTR);

for iiVoxel=1:length(all_idx_voxels)
    
    [idxx,idxy,idxz] = ind2sub(MNI_dim,all_idx_voxels(iiVoxel));
    
    MOT4Run1_voxels(iiVoxel,:) = squeeze(MOT4Run1(idxx,idxy,idxz,:));
    MOT4Run2_voxels(iiVoxel,:) = squeeze(MOT4Run2(idxx,idxy,idxz,:));
    RestingStateRun1_voxels(iiVoxel,:) = squeeze(RestingStateRun1(idxx,idxy,idxz,:));
    RestingStateRun2_voxels(iiVoxel,:) = squeeze(RestingStateRun2(idxx,idxy,idxz,:));
    
end
    
[MOT4Run1_GF,MOT4Run1_Gpval, GSig, morder, A, G, info] = getGrangerX(MOT4Run1_voxels); 

[MOT4Run2_GF,MOT4Run2_Gpval, GSig, morder, A, G, info] = getGrangerX(MOT4Run2_voxels); 

[RestingStateRun1_GF,RestingStateRun1_Gpval, GSig, morder, A, G, info] = getGrangerX(RestingStateRun1_voxels); 

[RestingStateRun2_GF,RestingStateRun2_Gpval, GSig, morder, A, G, info] = getGrangerX(RestingStateRun2_voxels); 

save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'.mat'),'all_networks_labels','ROI','all_idx_voxels','MOT4Run1_GF','MOT4Run2_GF','RestingStateRun1_GF','RestingStateRun2_GF','MOT4Run1_Gpval','MOT4Run2_Gpval','RestingStateRun1_Gpval','RestingStateRun2_Gpval');

end

function plotGrangerForFunctionalNetworksOnlyOnSeeds(settings)

analysis_label = 'Granger-Func-Seeds';

pcriterion = 0.01;

load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'.mat'));

nROI = length(ROI);
nVoxels = length(ROI(1).idx_voxels);
%nVoxels = 25;

MOT4Run1_fractions = zeros(nROI,nROI);
MOT4Run2_fractions = zeros(nROI,nROI);
RestingStateRun1_fractions = zeros(nROI,nROI);
RestingStateRun2_fractions = zeros(nROI,nROI);

for iROI=1:nROI
    
    start_line = (nVoxels*(iROI-1)) + 1;
    end_line =  (nVoxels*(iROI-1)) + nVoxels;
    
    for iiROI=1:nROI
        
        start_column = (nVoxels*(iiROI-1)) + 1;
        end_column =  (nVoxels*(iiROI-1)) + nVoxels;

        pvalues = MOT4Run1_Gpval(start_line:end_line,start_column:end_column);
        psig = length(find(pvalues<pcriterion))/(nVoxels*nVoxels);
        MOT4Run1_fractions(iROI,iiROI) = psig;
        
        pvalues = MOT4Run2_Gpval(start_line:end_line,start_column:end_column);
        psig = length(find(pvalues<pcriterion))/(nVoxels*nVoxels);
        MOT4Run2_fractions(iROI,iiROI) = psig;
        
        pvalues = RestingStateRun1_Gpval(start_line:end_line,start_column:end_column);
        psig = length(find(pvalues<pcriterion))/(nVoxels*nVoxels);
        RestingStateRun1_fractions(iROI,iiROI) = psig;
        
        pvalues = RestingStateRun2_Gpval(start_line:end_line,start_column:end_column);
        psig = length(find(pvalues<pcriterion))/(nVoxels*nVoxels);
        RestingStateRun2_fractions(iROI,iiROI) = psig;
        
    end
    
end

HighRestingStateRun1_fractions = MOT4Run1_fractions - RestingStateRun1_fractions;
HighRestingStateRun2_fractions = MOT4Run2_fractions - RestingStateRun2_fractions;
        
clrmp = colormap('jet');
clrmp(1,:) = [1 1 1];

label = 'High-Run1-Fractions-GF';
plotFractionsGFSignificant(settings,MOT4Run1_fractions,label,analysis_label,ROI,[0 0.4],clrmp);

% label = 'High-Run2-Fractions-GF';
% plotFractionsGFSignificant(settings,MOT4Run2_fractions,label,analysis_label,ROI,[0 0.8],clrmp);
% 
% label = 'RestingState-Run1-Fractions-GF';
% plotFractionsGFSignificant(settings,RestingStateRun1_fractions,label,analysis_label,ROI,[0 0.8],clrmp);
% 
% label = 'RestingState-Run2-Fractions-GF';
% plotFractionsGFSignificant(settings,RestingStateRun2_fractions,label,analysis_label,ROI,[0 0.8],clrmp);
% 
% clrmp = colormap('jet');
% 
% label = 'High-RestingState-Run1-Fractions-GF';
% plotFractionsGFSignificant(settings,HighRestingStateRun1_fractions,label,analysis_label,ROI,[-1 1],clrmp);
% 
% label = 'High-RestingState-Run2-Fractions-GF';
% plotFractionsGFSignificant(settings,HighRestingStateRun2_fractions,label,analysis_label,ROI,[-1 1],clrmp);

end

function plotGrangerForFunctionalNetworksOnlyOnSeedsOnlyNetworks(settings)

analysis_label = 'Granger-Func-Seeds';

pcriterion = 0.01;

load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'.mat'));

nROI = length(ROI);
nVoxels = length(ROI(1).idx_voxels);
%nVoxels = 25;

nNet = length(all_networks_labels);

MOT4Run1_fractions = zeros(nROI,nROI);
MOT4Run2_fractions = zeros(nROI,nROI);
RestingStateRun1_fractions = zeros(nROI,nROI);
RestingStateRun2_fractions = zeros(nROI,nROI);

for iROI=1:nROI
    
    start_line = (nVoxels*(iROI-1)) + 1;
    end_line =  (nVoxels*(iROI-1)) + nVoxels;
    
    for iiROI=1:nROI
        
        start_column = (nVoxels*(iiROI-1)) + 1;
        end_column =  (nVoxels*(iiROI-1)) + nVoxels;

        pvalues = MOT4Run1_Gpval(start_line:end_line,start_column:end_column);
        psig = length(find(pvalues<pcriterion))/(nVoxels*nVoxels);
        MOT4Run1_fractions(iROI,iiROI) = psig;
        
        pvalues = MOT4Run2_Gpval(start_line:end_line,start_column:end_column);
        psig = length(find(pvalues<pcriterion))/(nVoxels*nVoxels);
        MOT4Run2_fractions(iROI,iiROI) = psig;
        
        pvalues = RestingStateRun1_Gpval(start_line:end_line,start_column:end_column);
        psig = length(find(pvalues<pcriterion))/(nVoxels*nVoxels);
        RestingStateRun1_fractions(iROI,iiROI) = psig;
        
        pvalues = RestingStateRun2_Gpval(start_line:end_line,start_column:end_column);
        psig = length(find(pvalues<pcriterion))/(nVoxels*nVoxels);
        RestingStateRun2_fractions(iROI,iiROI) = psig;
        
    end
    
end

for iNet=1:nNet
    
    iiROI = 0;
    
    for iROI=1:nROI
       
        if strcmp(ROI(iROI).netLabel,all_networks_labels{iNet});
        
            iiROI = iiROI + 1;
            ROI_n(iNet).all_ROIs(iiROI) = iROI;
            
        end
        
    end
    
end

for iNet=1:nNet
    
    ROI_n(iNet).min_ROI = min(ROI_n(iNet).all_ROIs);
    ROI_n(iNet).max_ROI = max(ROI_n(iNet).all_ROIs);

end

for iNet=1:nNet
    
    for iiNet=1:nNet
        
        MOT4Run1_meanfractions(iNet,iiNet) = mean(mean(MOT4Run1_fractions(ROI_n(iNet).min_ROI:ROI_n(iNet).max_ROI,ROI_n(iiNet).min_ROI:ROI_n(iiNet).max_ROI)));
        MOT4Run2_meanfractions(iNet,iiNet) = mean(mean(MOT4Run2_fractions(ROI_n(iNet).min_ROI:ROI_n(iNet).max_ROI,ROI_n(iiNet).min_ROI:ROI_n(iiNet).max_ROI)));
        RestingStateRun1_meanfractions(iNet,iiNet) = mean(mean(RestingStateRun1_fractions(ROI_n(iNet).min_ROI:ROI_n(iNet).max_ROI,ROI_n(iiNet).min_ROI:ROI_n(iiNet).max_ROI)));
        RestingStateRun2_meanfractions(iNet,iiNet) = mean(mean(RestingStateRun2_fractions(ROI_n(iNet).min_ROI:ROI_n(iNet).max_ROI,ROI_n(iiNet).min_ROI:ROI_n(iiNet).max_ROI)));

    end
    
end

clrmp = colormap('jet');
clrmp(1,:) = [1 1 1];
CLIM = [0 0.15];

label = 'High-Run1-MeanFractions-GF';
plotMeanFractionsGFSignificant(settings,MOT4Run1_meanfractions,label,analysis_label,all_networks_labels,CLIM,clrmp);

% label = 'High-Run2-MeanFractions-GF';
% plotMeanFractionsGFSignificant(settings,MOT4Run2_meanfractions,label,analysis_label,all_networks_labels,CLIM,clrmp);
% 
% label = 'Rest-Run1-MeanFractions-GF';
% plotMeanFractionsGFSignificant(settings,RestingStateRun1_meanfractions,label,analysis_label,all_networks_labels,CLIM,clrmp);
% 
% label = 'Rest-Run2-MeanFractions-GF';
% plotMeanFractionsGFSignificant(settings,RestingStateRun2_meanfractions,label,analysis_label,all_networks_labels,CLIM,clrmp);

end

function plotMeanFractionsGFSignificant(settings,fractions_mat,label,analysis_label,ROI_labels,CLIM,clrmp)

nROI = length(ROI_labels);

rmin = CLIM(1);
rmax =  CLIM(2);

Tick = 1:1:nROI;
for iROI=1:nROI
    TickLabel{iROI} = ROI_labels{iROI};
end

f = figure;

hold on;
caxis([rmin rmax]);
h = imagesc( fractions_mat );
% set( h, 'EdgeColor', 'none');
colormap( clrmp);
% hold off;
% axis 'square'; 
xlabel('To');
ylabel('From');
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

xticklabel_rotate;

title(label);

h = colorbar;
set(h,'YLim',[rmin rmax],'YTick',[rmin rmax]);

print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-',label,'.jpg'));
print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-',label,'.eps'));
print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-',label,'.pdf'));


end

function plotFractionsGFSignificant(settings,fractions_mat,label,analysis_label,ROI,CLIM,clrmp)

nROI = length(ROI);

rmin = CLIM(1);
rmax =  CLIM(2);

Tick = 1:nROI;
for iROI=1:nROI
    TickLabel{iROI} = strcat(ROI(iROI).netLabel,'-',ROI(iROI).seedLabel);
end

f = figure;

hold on;
caxis([rmin rmax]);
h = pcolor( fractions_mat );
set( h, 'EdgeColor', 'none');
colormap( clrmp);
hold off;
axis 'square'; 
xlabel('To');
ylabel('From');
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );

xticklabel_rotate;

title(label);

h = colorbar;
set(h,'YLim',[rmin rmax],'YTick',[rmin rmax]);

print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-',label,'.jpg'));
print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-',label,'.eps'));
print(f,'-dpdf',strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-',label,'.pdf'));


end


