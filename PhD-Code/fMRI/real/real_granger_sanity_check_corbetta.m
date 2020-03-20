
function real_granger_sanity_check_corbetta

% extractVoxelsIndexes;

% computeGrangerBtwVoxels;

% computeMeanGrangerConsistency;

% plotConsistency;

% computeGrangerBtwMeanOfVoxels;

% computeMeanGrangerZScore;

plotZscore;

end


function extractVoxelsIndexes

%%% FUNCTIONAL NETWORKS LABELS

Corbetta.Net(1).label = 'DAN';
Corbetta.Net(2).label = 'VAN';
Corbetta.Net(3).label = 'VIS';

%%% SEEDS LABELS

Corbetta.Net(1).Seeds_Label{1} = 'L-dFEF';
Corbetta.Net(1).Seeds_Label{2} = 'R-dFEF';
Corbetta.Net(1).Seeds_Label{3} = 'L-pIPS';
Corbetta.Net(1).Seeds_Label{4} = 'R-pIPS';
Corbetta.Net(1).Seeds_Label{5} = 'L-SPL';
Corbetta.Net(1).Seeds_Label{6} = 'R-SPL';

Corbetta.Net(2).Seeds_Label{1} = 'PreCu';
Corbetta.Net(2).Seeds_Label{2} = 'R-vTPJ';

Corbetta.Net(3).Seeds_Label{1} = 'L-MT';
Corbetta.Net(3).Seeds_Label{2} = 'R-MT';
Corbetta.Net(3).Seeds_Label{3} = 'L-V3AV7';
Corbetta.Net(3).Seeds_Label{4} = 'R-V3AV7';
Corbetta.Net(3).Seeds_Label{5} = 'L-V4V8';
Corbetta.Net(3).Seeds_Label{6} = 'R-V4V8';

%%% SEEDS CORRESPONDENCE

% v1

% Corbetta.Net(1).Seeds_Correspondence(1).idx = [4 5 6];
% Corbetta.Net(1).Seeds_Correspondence(2).idx = [1 2 3];
% Corbetta.Net(1).Seeds_Correspondence(3).idx = [12 13 14];
% Corbetta.Net(1).Seeds_Correspondence(4).idx = [9 10 11];
% Corbetta.Net(1).Seeds_Correspondence(5).idx = [19];
% Corbetta.Net(1).Seeds_Correspondence(6).idx = [18];
% 
% Corbetta.Net(2).Seeds_Correspondence(1).idx = [10];
% Corbetta.Net(2).Seeds_Correspondence(2).idx = [3];
% 
% Corbetta.Net(3).Seeds_Correspondence(1).idx = [16 17];
% Corbetta.Net(3).Seeds_Correspondence(2).idx = [13 14 15];
% Corbetta.Net(3).Seeds_Correspondence(3).idx = [6 10];
% Corbetta.Net(3).Seeds_Correspondence(4).idx = [5 9];
% Corbetta.Net(3).Seeds_Correspondence(5).idx = [8];
% Corbetta.Net(3).Seeds_Correspondence(6).idx = [7];

% v2

Corbetta.Net(1).Seeds_Correspondence(1).idx = [4];
Corbetta.Net(1).Seeds_Correspondence(2).idx = [1];
Corbetta.Net(1).Seeds_Correspondence(3).idx = [12];
Corbetta.Net(1).Seeds_Correspondence(4).idx = [9];
Corbetta.Net(1).Seeds_Correspondence(5).idx = [19];
Corbetta.Net(1).Seeds_Correspondence(6).idx = [18];

Corbetta.Net(2).Seeds_Correspondence(1).idx = [10];
Corbetta.Net(2).Seeds_Correspondence(2).idx = [3];

Corbetta.Net(3).Seeds_Correspondence(1).idx = [16];
Corbetta.Net(3).Seeds_Correspondence(2).idx = [13];
Corbetta.Net(3).Seeds_Correspondence(3).idx = [6];
Corbetta.Net(3).Seeds_Correspondence(4).idx = [5];
Corbetta.Net(3).Seeds_Correspondence(5).idx = [8];
Corbetta.Net(3).Seeds_Correspondence(6).idx = [7];

idx_all_seeds_voxels = [];

for iNet=1:length(Corbetta.Net)
   
    seeds = getFunctionalSeeds_v6(Corbetta.Net(iNet).label);
    
    for iSeed=1:length(Corbetta.Net(iNet).Seeds_Correspondence)
        
        idx_voxels = [];
       
        for iiSeed=1:length(Corbetta.Net(iNet).Seeds_Correspondence(iSeed).idx)
            
            x = seeds.ROI(Corbetta.Net(iNet).Seeds_Correspondence(iSeed).idx(iiSeed)).x;
            y = seeds.ROI(Corbetta.Net(iNet).Seeds_Correspondence(iSeed).idx(iiSeed)).y;
            z = seeds.ROI(Corbetta.Net(iNet).Seeds_Correspondence(iSeed).idx(iiSeed)).z;
            
            MNI_dimension = [91 109 91];
            
            MNI_x_center = 45;
            MNI_y_center = 63;
            MNI_z_center = 36;

            size_voxels_mm = 2;
            ROI_radius_mm = 4 - size_voxels_mm;
            ROI_radius_voxels = ROI_radius_mm/size_voxels_mm;

            x_center = MNI_x_center + round(x/size_voxels_mm) + 1;
            y_center = MNI_y_center + round(y/size_voxels_mm) + 1;
            z_center = MNI_z_center + round(z/size_voxels_mm) + 1;

            xgv = (x_center-ROI_radius_voxels):(x_center+ROI_radius_voxels);
            ygv = (y_center-ROI_radius_voxels):(y_center+ROI_radius_voxels);
            zgv = (z_center-ROI_radius_voxels):(z_center+ROI_radius_voxels);

            [X,Y,Z] = meshgrid(xgv,ygv,zgv);
            
            tmp_idx = find(X);
            nVoxels = length(tmp_idx);
            
            for iVoxel=1:nVoxels

                idx_voxels = [idx_voxels, sub2ind(MNI_dimension,X(iVoxel),Y(iVoxel),Z(iVoxel))];
                
            end

        end
        
        Corbetta.Net(iNet).Seeds_Correspondence(iSeed).idx_voxels = idx_voxels;
        
    end
    
end

for iNet=1:length(Corbetta.Net)
    for iSeed=1:length(Corbetta.Net(iNet).Seeds_Correspondence)
        idx_all_seeds_voxels = [idx_all_seeds_voxels, Corbetta.Net(iNet).Seeds_Correspondence(iSeed).idx_voxels];
    end
end

Corbetta.idx_all_seeds_voxels = idx_all_seeds_voxels;

save('Corbetta-Granger-Sanity-Check-Seeds-Idx-Voxels.mat','Corbetta');


end

function computeGrangerBtwVoxels

load('Corbetta-Granger-Sanity-Check-Seeds-Idx-Voxels.mat');

idx_voxels = Corbetta.idx_all_seeds_voxels;

all_settings = getAllSettings;

nRuns = 4;
nTR = 150;
MNI_size = [91 109 91];
iiRun_t = 0;
iiRun_r = 0;

for iSettings=1:length(all_settings)
    
    settings = all_settings(iSettings).settings;

    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;

    real_load_all_data_FSL;
    
    for iRun=1:nRuns
        
        iiRun_t = iiRun_t + 1;
        
        X = zeros(length(idx_voxels),nTR);
        
        start = now;
        disp(strcat('SUBJ:',int2str(iSettings),':Track:',int2str(iRun),datestr(start)));
        
        disp('...extracting voxels');
        
        for iVoxel=1:length(idx_voxels)
            
           [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
           
           X(iVoxel,:) = squeeze(Track(iRun).run(idxx,idxy,idxz,:));
            
        end
        
        disp('...computing Granger');
        
        GC.track.GF = zeros(length(idx_voxels),length(idx_voxels));
        GC.track.Gpval = zeros(length(idx_voxels),length(idx_voxels));
        icomputation = 0;
        for iVoxel=1:length(idx_voxels)
            
            Voxel1 = squeeze(X(iVoxel,:));
            
            for iiVoxel=iVoxel:length(idx_voxels)
                
                Voxel2 = squeeze(X(iiVoxel,:));
        
                icomputation = icomputation + 1;
                if mod(icomputation,100) == 0; disp(strcat('...:',int2str(icomputation))); end
                
                [GF, Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(Voxel1,Voxel2);
        
                GC.track.GF(iiVoxel,iVoxel) = GF(2,1);
                GC.track.GF(iVoxel,iiVoxel) = GF(1,2);
                GC.track.Gpval(iiVoxel,iVoxel) = Gpval(2,1);
                GC.track.Gpval(iVoxel,iiVoxel) = Gpval(1,2);
            
            end
            
        end
        
        finish = now;
        disp(strcat('lasted:',datestr(finish-start)));
        
        save(strcat('Corbetta-Granger-Sanity-Check-Seeds-GC-Track-Run-',int2str(iiRun_t),'.mat'),'GC','iiRun_t');
        clear GC
 
    end
    
    for iRun=1:nRuns
        
        iiRun_r = iiRun_r + 1;
        
        X = zeros(length(idx_voxels),nTR);
        
        start = now;
        disp(strcat('SUBJ:',int2str(iSettings),':RestingState:',int2str(iRun),datestr(start)));
        
        disp('...extracting voxels');

        for iVoxel=1:length(idx_voxels)
            
           [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
           
           X(iVoxel,:) = squeeze(RestingState(iRun).run(idxx,idxy,idxz,:));
            
        end
        
        disp('...computing Granger');
        
        GC.rest.GF = zeros(length(idx_voxels),length(idx_voxels));
        GC.rest.Gpval = zeros(length(idx_voxels),length(idx_voxels));
        icomputation = 0;
        for iVoxel=1:length(idx_voxels)
            
            if mod(iVoxel,100) == 0; disp(strcat('...voxel:',int2str(iVoxel))); end
            
            Voxel1 = squeeze(X(iVoxel,:));
            
            for iiVoxel=iVoxel:length(idx_voxels)
                
                Voxel2 = squeeze(X(iiVoxel,:));
        
                icomputation = icomputation + 1;
                if mod(icomputation,100) == 0; disp(strcat('...:',int2str(icomputation))); end
                
                [GF, Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(Voxel1,Voxel2);
        
                GC.rest.GF(iiVoxel,iVoxel) = GF(2,1);
                GC.rest.GF(iVoxel,iiVoxel) = GF(1,2);
                GC.rest.Gpval(iiVoxel,iVoxel) = Gpval(2,1);
                GC.rest.Gpval(iVoxel,iiVoxel) = Gpval(1,2);
            
            end
            
        end
        
        finish = now;
        disp(strcat('lasted:',datestr(finish-start)));
        
        save(strcat('Corbetta-Granger-Sanity-Check-Seeds-GC-Rest-Run-',int2str(iiRun_r),'.mat'),'GC','iiRun_r');
        clear GC
 
    end

    clear Track
    clear Passive
    clear RestingState

end

end

function computeMeanGrangerConsistency

load('Corbetta-Granger-Sanity-Check-Seeds-Idx-Voxels.mat');

nRuns = 32;
nNet = length(Corbetta.Net);
nVoxelsPerSeed = length(Corbetta.Net(1).Seeds_Correspondence(1).idx_voxels);
conditions = {'Track' 'Rest'};
nConditions = length(conditions);
pcriterion = 0.01;

nSeeds = 0;
for iNet=1:nNet
    nSeeds = nSeeds + length(Corbetta.Net(iNet).Seeds_Correspondence);
end

all_samples = zeros(nConditions,nRuns,nSeeds,nSeeds);

for iCondition=1:nConditions
    
   disp(conditions{iCondition}); 
   
   for iRun=1:nRuns
   
        load(strcat('Corbetta-Granger-Sanity-Check-Seeds-GC-',conditions{iCondition},'-Run-',int2str(iRun),'.mat'));
        
        eval(strcat('GC.',lower(conditions{iCondition}),'.Gpval(isnan(GC.',lower(conditions{iCondition}),'.Gpval))=0'));
        
        for iSeed=1:nSeeds
           
            for iiSeed=1:nSeeds
                
                all_samples(iCondition,iRun,iSeed,iiSeed) = eval(strcat('length(find(GC.',lower(conditions{iCondition}),'.Gpval(((iSeed-1)*nVoxelsPerSeed+1):((iSeed-1)*nVoxelsPerSeed+nVoxelsPerSeed),((iiSeed-1)*nVoxelsPerSeed+1):((iiSeed-1)*nVoxelsPerSeed+nVoxelsPerSeed))<pcriterion)) / (nVoxelsPerSeed^2)'));
                
            end
            
        end
        
        clear GC
        clear iiRun_r
        clear iiRun_t
   
   end
    
end

track = squeeze(all_samples(1,:,:,:));
m_track = squeeze(mean(track,1));

rest = squeeze(all_samples(2,:,:,:));
m_rest = squeeze(mean(rest,1));

contrast_val = m_track - m_rest;

for iSeed=1:nSeeds
    
    for iiSeed=1:nSeeds
        
        for iRun=1:nRuns
            
            pair_track(iRun) = track(iRun,iSeed,iiSeed);
            pair_rest(iRun) = rest(iRun,iSeed,iiSeed);
            
        end
        
        contrast_sig(iSeed,iiSeed) = ttest(pair_track(:),pair_rest(:));
        
    end
    
end

save(strcat('Corbetta-Granger-Sanity-Check-Seeds-GC-Consistency','.mat'),'all_samples','m_track','m_rest','contrast_val','contrast_sig');

end

function plotConsistency

fs = 10;

load('Corbetta-Granger-Sanity-Check-Seeds-GC-Consistency.mat');
load('Corbetta-Granger-Sanity-Check-Seeds-Idx-Voxels.mat');

iiSeed = 0;
for iNet=1:length(Corbetta.Net)
    
    for iSeed=1:length(Corbetta.Net(iNet).Seeds_Correspondence)
        
        iiSeed = iiSeed + 1;
        seeds_label{iiSeed} = Corbetta.Net(iNet).Seeds_Label{iSeed};
        
    end
    
end

clrmp = colormap('jet');
min_C = min([m_track(:); m_rest(:)]);
max_C = max([m_track(:); m_rest(:)]);

figure;
imagesc(m_track);
caxis([min_C max_C]);
colormap(clrmp);
colorbar;

ax = gca;
set(ax,'XTick',1:length(seeds_label));
set(ax,'YTick',1:length(seeds_label));
set(ax,'XTickLabel',seeds_label);
set(ax,'YTickLabel',seeds_label);
set(ax,'FontSize',fs);

xticklabel_rotate([],90,[],'Fontsize',fs);

hold on
plot([0.5 iiSeed+0.5],[6.5 6.5],'k');
plot([0.5 iiSeed+0.5],[8.5 8.5],'k');
plot([6.5 6.5],[0.5 iiSeed+0.5],'k');
plot([8.5 8.5],[0.5 iiSeed+0.5],'k');

figure;
imagesc(m_rest);
caxis([min_C max_C]);
colormap(clrmp);
colorbar;

ax = gca;
set(ax,'XTick',1:length(seeds_label));
set(ax,'YTick',1:length(seeds_label));
set(ax,'XTickLabel',seeds_label);
set(ax,'YTickLabel',seeds_label);
set(ax,'FontSize',fs);

xticklabel_rotate([],90,[],'Fontsize',fs);

hold on
plot([0.5 iiSeed+0.5],[6.5 6.5],'k');
plot([0.5 iiSeed+0.5],[8.5 8.5],'k');
plot([6.5 6.5],[0.5 iiSeed+0.5],'k');
plot([8.5 8.5],[0.5 iiSeed+0.5],'k');

figure;
imagesc(contrast_val);
min_C = min(contrast_val(:));
max_C = max(contrast_val(:));
min_C = min([-1*abs(min_C) -1*abs(max_C)]);
max_C = max([abs(min_C) abs(max_C)]);
caxis([min_C max_C]);
colormap(clrmp);
colorbar;

ax = gca;
set(ax,'XTick',1:length(seeds_label));
set(ax,'YTick',1:length(seeds_label));
set(ax,'XTickLabel',seeds_label);
set(ax,'YTickLabel',seeds_label);
set(ax,'FontSize',fs);

xticklabel_rotate([],90,[],'Fontsize',fs);

hold on
plot([0.5 iiSeed+0.5],[6.5 6.5],'k');
plot([0.5 iiSeed+0.5],[8.5 8.5],'k');
plot([6.5 6.5],[0.5 iiSeed+0.5],'k');
plot([8.5 8.5],[0.5 iiSeed+0.5],'k');

for iSeed=1:iiSeed
   
    for iiiSeed=1:iiSeed

        if contrast_sig(iSeed,iiiSeed)
           
            plot(iSeed,iiiSeed,'k*','MarkerSize',15);
            
        end
        
    end
    
end

end

function computeGrangerBtwMeanOfVoxels

load('Corbetta-Granger-Sanity-Check-Seeds-Idx-Voxels.mat');

all_settings = getAllSettings;

nRuns = 4;
nTR = 150;
MNI_size = [91 109 91];
iiRun_t = 0;
iiRun_r = 0;
nNets = length(Corbetta.Net);
nSeeds = 14;
nVoxelsPerSeed = 27;

for iSettings=1:length(all_settings)
    
    settings = all_settings(iSettings).settings;

    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;

    real_load_all_data_FSL;
    
    for iRun=1:nRuns
        
        iiRun_t = iiRun_t + 1;
  
        start = now;
        disp(strcat('SUBJ:',int2str(iSettings),':Track:',int2str(iRun),datestr(start)));
        
        disp('...extracting voxels');
        
        iiSeed = 0;
        for iNet=1:nNets
            
            for iSeed=1:length(Corbetta.Net(iNet).Seeds_Correspondence)
                
                iiSeed = iiSeed + 1;
                
                idx_voxels = Corbetta.Net(iNet).Seeds_Correspondence(iSeed).idx_voxels;
                
                X = zeros(nVoxelsPerSeed,nTR);
            
                for iVoxel=1:length(idx_voxels)
            
                    [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
           
                    X(iVoxel,:) = squeeze(Track(iRun).run(idxx,idxy,idxz,:));
            
                end
                
                m_X = mean(X,1);
                
                Seeds_Mean(iiSeed).mean = m_X;
                
            end
        
        end
        
        disp('...computing Granger');
        
        GC.track.GF = zeros(nSeeds,nSeeds);
        GC.track.Gpval = zeros(nSeeds,nSeeds);
        
        for iSeed=1:nSeeds
            
            seed_ts1 = Seeds_Mean(iSeed).mean;
            
            for iiSeed=iSeed:nSeeds
                
                seed_ts2 = Seeds_Mean(iiSeed).mean;
                
                [GF, Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(seed_ts1,seed_ts2);
        
                GC.track.GF(iiSeed,iSeed) = GF(2,1);
                GC.track.GF(iSeed,iiSeed) = GF(1,2);
                GC.track.Gpval(iiSeed,iSeed) = Gpval(2,1);
                GC.track.Gpval(iSeed,iiSeed) = Gpval(1,2);
            
            end
            
        end
        
        finish = now;
        disp(strcat('lasted:',datestr(finish-start)));
        
        save(strcat('Corbetta-Granger-Sanity-Check-Seeds-GC-Mean-Track-Run-',int2str(iiRun_t),'.mat'),'GC','Seeds_Mean','iiRun_t');
        clear GC
 
    end
    
    for iRun=1:nRuns
        
        iiRun_r = iiRun_r + 1;

        start = now;
        disp(strcat('SUBJ:',int2str(iSettings),':RestingState:',int2str(iRun),datestr(start)));
        
        disp('...extracting voxels');
        
        iiSeed = 0;
        for iNet=1:nNets
            
            for iSeed=1:length(Corbetta.Net(iNet).Seeds_Correspondence)
                
                iiSeed = iiSeed + 1;
                
                idx_voxels = Corbetta.Net(iNet).Seeds_Correspondence(iSeed).idx_voxels;
                
                X = zeros(nVoxelsPerSeed,nTR);
            
                for iVoxel=1:length(idx_voxels)
            
                    [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
           
                    X(iVoxel,:) = squeeze(RestingState(iRun).run(idxx,idxy,idxz,:));
            
                end
                
                m_X = mean(X,1);
                
                Seeds_Mean(iiSeed).mean = m_X;
                
            end
        
        end
        
        disp('...computing Granger');
        
        GC.rest.GF = zeros(nSeeds,nSeeds);
        GC.rest.Gpval = zeros(nSeeds,nSeeds);
        
        for iSeed=1:nSeeds
            
            seed_ts1 = Seeds_Mean(iSeed).mean;
            
            for iiSeed=iSeed:nSeeds
                
                seed_ts2 = Seeds_Mean(iiSeed).mean;
                
                [GF, Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(seed_ts1,seed_ts2);
        
                GC.rest.GF(iiSeed,iSeed) = GF(2,1);
                GC.rest.GF(iSeed,iiSeed) = GF(1,2);
                GC.rest.Gpval(iiSeed,iSeed) = Gpval(2,1);
                GC.rest.Gpval(iSeed,iiSeed) = Gpval(1,2);
            
            end
            
        end
        
        finish = now;
        disp(strcat('lasted:',datestr(finish-start)));
        
        save(strcat('Corbetta-Granger-Sanity-Check-Seeds-GC-Mean-Rest-Run-',int2str(iiRun_r),'.mat'),'GC','Seeds_Mean','iiRun_r');
        clear GC
 
    end
    
    clear Track
    clear Passive
    clear RestingState

end

end

function computeMeanGrangerZScore

load('Corbetta-Granger-Sanity-Check-Seeds-Idx-Voxels.mat');

nRuns = 32;
nNet = length(Corbetta.Net);
nVoxelsPerSeed = length(Corbetta.Net(1).Seeds_Correspondence(1).idx_voxels);
conditions = {'Track' 'Rest'};
nConditions = length(conditions);
pcriterion = 0.01;

nSeeds = 0;
for iNet=1:nNet
    nSeeds = nSeeds + length(Corbetta.Net(iNet).Seeds_Correspondence);
end

for iCondition=1:nConditions
    
   disp(conditions{iCondition}); 
   
   for iRun=1:nRuns
   
        load(strcat('Corbetta-Granger-Sanity-Check-Seeds-GC-Mean-',conditions{iCondition},'-Run-',int2str(iRun),'.mat'));
        
        eval(strcat('GC.',lower(conditions{iCondition}),'.Gpval(isnan(GC.',lower(conditions{iCondition}),'.Gpval))=0'));
        
        for iSeed=1:nSeeds
           
            for iiSeed=1:nSeeds
                
                all_samples(iCondition,iRun,iSeed,iiSeed) = eval(strcat('GC.',lower(conditions{iCondition}),'.Gpval(iSeed,iiSeed);'));
                
            end
            
        end
        
        clear GC
        clear iiRun_r
        clear iiRun_t
   
   end
    
end

track = squeeze(all_samples(1,:,:,:));
track = -1*norminv(track);
m_track = squeeze(mean(track,1));

rest = squeeze(all_samples(2,:,:,:));
rest = -1*norminv(rest);
m_rest = squeeze(mean(rest,1));

contrast_val = m_track - m_rest;

for iSeed=1:nSeeds
    
    for iiSeed=1:nSeeds
        
        for iRun=1:nRuns
            
            pair_track(iRun) = track(iRun,iSeed,iiSeed);
            pair_rest(iRun) = rest(iRun,iSeed,iiSeed);
            
        end
        
        contrast_sig(iSeed,iiSeed) = ttest(pair_track(:),pair_rest(:));
        
    end
    
end

save(strcat('Corbetta-Granger-Sanity-Check-Seeds-GC-Mean-Consistency','.mat'),'all_samples','m_track','m_rest','contrast_val','contrast_sig');

end

function plotZscore

fs = 10;

load('Corbetta-Granger-Sanity-Check-Seeds-GC-Mean-Consistency.mat');
load('Corbetta-Granger-Sanity-Check-Seeds-Idx-Voxels.mat');

iiSeed = 0;
for iNet=1:length(Corbetta.Net)
    
    for iSeed=1:length(Corbetta.Net(iNet).Seeds_Correspondence)
        
        iiSeed = iiSeed + 1;
        seeds_label{iiSeed} = Corbetta.Net(iNet).Seeds_Label{iSeed};
        
    end
    
end

m_track(isnan(m_track)) = 0;
m_rest(isnan(m_rest)) = 0;

m_track(isinf(m_track)) = 0;
m_rest(isinf(m_rest)) = 0;

clrmp = colormap('jet');
min_C = min([m_track(:); m_rest(:)]);
max_C = max([m_track(:); m_rest(:)]);

figure;
imagesc(m_track);
caxis([min_C max_C]);
colormap(clrmp);
colorbar;

ax = gca;
set(ax,'XTick',1:length(seeds_label));
set(ax,'YTick',1:length(seeds_label));
set(ax,'XTickLabel',seeds_label);
set(ax,'YTickLabel',seeds_label);
set(ax,'FontSize',fs);

xticklabel_rotate([],90,[],'Fontsize',fs);

hold on
plot([0.5 iiSeed+0.5],[6.5 6.5],'k');
plot([0.5 iiSeed+0.5],[8.5 8.5],'k');
plot([6.5 6.5],[0.5 iiSeed+0.5],'k');
plot([8.5 8.5],[0.5 iiSeed+0.5],'k');

figure;
imagesc(m_rest);
caxis([min_C max_C]);
colormap(clrmp);
colorbar;

ax = gca;
set(ax,'XTick',1:length(seeds_label));
set(ax,'YTick',1:length(seeds_label));
set(ax,'XTickLabel',seeds_label);
set(ax,'YTickLabel',seeds_label);
set(ax,'FontSize',fs);

xticklabel_rotate([],90,[],'Fontsize',fs);

hold on
plot([0.5 iiSeed+0.5],[6.5 6.5],'k');
plot([0.5 iiSeed+0.5],[8.5 8.5],'k');
plot([6.5 6.5],[0.5 iiSeed+0.5],'k');
plot([8.5 8.5],[0.5 iiSeed+0.5],'k');

figure;
contrast_val(isnan(contrast_val)) = 0;
contrast_val(isinf(contrast_val)) = 0;
imagesc(contrast_val);
min_C = min(contrast_val(:));
max_C = max(contrast_val(:));
min_C = min([-1*abs(min_C) -1*abs(max_C)]);
max_C = max([abs(min_C) abs(max_C)]);
caxis([min_C max_C]);
colormap(clrmp);
colorbar;

ax = gca;
set(ax,'XTick',1:length(seeds_label));
set(ax,'YTick',1:length(seeds_label));
set(ax,'XTickLabel',seeds_label);
set(ax,'YTickLabel',seeds_label);
set(ax,'FontSize',fs);

xticklabel_rotate([],90,[],'Fontsize',fs);

hold on
plot([0.5 iiSeed+0.5],[6.5 6.5],'k');
plot([0.5 iiSeed+0.5],[8.5 8.5],'k');
plot([6.5 6.5],[0.5 iiSeed+0.5],'k');
plot([8.5 8.5],[0.5 iiSeed+0.5],'k');

contrast_sig(isnan(contrast_sig)) = 0;
for iSeed=1:iiSeed
   
    for iiiSeed=1:iiSeed

        if contrast_sig(iSeed,iiiSeed)
           
            plot(iSeed,iiiSeed,'k*','MarkerSize',15);
            
        end
        
    end
    
end

end