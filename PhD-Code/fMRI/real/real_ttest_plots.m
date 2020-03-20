
function real_ttest_plots

% getTTESTtsPerAAL;

% getTTESTtsPerSeedNetwork;

plotTTESTtsPerAAL;
 
% plotTTESTtsPerSeedNetwork;

end

function getTTESTtsPerAAL

[idx_attention, idx_resting] = getActivationIdx;

ROI = getROIIdxVoxels(idx_attention,idx_resting);

[ROI_mean_ts_att,ROI_mean_ts_rest] = getMeanTSAllSubjects(ROI);

save(strcat('LHR','-','All-Subjects','-','T-Test','-','Mean','Time-Series','.mat'),'idx_attention','idx_resting','ROI','ROI_mean_ts_att','ROI_mean_ts_rest');

end

function [idx_attention, idx_resting] = getActivationIdx

pcriterion = 1.6;

%%% LOAD T-TEST

load_img = nifti('ttest-activation.nii');
load_img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FEAT\FEAT-ATTENTION-RESTING\',load_img.dat.fname);

TTEST_img = load_img.dat(:,:,:);

idx_attention = find(TTEST_img > pcriterion);
idx_resting = find(TTEST_img < -pcriterion);

end

function ROI = getROIIdxVoxels(idx_attention,idx_resting)

load_img = nifti('ttest-activation.nii');
load_img.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FEAT\FEAT-ATTENTION-RESTING\',load_img.dat.fname);

TTEST_img = load_img.dat(:,:,:);

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

%%% CHECK WHICH REGIONS

ROI.att_exist = zeros(1,nNodes);
ROI.rest_exist = zeros(1,nNodes);

for iROI=1:nNodes
    
    label_ROI = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    idx_att_voxels = idx_attention(ismember(idx_attention,idx_voxels));
    
    ROI.label{iROI} = label_ROI;
    
    if ~isempty(idx_att_voxels)
        
        ROI.att_exist(iROI) = 1;
        ROI.voxels(iROI).idx_att = idx_att_voxels;
        ROI.att_nVoxels(iROI) = length(idx_att_voxels);
        ROI.max_zscore_att(iROI).zscore = max(TTEST_img(idx_att_voxels));
        
        att_voxels_zscore = TTEST_img(idx_att_voxels);
        [s,i] = sort(abs(att_voxels_zscore),'descend');
        
        idx = idx_att_voxels(i(1));
        [x,y,z] = ind2sub(size(AAL_img),idx);
        
        ROI.max_zscore_att(iROI).x = x;
        ROI.max_zscore_att(iROI).y = y;
        ROI.max_zscore_att(iROI).z = z;
        
        [s,ROI.max_zscore_att_sorted_idx] = sort(ROI.max_zscore_att(iROI).zscore,'descend');
          
    end
    
    [s,ROI.max_amount_att_sorted_idx] = sort(ROI.att_nVoxels,'descend');
        
    idx_rest_voxels = idx_resting(ismember(idx_resting,idx_voxels));
    
    if ~isempty(idx_rest_voxels)
        
        ROI.rest_exist(iROI) = 1;
        ROI.voxels(iROI).idx_rest = idx_rest_voxels;
        ROI.rest_nVoxels(iROI) = length(idx_rest_voxels);
        ROI.max_zscore_rest(iROI).zscore = max(TTEST_img(idx_rest_voxels));
        
        rest_voxels_zscore = TTEST_img(idx_rest_voxels);
        [s,i] = sort(abs(rest_voxels_zscore),'descend');
        
        idx = idx_rest_voxels(i(1));
        [x,y,z] = ind2sub(size(AAL_img),idx);
        
        ROI.max_zscore_rest(iROI).x = x;
        ROI.max_zscore_rest(iROI).y = y;
        ROI.max_zscore_rest(iROI).z = z;
        
        [s,ROI.max_zscore_rest_sorted_idx] = sort(ROI.max_zscore_rest(iROI).zscore,'ascend');
        
    end
    
    [s,ROI.max_amount_rest_sorted_idx] = sort(ROI.rest_nVoxels,'descend');
    
end

end

function [ROI_mean_ts_att,ROI_mean_ts_rest] = getMeanTSAllSubjects(ROI)

all_settings = getAllSettings;

nRuns = 4;
nSets = length(all_settings);
nROI = 90;

MNI_size = [91 109 91];

iirun = 0;
for iSet=1:nSets
    
    settings = all_settings(iSet).settings;

    %% LOAD DATA

    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;

    kind = 'Trials';
    for irun=1:nRuns;

        iirun = iirun + 1;

        [Trials(iirun).run, Trials(iirun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);

    end
    
end

nTR = 170;

ROI_mean_ts_att = zeros(nROI,nTR);

for iROI=1:nROI
   
    if ROI.att_exist(iROI)
        
        idx_att_voxels = ROI.voxels(iROI).idx_att;
        
        nVoxels = length(idx_att_voxels)*nSets*nRuns;
        
        this_roi_voxels_ts = zeros(nVoxels,nTR);
        
        iiVoxel = 0;
        for iVoxel=1:length(idx_att_voxels)
            
            [x,y,z] = ind2sub(MNI_size,idx_att_voxels(iVoxel));
            
            iirun = 0;
            for iSet=1:nSets
                
                for irun=1:nRuns
                    
                    iiVoxel = iiVoxel + 1;
                    
                    iirun = iirun + 1;
                    
                    if mod(iirun,2) == 0; phase = 17; else phase = 1; end
                    
                    this_roi_voxels_ts(iiVoxel,:) = squeeze(Trials(iirun).run(x,y,z,(phase+1):(nTR+phase)));
                    
                end
                
            end
            
        end
        
        ROI_mean_ts_att(iROI,:) = mean(this_roi_voxels_ts,1);
        
    end
    
end

ROI_mean_ts_rest = zeros(nROI,nTR);

for iROI=1:nROI
   
    if ROI.rest_exist(iROI)
        
        idx_rest_voxels = ROI.voxels(iROI).idx_rest;
        
        nVoxels = length(idx_rest_voxels)*nSets*nRuns;
        
        this_roi_voxels_ts = zeros(nVoxels,nTR);
        
        iiVoxel = 0;
        for iVoxel=1:length(idx_rest_voxels)
            
            [x,y,z] = ind2sub(MNI_size,idx_rest_voxels(iVoxel));
            
            iirun = 0;
            for iSet=1:nSets
                
                for irun=1:nRuns
                    
                    iiVoxel = iiVoxel + 1;
                    
                     iirun = iirun + 1;
                     
                     if mod(iirun,2) ~= 0; phase = 17; else phase = 1; end
                    
                    this_roi_voxels_ts(iiVoxel,:) = squeeze(Trials(iirun).run(x,y,z,(phase+1):(nTR+phase)));
                    
                end
                
            end
            
        end
        
        ROI_mean_ts_rest(iROI,:) = mean(this_roi_voxels_ts,1);
        
    end
    
end

end

function plotTTESTtsPerAAL

load('LHR-All-Subjects-T-Test-MeanTime-Series.mat');

nROI = 8;

idx_attention = ROI.max_amount_att_sorted_idx(1:nROI);
idx_resting = ROI.max_amount_rest_sorted_idx(1:nROI);

nTRs=170;

t=1:1:nTRs;
box=[ones(1,17),zeros(1,17),ones(1,17),zeros(1,17),ones(1,17),zeros(1,17),ones(1,17),zeros(1,17),ones(1,17),zeros(1,17)];

T0=0; n=4; lamda=0.5;
hrf=((t-T0).^(n-1)).*exp(-(t-T0)/lamda)/((lamda^n)*factorial(n-1));

B=conv(hrf,box)/10;

for iROI=1:nROI
    
    %%% ATTENTION
    
    condition = 'Attention';
    
    area_label = ROI.label{idx_attention(iROI)};
    area_label = strrep(area_label,'_','-');
    
    ts = ROI_mean_ts_att(idx_attention(iROI),:);
    
    x = ROI.max_zscore_att(idx_attention(iROI)).x;
    y = ROI.max_zscore_att(idx_attention(iROI)).y;
    z = ROI.max_zscore_att(idx_attention(iROI)).z;
    
    f = figure;

    h = area([1 17],[max(zscore(ts)) max(zscore(ts))]);
    hold on
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([35 51],[max(zscore(ts)) max(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([69 85],[max(zscore(ts)) max(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([103 119],[max(zscore(ts)) max(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([137 153],[max(zscore(ts)) max(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    
    h = area([1 17],[min(zscore(ts)) min(zscore(ts))]);
    hold on
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([35 51],[min(zscore(ts)) min(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([69 85],[min(zscore(ts)) min(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([103 119],[min(zscore(ts)) min(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([137 153],[min(zscore(ts)) min(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);

    plot(zscore(ts));
    hold on
    plot(zscore(B(1:nTRs)),'r');
    
    xlim([1 nTRs]);
    ylim([min(zscore(ts)) max(zscore(ts))]);
    xlabel('TRs');
    ylabel('BOLD (z-scored)');
    ax = gca;
    %set(ax,'XTickLabel',{'1' '17' '34' '51' '68' '85' '102' '119' '136' '153' '170'});
    set(ax,'XTick',[1 17 34 51 68 85 102 119 136 153 170]);
    %legend({'BOLD' 'HRF'});
    
    title(strcat(area_label,':x:',int2str(x),':y:',int2str(y),':z:',int2str(z)));
    
    print(f,'-depsc',strcat('LHR','-','All-Subjects','-',condition,'-',area_label,'.eps'));
    print(f,'-djpeg',strcat('LHR','-','All-Subjects','-',condition,'-',area_label,'.jpg'));
    
    %%% RESTING
    
    condition = 'RestingState';
    
    area_label = ROI.label{idx_resting(iROI)};
    area_label = strrep(area_label,'_','-');
    
    ts = ROI_mean_ts_rest(idx_resting(iROI),:);
    
    x = ROI.max_zscore_rest(idx_resting(iROI)).x;
    y = ROI.max_zscore_rest(idx_resting(iROI)).y;
    z = ROI.max_zscore_rest(idx_resting(iROI)).z;
    
    f = figure;
    
    h = area([1 17]+17,[max(zscore(ts)) max(zscore(ts))]);
    hold on
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([35 51]+17,[max(zscore(ts)) max(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([69 85]+17,[max(zscore(ts)) max(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([103 119]+17,[max(zscore(ts)) max(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([137 153]+17,[max(zscore(ts)) max(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    
    h = area([1 17]+17,[min(zscore(ts)) min(zscore(ts))]);
    hold on
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([35 51]+17,[min(zscore(ts)) min(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([69 85]+17,[min(zscore(ts)) min(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([103 119]+17,[min(zscore(ts)) min(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    h = area([137 153]+17,[min(zscore(ts)) min(zscore(ts))]);
    set(h,'FaceColor',[0.75 0.75 0.75]);
    
    plot(zscore(ts));
    hold on
    plot(zscore(B(1:nTRs)),'r');
    
    xlim([1 nTRs]);
    ylim([min(zscore(ts)) max(zscore(ts))]);
    xlabel('TRs');
    ylabel('BOLD (z-scored)');
    ax = gca;
    %set(ax,'XTickLabel',{'1' '17' '34' '51' '68' '85' '102' '119' '136' '153' '170'});
    set(ax,'XTick',[1 17 34 51 68 85 102 119 136 153 170]);
    %legend({'BOLD' 'HRF'});
    
    title(strcat(area_label,':x:',int2str(x),':y:',int2str(y),':z:',int2str(z)));
    
    print(f,'-depsc',strcat('LHR','-','All-Subjects','-',condition,'-',area_label,'.eps'));
    print(f,'-djpeg',strcat('LHR','-','All-Subjects','-',condition,'-',area_label,'.jpg'));
    
end


end



    
    