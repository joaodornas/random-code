%function real_FC_voxel_AAL_ROI_densities_ANOVA

%%% GET SUBJECTS SETTINGS

all_settings = getAllSettings;
nSubjects = length(all_settings);

%%% LOAD AAL INFO

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nAllVoxels = length(find(AAL_img));

nRuns = 4;

nROI = 90;

contrast_TP_all = zeros(size(AAL_img));
contrast_PR_all = zeros(size(AAL_img));
contrast_TP_neg = zeros(size(AAL_img));
contrast_PR_neg = zeros(size(AAL_img));
contrast_TP_pos = zeros(size(AAL_img));
contrast_PR_pos = zeros(size(AAL_img));

contrast_PT_all = zeros(size(AAL_img));
contrast_RP_all = zeros(size(AAL_img));
contrast_PT_neg = zeros(size(AAL_img));
contrast_RP_neg = zeros(size(AAL_img));
contrast_PT_pos = zeros(size(AAL_img));
contrast_RP_pos = zeros(size(AAL_img));

for iROI=1:nROI
   
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(strcat(int2str(iROI),':',label_ROI{iROI},':',int2str(length(idx_voxels)),':voxels',':',datestr(now)));
    
    nVoxels = length(idx_voxels);
    
    for iVoxel=1:nVoxels
       
        if mod(iVoxel,1000) == 0
            
            disp(strcat('iVoxel:',int2str(iVoxel),':',datestr(now)));
        
        end
        
        samples_track_all = zeros(1,nRuns*nSubjects);
        samples_passive_all = zeros(1,nRuns*nSubjects);
        samples_restingstate_all = zeros(1,nRuns*nSubjects);
        
        samples_track_pos = zeros(1,nRuns*nSubjects);
        samples_passive_pos = zeros(1,nRuns*nSubjects);
        samples_restingstate_pos = zeros(1,nRuns*nSubjects);
        
        samples_track_neg = zeros(1,nRuns*nSubjects);
        samples_passive_neg = zeros(1,nRuns*nSubjects);
        samples_restingstate_neg = zeros(1,nRuns*nSubjects);
        
        it = 0;
        ip = 0;
        ir = 0;
        
        for isubject=1:nSubjects
            
            settings = all_settings(isubject).settings;
            
            %%% TRACK
            load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-density-Track-',area_label,'.mat'));
            
            for irun=1:4
            
                it = it + 1;
                samples_track_all(it) = FC_Density.run(irun).voxels_all_density(iVoxel);
                samples_track_pos(it) = FC_Density.run(irun).voxels_pos_density(iVoxel);
                samples_track_neg(it) = FC_Density.run(irun).voxels_neg_density(iVoxel);
            
            end
            
            clear FC_Density
            
            %%% PASSIVE
            load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-density-Passive-',area_label,'.mat'));
            
            for irun=1:4
            
                ip = ip + 1;
                samples_passive_all(ip) = FC_Density.run(irun).voxels_all_density(iVoxel);
                samples_passive_pos(ip) = FC_Density.run(irun).voxels_pos_density(iVoxel);
                samples_passive_neg(ip) = FC_Density.run(irun).voxels_neg_density(iVoxel);
            
            end
            
            clear FC_Density
            
            %%% RESTING STATE
            load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr-density-RestingState-',area_label,'.mat'));
            
            for irun=1:4
            
                ir = ir + 1;
                samples_restingstate_all(ir) = FC_Density.run(irun).voxels_all_density(iVoxel);
                samples_restingstate_pos(ir) = FC_Density.run(irun).voxels_pos_density(iVoxel);
                samples_restingstate_neg(ir) = FC_Density.run(irun).voxels_neg_density(iVoxel);
            
            end
            
            clear FC_Density
            
        end
        
        [H_TP_all,P_TP_all] = ttest(samples_track_all,samples_passive_all);
        [H_PR_all,P_PR_all] = ttest(samples_passive_all,samples_restingstate_all);
        
        [H_TP_pos,P_TP_pos] = ttest(samples_track_pos,samples_passive_pos);
        [H_PR_pos,P_PR_pos] = ttest(samples_passive_pos,samples_restingstate_pos);
        
        [H_TP_neg,P_TP_neg] = ttest(samples_track_neg,samples_passive_neg);
        [H_PR_neg,P_PR_neg] = ttest(samples_passive_neg,samples_restingstate_neg);
        
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
        
        if isnan(H_TP_all), H_TP_all = 0; end
        if isnan(H_PR_all), H_PR_all = 0; end
        if isnan(H_TP_pos), H_TP_pos = 0; end
        if isnan(H_PR_pos), H_PR_pos = 0; end
        if isnan(H_TP_neg), H_TP_neg = 0; end
        if isnan(H_PR_neg), H_PR_neg = 0; end
        
        if H_TP_all && (nanmean(samples_track_all)>nanmean(samples_passive_all))
            
            contrast_TP_all(idxx,idxy,idxz) = P_TP_all;
            
        elseif H_TP_all && (nanmean(samples_track_all)<nanmean(samples_passive_all))

            contrast_PT_all(idxx,idxy,idxz) = P_TP_all;
            
        end
        
        if H_PR_all && (nanmean(samples_passive_all)>nanmean(samples_restingstate_all))
            
            contrast_PR_all(idxx,idxy,idxz) = P_PR_all;
            
        elseif H_PR_all && (nanmean(samples_passive_all)<nanmean(samples_restingstate_all))

            contrast_RP_all(idxx,idxy,idxz) = P_PR_all;
            
        end
        
        if H_TP_pos && (nanmean(samples_track_pos)>nanmean(samples_passive_pos))
            
            contrast_TP_pos(idxx,idxy,idxz) = P_TP_pos;
            
        elseif H_TP_pos && (nanmean(samples_track_pos)<nanmean(samples_passive_pos))

            contrast_PT_pos(idxx,idxy,idxz) = P_TP_pos;
            
        end
            
        if H_PR_pos && (nanmean(samples_passive_pos)>nanmean(samples_restingstate_pos))
            
            contrast_PR_pos(idxx,idxy,idxz) = P_PR_pos;
            
        elseif H_PR_pos && (nanmean(samples_passive_pos)<nanmean(samples_restingstate_pos))

            contrast_RP_pos(idxx,idxy,idxz) = P_PR_pos;
            
        end
        
        if H_TP_neg && (nanmean(samples_track_neg)>nanmean(samples_passive_neg))
            
            contrast_TP_neg(idxx,idxy,idxz) = P_TP_neg;
            
        elseif H_TP_neg && (nanmean(samples_track_neg)<nanmean(samples_passive_neg))

            contrast_PT_neg(idxx,idxy,idxz) = P_TP_neg;
            
        end
            
        if H_PR_neg && (nanmean(samples_passive_neg)>nanmean(samples_restingstate_neg))
            
            contrast_PR_neg(idxx,idxy,idxz) = P_PR_neg;
            
        elseif H_PR_neg && (nanmean(samples_passive_neg)<nanmean(samples_restingstate_neg))

            contrast_RP_neg(idxx,idxy,idxz) = P_PR_neg;
            
        end
    
    end
    
end

nifti_file = load_aal;
offset = load_aal.dat.offset;
scl_slope = load_aal.dat.scl_slope;
scl_inter = load_aal.dat.scl_inter;
dtype = 'FLOAT32';
offset = 0;
dim = load_aal.dat.dim;

descrip = 'Contrast-Correlations';
experiment_label = 'LHR';

%%% CONTRAST

fname = strcat(experiment_label,'-','Correlation-Contrast-TP-All','.nii');
input_data = contrast_TP_all;
real_save_image;

fname = strcat(experiment_label,'-','Correlation-Contrast-PR-All','.nii');
input_data = contrast_PR_all;
real_save_image;

fname = strcat(experiment_label,'-','Correlation-Contrast-TP-Neg','.nii');
input_data = contrast_TP_neg;
real_save_image;

fname = strcat(experiment_label,'-','Correlation-Contrast-PR-Neg','.nii');
input_data = contrast_PR_neg;
real_save_image;

fname = strcat(experiment_label,'-','Correlation-Contrast-TP-Pos','.nii');
input_data = contrast_TP_pos;
real_save_image;

fname = strcat(experiment_label,'-','Correlation-Contrast-PR-Pos','.nii');
input_data = contrast_PR_pos;
real_save_image;

%%% INVERTED CONTRAST

fname = strcat(experiment_label,'-','Correlation-Contrast-PT-All','.nii');
input_data = contrast_PT_all;
real_save_image;

fname = strcat(experiment_label,'-','Correlation-Contrast-RP-All','.nii');
input_data = contrast_RP_all;
real_save_image;

fname = strcat(experiment_label,'-','Correlation-Contrast-PT-Neg','.nii');
input_data = contrast_PT_neg;
real_save_image;

fname = strcat(experiment_label,'-','Correlation-Contrast-RP-Neg','.nii');
input_data = contrast_RP_neg;
real_save_image;

fname = strcat(experiment_label,'-','Correlation-Contrast-PT-Pos','.nii');
input_data = contrast_PT_pos;
real_save_image;

fname = strcat(experiment_label,'-','Correlation-Contrast-RP-Pos','.nii');
input_data = contrast_RP_pos;
real_save_image;


%end

