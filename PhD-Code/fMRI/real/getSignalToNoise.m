
function getSignalToNoise

settings_subj1_2210;
doTheMath(settings);

settings_subj2_2610;
doTheMath(settings);

settings_subj3_0311;
doTheMath(settings);

settings_subj4_0211;
doTheMath(settings);
 
settings_subj5_0211;
doTheMath(settings);
 
settings_subj6_2411;
doTheMath(settings);

settings_subj7_1401;
doTheMath(settings);

settings_subj8_1401;
doTheMath(settings);


end


function doTheMath(settings)

subject = settings.folders.subject;

disp(strcat(subject,'-',datestr(now)));

%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

real_load_all_data_FSL;

nRuns = 4;
nConditions = 3;
nTotalVoxels = 160990;

idx_resting = 1;
idx_passive = 2;
idx_attention = 3;

ms_ts = zeros(nConditions,nRuns,nTotalVoxels);
std_ts = zeros(nConditions,nRuns,nTotalVoxels);

for iRun=1:nRuns

    [m_ts(idx_resting,iRun,:), std_ts(idx_resting,iRun,:)] = getSNR(RestingState(iRun).run);
    [m_ts(idx_passive,iRun,:), std_ts(idx_passive,iRun,:)] = getSNR(Passive(iRun).run);
    [m_ts(idx_attention,iRun,:), std_ts(idx_attention,iRun,:)] = getSNR(Track(iRun).run);
    
end

save(strcat('Mean-STD-Voxels','-',settings.codes.subject,'.mat'),'m_ts','std_ts');

end

function [m_ts, std_ts] = getSNR(run)

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\_ATLAS\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

MNI_size = [91 109 91];

AAL_ROI = load_roi.ROI;

nRuns = 4;
nROIs = 90;
iiVoxel = 0;
nTotalVoxels = 160990;

m_ts = zeros(1,nTotalVoxels);
std_ts = zeros(1,nTotalVoxels);

for iROI=1:nROIs
    
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    nVoxels = length(idx_voxels);
    
    disp(strcat(int2str(iROI),':',label_ROI{iROI},':',int2str(nVoxels),':voxels'));
    
    for iVoxel=1:nVoxels
        
       iiVoxel = iiVoxel + 1; 
       
       [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
       
       ts = squeeze(run(idxx,idxy,idxz,:));
       
       m_ts(iiVoxel) = mean(ts);
       std_ts(iiVoxel) = std(ts,0,1);
        
    end
    
end

end
