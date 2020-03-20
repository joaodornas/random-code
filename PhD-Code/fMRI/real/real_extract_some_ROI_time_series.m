function real_extract_some_ROI_time_series

% idx_ROI = [4 28 68 61 86]; %% 4 - Frontal Sup R; %% 28 - Rectus R; %% 68 - Precuneus - R; %% 61 - Parietal Inf L; %% 86 - Temporal Mid R

% settings_subj1_2210;
% doTheMath(settings,idx_ROI);
% 
% settings_subj2_2610;
% doTheMath(settings,idx_ROI);
% 
% settings_subj3_0311;
% doTheMath(settings,idx_ROI);
% 
% settings_subj4_0211;
% doTheMath(settings,idx_ROI);
% 
% settings_subj5_0211;
% doTheMath(settings,idx_ROI);
% 
% settings_subj6_2411;
% doTheMath(settings,idx_ROI);
% 
% settings_subj7_1401;
% doTheMath(settings,idx_ROI);
% 
% settings_subj8_1401;
% doTheMath(settings,idx_ROI);




settings_subj8_1401;
doTheMathMean(settings);

end

function doTheMath(settings,idx_ROI)

%%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.filtered;
mask = settings.FSL.files.mask.custom;

real_load_all_data_FSL;

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRuns = 4;
nTR = size(RestingState(1).run,4);

subject = settings.folders.subject;

for iROI=1:length(idx_ROI)
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    idx_region = AAL_ROI(idx_ROI(iROI)).ID;
    idx_voxels = find(AAL_img == idx_region);
    nVoxels = length(idx_voxels);
    
    disp(strcat('nVoxels:',int2str(nVoxels)));
    
    for iRun=1:nRuns
        
        ROI_RestingState.run(iRun).voxels = zeros(nVoxels,nTR);
        ROI_Passive.run(iRun).voxels = zeros(nVoxels,nTR);
        ROI_Track.run(iRun).voxels = zeros(nVoxels,nTR);
        
    end
    
    for iVoxel=1:nVoxels
        
       [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
       
       for iRun=1:nRuns
           
            ROI_RestingState.run(iRun).voxels(iVoxel,:) = squeeze(RestingState(iRun).run(idxx,idxy,idxz,:));
            ROI_Passive.run(iRun).voxels(iVoxel,:) = squeeze(Passive(iRun).run(idxx,idxy,idxz,:));
            ROI_Track.run(iRun).voxels(iVoxel,:) = squeeze(Track(iRun).run(idxx,idxy,idxz,:));
            
       end
       
    end
    
    save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',area_label,'.mat'),'subject','area_label','nVoxels','ROI_Track','ROI_Passive','ROI_RestingState');
    
    clear ROI_RestingState
    clear ROI_Passive
    clear ROI_Track
    
end

end

function doTheMathMean(settings)

%%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.filtered;
mask = settings.FSL.files.mask.custom;

real_load_all_data_FSL;

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRuns = 4;
nROI = 90;
idx_ROI = 1:nROI;
nTR = size(RestingState(1).run,4);

subject = settings.folders.subject;

for iROI=1:nROI
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    idx_region = AAL_ROI(idx_ROI(iROI)).ID;
    idx_voxels = find(AAL_img == idx_region);
    nVoxels = length(idx_voxels);
    
    disp(strcat('nVoxels:',int2str(nVoxels)));
    
    for iRun=1:nRuns
        
        ROI_RestingState.run(iRun).voxels = zeros(nVoxels,nTR);
        ROI_Passive.run(iRun).voxels = zeros(nVoxels,nTR);
        ROI_Track.run(iRun).voxels = zeros(nVoxels,nTR);
        
    end
    
    for iVoxel=1:nVoxels
        
       [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
       
       for iRun=1:nRuns
           
            ROI_RestingState.run(iRun).voxels(iVoxel,:) = squeeze(RestingState(iRun).run(idxx,idxy,idxz,:));
            ROI_Passive.run(iRun).voxels(iVoxel,:) = squeeze(Passive(iRun).run(idxx,idxy,idxz,:));
            ROI_Track.run(iRun).voxels(iVoxel,:) = squeeze(Track(iRun).run(idxx,idxy,idxz,:));
            
       end
       
    end
    
    for iRun=1:nRuns
        
       ROI(iROI).Mean_RestingState.run(iRun,:) = mean(ROI_RestingState.run(iRun).voxels,1);
       ROI(iROI).Mean_Passive.run(iRun,:) = mean(ROI_Passive.run(iRun).voxels,1);
       ROI(iROI).Mean_Track.run(iRun,:) = mean(ROI_Track.run(iRun).voxels,1);
        
    end

    clear ROI_RestingState
    clear ROI_Passive
    clear ROI_Track
      
end

save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','Mean-AAL-Time-Series','.mat'),'subject','ROI');    

end

