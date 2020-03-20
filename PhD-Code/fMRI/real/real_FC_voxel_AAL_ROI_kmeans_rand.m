function real_FC_voxel_AAL_ROI_kmeans_rand

%settings_subj1_2210;
%settings_subj2_2610;
%settings_subj3_0311;
%settings_subj4_0211;
%settings_subj5_0211;
%settings_subj6_2411;

% settings_subj1_2210;
% doTheMath(settings);
% clear settings
% 
% settings_subj2_2610;
% doTheMath(settings);
% clear settings
% 
% settings_subj3_0311;
% doTheMath(settings);
% clear settings
% 
% settings_subj4_0211;
% doTheMath(settings);
% clear settings
% 
% settings_subj5_0211;
% doTheMath(settings);
% clear settings
% 
% settings_subj6_2411;
% doTheMath(settings);
% clear settings

doTheMathAllSubjects;

end


function doTheMath(settings)

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_ROI = 1:90;
nRuns = 4;
nROI = length(idx_ROI);

for iROI=1:nROI
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    TrackRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-KMeans','.mat'));
    PassiveRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-KMeans','.mat'));
    RestingStateRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-KMeans','.mat'));
        
end

end

function doTheMathAllSubjects

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

idx_ROI = 1:90;
nROI = length(idx_ROI);

RandIndex.nVoxels = zeros(1,nROI);
RandIndex.n11_TrackPassive = zeros(1,nROI);
RandIndex.n00_TrackPassive = zeros(1,nROI);
RandIndex.n11_TrackRestingState = zeros(1,nROI);
RandIndex.n00_TrackRestingState = zeros(1,nROI);
RandIndex.n11_PassiveRestingState = zeros(1,nROI);
RandIndex.n00_PassiveRestingState = zeros(1,nROI);

for iROI=1:nROI
    
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-KMeans','.mat'));
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-KMeans','.mat'));
    load(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-KMeans','.mat'));
        
    RandIndex.nVoxels(iROI) = length(FC_Track_KMeans.IdxClusters);
    
    for iVoxel=1:RandIndex.nVoxels(iROI)
        
       for iiVoxel=iVoxel:RandIndex.nVoxels(iROI)
           
          if (FC_Track_KMeans.IdxClusters(iVoxel) == FC_Track_KMeans.IdxClusters(iiVoxel)) && (FC_Passive_KMeans.IdxClusters(iVoxel) == FC_Passive_KMeans.IdxClusters(iiVoxel)) 
              
              RandIndex.n11_TrackPassive(iROI) = RandIndex.n11_TrackPassive(iROI) + 1;
              
          elseif (FC_Track_KMeans.IdxClusters(iVoxel) ~= FC_Track_KMeans.IdxClusters(iiVoxel)) && (FC_Passive_KMeans.IdxClusters(iVoxel) ~= FC_Passive_KMeans.IdxClusters(iiVoxel))
              
              RandIndex.n00_TrackPassive(iROI) = RandIndex.n00_TrackPassive(iROI) + 1;
              
          end
          
          if (FC_Passive_KMeans.IdxClusters(iVoxel) == FC_Passive_KMeans.IdxClusters(iiVoxel)) && (FC_RestingState_KMeans.IdxClusters(iVoxel) == FC_RestingState_KMeans.IdxClusters(iiVoxel)) 
              
              RandIndex.n11_PassiveRestingState(iROI) = RandIndex.n11_PassiveRestingState(iROI) + 1;
              
          elseif (FC_Passive_KMeans.IdxClusters(iVoxel) ~= FC_Passive_KMeans.IdxClusters(iiVoxel)) && (FC_RestingState_KMeans.IdxClusters(iVoxel) ~= FC_RestingState_KMeans.IdxClusters(iiVoxel))
              
              RandIndex.n00_PassiveRestingState(iROI) = RandIndex.n00_PassiveRestingState(iROI) + 1;
              
          end
          
          if (FC_Track_KMeans.IdxClusters(iVoxel) == FC_Track_KMeans.IdxClusters(iiVoxel)) && (FC_RestingState_KMeans.IdxClusters(iVoxel) == FC_RestingState_KMeans.IdxClusters(iiVoxel)) 
              
              RandIndex.n11_TrackRestingState(iROI) = RandIndex.n11_TrackRestingState(iROI) + 1;
              
          elseif (FC_Track_KMeans.IdxClusters(iVoxel) ~= FC_Track_KMeans.IdxClusters(iiVoxel)) && (FC_RestingState_KMeans.IdxClusters(iVoxel) ~= FC_RestingState_KMeans.IdxClusters(iiVoxel))
              
              RandIndex.n00_TrackRestingState(iROI) = RandIndex.n00_TrackRestingState(iROI) + 1;
              
          end
           
       end
        
        
    end
    
end

save(strcat('All-Subjects','-','FC-Voxels-AAL-ROI','-','Rand-Index','.mat'),'RandIndex');

end


