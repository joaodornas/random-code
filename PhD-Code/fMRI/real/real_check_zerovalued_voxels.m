
%settings_subj1_2210;
%settings_subj2_2610;
%settings_subj3_0311;
%settings_subj4_0211;
%settings_subj5_0211;
settings_subj6_2411;

%%% LOAD DATA

subject = settings.folders.subject;

%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

real_load_all_data_FSL;

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;
nRun = length(Track);

AllZerosVoxels{1,2} = 'Track-1';
AllZerosVoxels{1,3} = 'Track-2';
AllZerosVoxels{1,4} = 'Track-3';
AllZerosVoxels{1,5} = 'Track-4';
AllZerosVoxels{1,6} = 'Passive-1';
AllZerosVoxels{1,7} = 'Passive-2';
AllZerosVoxels{1,8} = 'Passive-3';
AllZerosVoxels{1,9} = 'Passive-4';
AllZerosVoxels{1,10} = 'RestingState-1';
AllZerosVoxels{1,11} = 'RestingState-2';
AllZerosVoxels{1,12} = 'RestingState-3';
AllZerosVoxels{1,13} = 'RestingState-4';
AllZerosVoxels{1,14} = 'Trials-1';
AllZerosVoxels{1,15} = 'Trials-2';
AllZerosVoxels{1,16} = 'Trials-3';
AllZerosVoxels{1,17} = 'Trials-4';

for iROI=1:nNodes
    
    label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
    idx_ROI = AAL_ROI(iROI).ID;
    
    idx_voxels = find(AAL_img==idx_ROI);
    nVoxels = length(idx_voxels);
    
    area_label = strrep(label_ROI{iROI},'_','-');
    
    disp(strcat(int2str(iROI),':',label_ROI{iROI},':',int2str(length(idx_voxels)),':voxels'));
    
    %%% Track
    for iRun=1:nRun
        
        nZeros = 0;
        
        for iVoxel=1:nVoxels
    
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
            
            voxel_ts = squeeze(Track(iRun).run(idxx,idxy,idxz,:));
            
            if sum(voxel_ts) == 0;
                
                nZeros = nZeros + 1;
                
            end
            
            TrackZeroVoxels(iRun) = nZeros;
        
        end
        
    end
        
    %%% Passive
    for iRun=1:nRun
        
        nZeros = 0;
        
        for iVoxel=1:nVoxels
    
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
            
            voxel_ts = squeeze(Passive(iRun).run(idxx,idxy,idxz,:));
            
            if sum(voxel_ts) == 0;
                
                nZeros = nZeros + 1;
                
            end
            
            PassiveZeroVoxels(iRun) = nZeros;
        
        end
        
    end
    
    %%% RestingState
    for iRun=1:nRun
        
        nZeros = 0;
        
        for iVoxel=1:nVoxels
    
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
            
            voxel_ts = squeeze(RestingState(iRun).run(idxx,idxy,idxz,:));
            
            if sum(voxel_ts) == 0;
                
                nZeros = nZeros + 1;
                
            end
            
            RestingStateZeroVoxels(iRun) = nZeros;
        
        end
        
    end
    
    %%% Trials
    for iRun=1:nRun
        
        nZeros = 0;
        
        for iVoxel=1:nVoxels
    
            [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
            
            voxel_ts = squeeze(Trials(iRun).run(idxx,idxy,idxz,:));
            
            if sum(voxel_ts) == 0;
                
                nZeros = nZeros + 1;
                
            end
            
            TrialsZeroVoxels(iRun) = nZeros;
        
        end
        
    end

    AllZerosVoxels{iROI+1,1} = area_label;
    AllZerosVoxels{iROI+1,2} = TrackZeroVoxels(1);
    AllZerosVoxels{iROI+1,3} = TrackZeroVoxels(2);
    AllZerosVoxels{iROI+1,4} = TrackZeroVoxels(3);
    AllZerosVoxels{iROI+1,5} = TrackZeroVoxels(4);
    AllZerosVoxels{iROI+1,6} = PassiveZeroVoxels(1);
    AllZerosVoxels{iROI+1,7} = PassiveZeroVoxels(2);
    AllZerosVoxels{iROI+1,8} = PassiveZeroVoxels(3);
    AllZerosVoxels{iROI+1,9} = PassiveZeroVoxels(4);
    AllZerosVoxels{iROI+1,10} = RestingStateZeroVoxels(1);
    AllZerosVoxels{iROI+1,11} = RestingStateZeroVoxels(2);
    AllZerosVoxels{iROI+1,12} = RestingStateZeroVoxels(3);
    AllZerosVoxels{iROI+1,13} = RestingStateZeroVoxels(4);
    AllZerosVoxels{iROI+1,14} = TrialsZeroVoxels(1);
    AllZerosVoxels{iROI+1,15} = TrialsZeroVoxels(2);
    AllZerosVoxels{iROI+1,16} = TrialsZeroVoxels(3);
    AllZerosVoxels{iROI+1,17} = TrialsZeroVoxels(4);
        
end