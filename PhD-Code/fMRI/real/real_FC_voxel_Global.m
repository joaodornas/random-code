
function real_FC_voxel_Global

% extractDataFromCluster;

% identifyMissingVoxels;

% doContrasts;

% doContrastsGlobalNormalization;

doContrastsGlobalNormalizationPerRun;


end

function extractDataFromCluster

nProcesses = 1872;
nRunsPerCondition = 32;
nRuns = 96;
condition_label = {'Track' 'Passive' 'Rest'};

folder = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\GLOBAL\densities\densities-outside';

for iCondition=1:length(condition_label)
    
   disp(condition_label{iCondition});
   
   for iRun=((iCondition-1)*nRunsPerCondition+1):(iCondition)*nRunsPerCondition
       
       disp(strcat('run:',int2str(iRun)));
 
       density_label = 'pos';
       [error_pos dataFromFiles_pos] = extractDataFromFile(folder,iRun,density_label,nProcesses);
       
       density_label = 'neg';
       [error_neg dataFromFiles_neg] = extractDataFromFile(folder,iRun,density_label,nProcesses);
       
       all_runs_per_condition_pos = dataFromFiles_pos;
       all_runs_per_condition_neg = dataFromFiles_neg;
       
       idx_globals_pos = squeeze(all_runs_per_condition_pos(:,1));
   
       [s_idx_globals_pos, i_idx_globals_pos] = sort(idx_globals_pos,'ascend');

       eval(strcat(condition_label{iCondition},'.run_pos = all_runs_per_condition_pos(i_idx_globals_pos,1:4)'));

       idx_globals_neg = squeeze(all_runs_per_condition_neg(:,1));

       [s_idx_globals_neg, i_idx_globals_neg] = sort(idx_globals_neg,'ascend');

       eval(strcat(condition_label{iCondition},'.run_neg = all_runs_per_condition_neg(i_idx_globals_neg,1:4)'));
       
       eval(strcat(condition_label{iCondition},'.error_neg = error_neg'));
       eval(strcat(condition_label{iCondition},'.error_pos = error_pos'));

       save(strcat('Global-',condition_label{iCondition},'-Run-',int2str(iRun),'.mat'),condition_label{iCondition});
       
       eval(sprintf('clear %s',condition_label{iCondition}));
          
   end
      
end


end

function [errorVoxels dataFromFile] = extractDataFromFile(folder,iRun,density_label,nProcesses)

iVoxel = 0;
iError = 0;
errorVoxels = cell.empty;
for iProcess=0:(nProcesses-1)
    
    disp(strcat('process:',int2str(iProcess)));
    
    filename = strcat(folder,'\','run-',int2str(iRun),'-',int2str(iProcess),'-','density','-',density_label,'.txt');
    fid  = fopen(filename,'r');

    tline = fgets(fid);

    while ischar(tline)

        iVoxel = iVoxel + 1;
        line = tline(:);

        idx_comma = strfind(tline,',');
        
        if length(idx_comma) <= 3

            try
                
                dataFromFile(iVoxel,1) = str2num(tline(1:idx_comma(1)-1));
                
            catch
                
                dataFromFile(iVoxel,1) = NaN;
                
            end
            
            try
                
                dataFromFile(iVoxel,2) = str2num(tline(idx_comma(1)+1:idx_comma(2)-1));
        
            catch
                
                dataFromFIle(iVoxel,2) = NaN;
                
            end
            
            try
                
                dataFromFile(iVoxel,3) = str2num(tline(idx_comma(2)+1:idx_comma(3)-1));
        
            catch
                
                dataFromFile(iVoxel,3) = NaN;
                
            end
               
            try
               
                dataFromFile(iVoxel,4) = str2num(tline(idx_comma(3)-1:end));
        
            catch
                
                dataFromFile(iVoxel,4) = NaN;
                
            end
            
        else
            
            iError = iError + 1;
            errorVoxels{iError} = tline(:);
            
        end

        tline = fgets(fid);

    end

    fclose(fid);

end


end

function identifyMissingVoxels

nRuns = 96;
nRunsPerCondition = 32;
condition_label = {'Track' 'Passive' 'Rest'};

pos = [];
neg = [];
for iCondition=1:length(condition_label)
    
   disp(condition_label{iCondition});
   
   for iRun=((iCondition-1)*nRunsPerCondition+1):(iCondition)*nRunsPerCondition
       
       disp(strcat('run:',int2str(iRun)));
       
       eval(strcat('load(''Global-',condition_label{iCondition},'-Run-',int2str(iRun),'.mat'');'));
       
       eval(strcat('pos=[pos;',condition_label{iCondition},'.run_pos];'));
       eval(strcat('neg=[neg;',condition_label{iCondition},'.run_neg];'));
       
   end
   
end

idx_pos = squeeze(pos(:,1));
idx_neg = squeeze(neg(:,1));

idx_pos(idx_pos==0) = [];
idx_neg(idx_neg==0) = [];

idx_pos = sort(idx_pos);
idx_neg = sort(idx_neg);

%isequal(idx_pos,idx_neg)

totalVoxels = 160990*96;
idx_all = 1:totalVoxels;

members = ismember(idx_all,idx_pos);

idx_missing = find(members==0);

fid = fopen('missing_voxels.txt','w');

for iline=1:length(idx_missing)
   
    fprintf(fid,'%s\n',int2str(idx_missing(iline)));
    
end

fclose(fid);

end

function doContrasts

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nROI = 90;
nRuns = 32;
imissing = 0;
idx_ROI = 1:nROI;

disp('Loading Data');

for iRun=1:nRuns
     load(strcat('Global-Track-Run-',int2str(iRun),'.mat'));
     run(iRun).condition = Track;
     clear Track
end

for iRun=(nRuns+1):nRuns*2
     load(strcat('Global-Passive-Run-',int2str(iRun),'.mat'));
     run(iRun).condition = Passive;
     clear Passive
end

for iRun=(nRuns*2+1):nRuns*3
     load(strcat('Global-Rest-Run-',int2str(iRun),'.mat'));
     run(iRun).condition = Rest;
     clear Rest
end

attention_contrast_pos_decrease = zeros(size(AAL_img));
attention_contrast_neg_decrease = zeros(size(AAL_img));

attention_contrast_pos_increase = zeros(size(AAL_img));
attention_contrast_neg_increase = zeros(size(AAL_img));

attention_contrast_pos_equal = zeros(size(AAL_img));
attention_contrast_neg_equal = zeros(size(AAL_img));

stimulus_contrast_pos_decrease = zeros(size(AAL_img));
stimulus_contrast_neg_decrease = zeros(size(AAL_img));

stimulus_contrast_pos_increase = zeros(size(AAL_img));
stimulus_contrast_neg_increase = zeros(size(AAL_img));

stimulus_contrast_pos_equal = zeros(size(AAL_img));
stimulus_contrast_neg_equal = zeros(size(AAL_img));

taskrest_contrast_pos_decrease = zeros(size(AAL_img));
taskrest_contrast_neg_decrease = zeros(size(AAL_img));

taskrest_contrast_pos_increase = zeros(size(AAL_img));
taskrest_contrast_neg_increase = zeros(size(AAL_img));

taskrest_contrast_pos_equal = zeros(size(AAL_img));
taskrest_contrast_neg_equal = zeros(size(AAL_img));

disp('...looping over voxels');

iiVoxel = 0;

allVoxelsDensities.track_pos = zeros(1,nRuns+1);
allVoxelsDensities.track_neg = zeros(1,nRuns+1);
allVoxelsDensities.passive_pos = zeros(1,nRuns+1);
allVoxelsDensities.passive_neg = zeros(1,nRuns+1);
allVoxelsDensities.rest_pos = zeros(1,nRuns+1);
allVoxelsDensities.rest_neg = zeros(1,nRuns+1);

for iROI=1:nROI
   
    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    idx_voxels = find(AAL_img==AAL_ROI(idx_ROI(iROI)).ID);
    
    nVoxels = length(idx_voxels);
    
    disp(strcat('nVoxes=',int2str(nVoxels)));
    
    for iVoxel=1:nVoxels
        
        iiVoxel = iiVoxel + 1;
        
        if mod(iVoxel,1000) == 0; disp(strcat('...voxels so far:',int2str(iVoxel))); end
       
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
        
        track_pos = zeros(1,nRuns);
        track_neg = zeros(1,nRuns);
        
        passive_pos = zeros(1,nRuns);
        passive_neg = zeros(1,nRuns);
        
        rest_pos = zeros(1,nRuns);
        rest_neg = zeros(1,nRuns);
        
        iiRun = 0;
        for iRun=1:nRuns
            
            iiRun = iiRun + 1;
            
            %load(strcat('Global-Track-Run-',int2str(iRun),'.mat'));
            Track = run(iRun).condition;
            
            idx_global = squeeze(Track.run_pos(:,1));
            idx_roi = squeeze(Track.run_pos(:,2));
            idx_voxel = squeeze(Track.run_pos(:,3));
            
            %idx_voxel_on_run = find(idx_global==iiAll_Voxels);
            idx_roi_on_run = find(idx_roi==iROI);
            idx_local_voxel_on_run = find(idx_voxel==iVoxel);
            
            idx_voxel_on_run = idx_roi_on_run(ismember(idx_roi_on_run,idx_local_voxel_on_run));
            
            if ~isempty(idx_voxel_on_run)
                
                track_pos(iiRun) = squeeze(Track.run_pos(idx_voxel_on_run,4));
                track_neg(iiRun) = squeeze(Track.run_neg(idx_voxel_on_run,4));
                
                imissing = imissing + 1;
                not_missing(imissing) = squeeze(Track.run_pos(idx_voxel_on_run,1));
            
            else
                
                track_pos(iiRun) = NaN;
                track_neg(iiRun) = NaN;
                
            end
            
            clear Track
            
        end
        
        iiRun = 0;
        for iRun=(nRuns+1):nRuns*2
            
            iiRun = iiRun + 1;
            
            %load(strcat('Global-Passive-Run-',int2str(iRun),'.mat'));
            Passive = run(iRun).condition;
            
            idx_global = squeeze(Passive.run_pos(:,1));
            idx_roi = squeeze(Passive.run_pos(:,2));
            idx_voxel = squeeze(Passive.run_pos(:,3));
            
            %idx_voxel_on_run = find(idx_global==iiAll_Voxels);
            idx_roi_on_run = find(idx_roi==iROI);
            idx_local_voxel_on_run = find(idx_voxel==iVoxel);
            
            idx_voxel_on_run = idx_roi_on_run(ismember(idx_roi_on_run,idx_local_voxel_on_run));
            
            if ~isempty(idx_voxel_on_run)
                
                passive_pos(iiRun) = squeeze(Passive.run_pos(idx_voxel_on_run,4));
                passive_neg(iiRun) = squeeze(Passive.run_neg(idx_voxel_on_run,4));
                
                imissing = imissing + 1;
                not_missing(imissing) = squeeze(Passive.run_pos(idx_voxel_on_run,1));
            
            else
                
                passive_pos(iiRun) = NaN;
                passive_neg(iiRun) = NaN;
                
            end
            
            clear Passive
            
        end
        
        iiRun = 0;
        for iRun=(nRuns*2+1):nRuns*3
            
            iiRun = iiRun + 1;
            
            %load(strcat('Global-Rest-Run-',int2str(iRun),'.mat'));
            Rest = run(iRun).condition;
            
            idx_global = squeeze(Rest.run_pos(:,1));
            idx_roi = squeeze(Rest.run_pos(:,2));
            idx_voxel = squeeze(Rest.run_pos(:,3));
            
            %idx_voxel_on_run = find(idx_global==iiAll_Voxels);
            idx_roi_on_run = find(idx_roi==iROI);
            idx_local_voxel_on_run = find(idx_voxel==iVoxel);
            
            idx_voxel_on_run = idx_roi_on_run(ismember(idx_roi_on_run,idx_local_voxel_on_run));
            
            if ~isempty(idx_voxel_on_run)
                
                rest_pos(iiRun) = squeeze(Rest.run_pos(idx_voxel_on_run,4));
                rest_neg(iiRun) = squeeze(Rest.run_neg(idx_voxel_on_run,4));
                
                imissing = imissing + 1;
                not_missing(imissing) = squeeze(Rest.run_pos(idx_voxel_on_run,1));
            
            else
                
                rest_pos(iiRun) = NaN;
                rest_neg(iiRun) = NaN;
                
            end
            
            clear Rest
            
        end
        
       [H,P] = ttest(track_pos,passive_pos);
       if nanmean(track_pos) > nanmean(passive_pos)
           attention_contrast_pos_increase(idxx,idxy,idxz) = P;
       elseif nanmean(track_pos) < nanmean(passive_pos)
           attention_contrast_pos_decrease(idxx,idxy,idxz) = P;
       elseif nanmean(track_pos) == nanmean(passive_pos)
           attention_contrast_pos_equal(idxx,idxy,idxz) = 1;
       end
       
       [H,P] = ttest(track_neg,passive_neg);
       if nanmean(track_neg) > nanmean(passive_neg)
           attention_contrast_neg_increase(idxx,idxy,idxz) = P;
       elseif nanmean(track_neg) < nanmean(passive_neg)
           attention_contrast_neg_decrease(idxx,idxy,idxz) = P;
       elseif nanmean(track_neg) == nanmean(passive_neg)
           attention_contrast_neg_equal(idxx,idxy,idxz) = 1;
       end
       
       [H,P] = ttest(passive_pos,rest_pos);
       if nanmean(passive_pos) > nanmean(rest_pos)
           stimulus_contrast_pos_increase(idxx,idxy,idxz) = P;
       elseif nanmean(passive_pos) < nanmean(rest_pos)
           stimulus_contrast_pos_decrease(idxx,idxy,idxz) = P;
       elseif nanmean(passive_pos) == nanmean(rest_pos)
           stimulus_contrast_pos_equal(idxx,idxy,idxz) = 1;
       end
       
       [H,P] = ttest(passive_neg,rest_neg);
       if nanmean(passive_neg) > nanmean(rest_neg)
           stimulus_contrast_neg_increase(idxx,idxy,idxz) = P;
       elseif nanmean(passive_neg) < nanmean(rest_neg)
           stimulus_contrast_neg_decrease(idxx,idxy,idxz) = P;
       elseif nanmean(passive_neg) == nanmean(rest_neg)
           stimulus_contrast_neg_equal(idxx,idxy,idxz) = 1;
       end
       
       [H,P] = ttest(track_pos,rest_pos);
       if nanmean(track_pos) > nanmean(rest_pos)
           taskrest_contrast_pos_increase(idxx,idxy,idxz) = P;
       elseif nanmean(track_pos) < nanmean(rest_pos)
           taskrest_contrast_pos_decrease(idxx,idxy,idxz) = P;
       elseif nanmean(track_pos) == nanmean(rest_pos)
           taskrest_contrast_pos_equal(idxx,idxy,idxz) = 1;
       end
       
       [H,P] = ttest(track_neg,rest_neg);
       if nanmean(track_neg) > nanmean(rest_neg)
           taskrest_contrast_neg_increase(idxx,idxy,idxz) = P;
       elseif nanmean(track_neg) < nanmean(rest_neg)
           taskrest_contrast_neg_decrease(idxx,idxy,idxz) = P;
       elseif nanmean(track_neg) == nanmean(rest_neg)
           taskrest_contrast_neg_equal(idxx,idxy,idxz) = 1;
       end
       
       allVoxelsDensities.track_pos(iiVoxel,1:end) = [track_pos(:)',nanmean(track_pos)];
       allVoxelsDensities.track_neg(iiVoxel,1:end) = [track_neg(:)',nanmean(track_neg)];
       allVoxelsDensities.passive_pos(iiVoxel,1:end) = [passive_pos(:)',nanmean(passive_pos)];
       allVoxelsDensities.passive_neg(iiVoxel,1:end) = [passive_neg(:)',nanmean(passive_neg)];
       allVoxelsDensities.rest_pos(iiVoxel,1:end) = [rest_pos(:)',nanmean(rest_pos)];
       allVoxelsDensities.rest_neg(iiVoxel,1:end) = [rest_neg(:)',nanmean(rest_neg)];
        
    end
    
end

save('not_missing.mat','not_missing');
save('allVoxelsDensities.mat','allVoxelsDensities');

scl_slope = 1;
scl_inter = 0; 
dim = [size(AAL_img,1),size(AAL_img,2),size(AAL_img,3)];
dtype = 'FLOAT32';
offset = 0;
nifti_file.mat = load_aal.mat;
nifti_file.mat_intent = load_aal.mat_intent;
nifti_file.mat0 = load_aal.mat0;
nifti_file.mat0_intent = load_aal.mat0_intent;

fname = strcat('Global-Attention-Pos-Increase','.nii');
descrip = 'Global-Attention-Pos-Increase';
input_data = attention_contrast_pos_increase; 
real_save_image;

fname = strcat('Global-Attention-Neg-Increase','.nii');
descrip = 'Global-Attention-Neg-Increase';
input_data = attention_contrast_neg_increase; 
real_save_image;

fname = strcat('Global-Stimululs-Pos-Increase','.nii');
descrip = 'Global-Stimulus-Pos-Increase';
input_data = stimulus_contrast_pos_increase; 
real_save_image;

fname = strcat('Global-Stimulus-Neg-Increase','.nii');
descrip = 'Global-Stimulus-Neg-Increase';
input_data = stimulus_contrast_neg_increase; 
real_save_image;

fname = strcat('Global-TaskRest-Pos-Increase','.nii');
descrip = 'Global-TaskRest-Pos-Increase';
input_data = taskrest_contrast_pos_increase; 
real_save_image;

fname = strcat('Global-TaskRest-Neg-Increase','.nii');
descrip = 'Global-TaskRest-Neg-Increase';
input_data = taskrest_contrast_neg_increase; 
real_save_image;

fname = strcat('Global-Attention-Pos-Decrease','.nii');
descrip = 'Global-Attention-Pos-Decrease';
input_data = attention_contrast_pos_decrease; 
real_save_image;

fname = strcat('Global-Attention-Neg-Decrease','.nii');
descrip = 'Global-Attention-Neg-Decrease';
input_data = attention_contrast_neg_decrease; 
real_save_image;

fname = strcat('Global-Stimululs-Pos-Decrease','.nii');
descrip = 'Global-Stimulus-Pos-Decrease';
input_data = stimulus_contrast_pos_decrease; 
real_save_image;

fname = strcat('Global-Stimulus-Neg-Decrease','.nii');
descrip = 'Global-Stimulus-Neg-Decrease';
input_data = stimulus_contrast_neg_decrease; 
real_save_image;

fname = strcat('Global-TaskRest-Pos-Decrease','.nii');
descrip = 'Global-TaskRest-Pos-Decrease';
input_data = taskrest_contrast_pos_decrease; 
real_save_image;

fname = strcat('Global-TaskRest-Neg-Decrease','.nii');
descrip = 'Global-TaskRest-Neg-Decrease';
input_data = taskrest_contrast_neg_decrease; 
real_save_image;

fname = strcat('Global-Attention-Pos-Equal','.nii');
descrip = 'Global-Attention-Pos-Equal';
input_data = attention_contrast_pos_equal; 
real_save_image;

fname = strcat('Global-Attention-Neg-Equal','.nii');
descrip = 'Global-Attention-Neg-Equal';
input_data = attention_contrast_neg_equal; 
real_save_image;

fname = strcat('Global-Stimululs-Pos-Equal','.nii');
descrip = 'Global-Stimulus-Pos-Equal';
input_data = stimulus_contrast_pos_equal; 
real_save_image;

fname = strcat('Global-Stimulus-Neg-Equal','.nii');
descrip = 'Global-Stimulus-Neg-Equal';
input_data = stimulus_contrast_neg_equal; 
real_save_image;

fname = strcat('Global-TaskRest-Pos-Equal','.nii');
descrip = 'Global-TaskRest-Pos-Equal';
input_data = taskrest_contrast_pos_equal; 
real_save_image;

fname = strcat('Global-TaskRest-Neg-Equal','.nii');
descrip = 'Global-TaskRest-Neg-Equal';
input_data = taskrest_contrast_neg_equal; 
real_save_image;

end

function doContrastsGlobalNormalization

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nROI = 90;
idx_ROI = 1:nROI;

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\GLOBAL\densities\outside-volumes\2nd\allVoxelsDensities.mat');

allVoxelsDensities.track_neg = allVoxelsDensities.track_neg(:,1:32);
allVoxelsDensities.track_pos = allVoxelsDensities.track_pos(:,1:32);

allVoxelsDensities.passive_neg = allVoxelsDensities.passive_neg(:,1:32);
allVoxelsDensities.passive_pos = allVoxelsDensities.passive_pos(:,1:32);

allVoxelsDensities.rest_neg = allVoxelsDensities.rest_neg(:,1:32);
allVoxelsDensities.rest_pos = allVoxelsDensities.rest_pos(:,1:32);

nSubjects = 8;
nRunsPerSubject = 4;

allVoxelsDensities.track_neg_norm = [];
allVoxelsDensities.track_pos_norm = [];
allVoxelsDensities.passive_neg_norm = [];
allVoxelsDensities.passive_pos_norm = [];
allVoxelsDensities.rest_neg_norm = [];
allVoxelsDensities.rest_pos_norm = [];

for iSubject=1:nSubjects
    
    all_track_neg = allVoxelsDensities.track_neg(:,((iSubject-1)*nRunsPerSubject+1):(iSubject*nRunsPerSubject));
    mean_all_track_neg = mean(all_track_neg(:));
    all_track_neg_norm = all_track_neg ./ mean_all_track_neg;
    
    allVoxelsDensities.track_neg_norm = [allVoxelsDensities.track_neg_norm,all_track_neg_norm];
    
    all_track_pos = allVoxelsDensities.track_pos(:,((iSubject-1)*nRunsPerSubject+1):(iSubject*nRunsPerSubject));
    mean_all_track_pos = mean(all_track_pos(:));
    all_track_pos_norm = all_track_pos ./ mean_all_track_pos;
    
    allVoxelsDensities.track_pos_norm = [allVoxelsDensities.track_pos_norm,all_track_pos_norm];
    
    all_passive_neg = allVoxelsDensities.passive_neg(:,((iSubject-1)*nRunsPerSubject+1):(iSubject*nRunsPerSubject));
    mean_all_passive_neg = mean(all_passive_neg(:));
    all_passive_neg_norm = all_passive_neg ./ mean_all_passive_neg;
    
    allVoxelsDensities.passive_neg_norm = [allVoxelsDensities.passive_neg_norm,all_passive_neg_norm];
    
    all_passive_pos = allVoxelsDensities.passive_pos(:,((iSubject-1)*nRunsPerSubject+1):(iSubject*nRunsPerSubject));
    mean_all_passive_pos = mean(all_passive_pos(:));
    all_passive_pos_norm = all_passive_pos ./ mean_all_passive_pos;
    
    allVoxelsDensities.passive_pos_norm = [allVoxelsDensities.passive_pos_norm,all_passive_pos_norm];
    
    all_rest_neg = allVoxelsDensities.rest_neg(:,((iSubject-1)*nRunsPerSubject+1):(iSubject*nRunsPerSubject));
    mean_all_rest_neg = mean(all_rest_neg(:));
    all_rest_neg_norm = all_rest_neg ./ mean_all_rest_neg;
    
    allVoxelsDensities.rest_neg_norm = [allVoxelsDensities.rest_neg_norm,all_rest_neg_norm];
    
    all_rest_pos = allVoxelsDensities.rest_pos(:,((iSubject-1)*nRunsPerSubject+1):(iSubject*nRunsPerSubject));
    mean_all_rest_pos = mean(all_rest_pos(:));
    all_rest_pos_norm = all_rest_pos ./ mean_all_rest_pos;
    
    allVoxelsDensities.rest_pos_norm = [allVoxelsDensities.rest_pos_norm,all_rest_pos_norm];
    
end

attention_contrast_pos_decrease = zeros(size(AAL_img));
attention_contrast_neg_decrease = zeros(size(AAL_img));

attention_contrast_pos_increase = zeros(size(AAL_img));
attention_contrast_neg_increase = zeros(size(AAL_img));

stimulus_contrast_pos_decrease = zeros(size(AAL_img));
stimulus_contrast_neg_decrease = zeros(size(AAL_img));

stimulus_contrast_pos_increase = zeros(size(AAL_img));
stimulus_contrast_neg_increase = zeros(size(AAL_img));

nTotalVoxels = 160990;

H_attention_neg = zeros(1,nTotalVoxels);
H_attention_pos = zeros(1,nTotalVoxels);

H_stimulus_neg = zeros(1,nTotalVoxels);
H_stimulus_pos = zeros(1,nTotalVoxels);

p_attention_neg = zeros(1,nTotalVoxels);
p_attention_pos = zeros(1,nTotalVoxels);

p_stimulus_neg = zeros(1,nTotalVoxels);
p_stimulus_pos = zeros(1,nTotalVoxels);

iiVoxel = 0;

for iROI=1:nROI

    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    idx_voxels = find(AAL_img==AAL_ROI(idx_ROI(iROI)).ID);
    
    nVoxels = length(idx_voxels);
    
    disp(strcat('nVoxes=',int2str(nVoxels)));
    
    for iVoxel=1:nVoxels
        
        iiVoxel = iiVoxel + 1;
        
        if mod(iVoxel,1000) == 0; disp(strcat('...voxels so far:',int2str(iVoxel))); end
       
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
    
        [H_attention_neg(iiVoxel),p_attention_neg(iiVoxel)] = ttest(squeeze(allVoxelsDensities.track_neg_norm(iiVoxel,:)),squeeze(allVoxelsDensities.passive_neg_norm(iiVoxel,:)));
        
        if nanmean(allVoxelsDensities.track_neg_norm(iiVoxel,:)) > nanmean(allVoxelsDensities.passive_neg_norm(iiVoxel,:))
           attention_contrast_neg_increase(idxx,idxy,idxz) = p_attention_neg(iiVoxel);
        elseif nanmean(allVoxelsDensities.track_neg_norm(iiVoxel,:)) < nanmean(allVoxelsDensities.passive_neg_norm(iiVoxel,:))
           attention_contrast_neg_decrease(idxx,idxy,idxz) = p_attention_neg(iiVoxel);
        end
        
        [H_attention_pos(iiVoxel),p_attention_pos(iiVoxel)] = ttest(squeeze(allVoxelsDensities.track_pos_norm(iiVoxel,:)),squeeze(allVoxelsDensities.passive_pos_norm(iiVoxel,:))); 
   
        if nanmean(allVoxelsDensities.track_pos_norm(iiVoxel,:)) > nanmean(allVoxelsDensities.passive_pos_norm(iiVoxel,:))
           attention_contrast_pos_increase(idxx,idxy,idxz) = p_attention_pos(iiVoxel);
        elseif nanmean(allVoxelsDensities.track_pos_norm(iiVoxel,:)) < nanmean(allVoxelsDensities.passive_pos_norm(iiVoxel,:))
           attention_contrast_pos_decrease(idxx,idxy,idxz) = p_attention_pos(iiVoxel);
        end
        
        [H_stimulus_neg(iiVoxel),p_stimulus_neg(iiVoxel)] = ttest(squeeze(allVoxelsDensities.passive_neg_norm(iiVoxel,:)),squeeze(allVoxelsDensities.rest_neg_norm(iiVoxel,:)));
        
        if nanmean(allVoxelsDensities.passive_neg_norm(iiVoxel,:)) > nanmean(allVoxelsDensities.rest_neg_norm(iiVoxel,:))
           stimulus_contrast_neg_increase(idxx,idxy,idxz) = p_stimulus_neg(iiVoxel);
        elseif nanmean(allVoxelsDensities.passive_neg_norm(iiVoxel,:)) < nanmean(allVoxelsDensities.rest_neg_norm(iiVoxel,:))
           stimulus_contrast_neg_decrease(idxx,idxy,idxz) = p_stimulus_neg(iiVoxel);
        end
        
        [H_stimulus_pos(iiVoxel),p_stimulus_pos(iiVoxel)] = ttest(squeeze(allVoxelsDensities.passive_pos_norm(iiVoxel,:)),squeeze(allVoxelsDensities.rest_pos_norm(iiVoxel,:))); 

        if nanmean(allVoxelsDensities.passive_pos_norm(iiVoxel,:)) > nanmean(allVoxelsDensities.rest_pos_norm(iiVoxel,:))
           stimulus_contrast_pos_increase(idxx,idxy,idxz) = p_stimulus_pos(iiVoxel);
        elseif nanmean(allVoxelsDensities.passive_pos_norm(iiVoxel,:)) < nanmean(allVoxelsDensities.rest_pos_norm(iiVoxel,:))
           stimulus_contrast_pos_decrease(idxx,idxy,idxz) = p_stimulus_pos(iiVoxel);
        end
        
    end
    
end

save('allVoxelsDensitiesGlobalNormalization.mat','allVoxelsDensities','H_attention_neg','H_attention_pos','H_stimulus_neg','H_stimulus_pos','p_attention_neg','p_attention_pos','p_stimulus_neg','p_stimulus_pos');

scl_slope = 1;
scl_inter = 0; 
dim = [size(AAL_img,1),size(AAL_img,2),size(AAL_img,3)];
dtype = 'FLOAT32';
offset = 0;
nifti_file.mat = load_aal.mat;
nifti_file.mat_intent = load_aal.mat_intent;
nifti_file.mat0 = load_aal.mat0;
nifti_file.mat0_intent = load_aal.mat0_intent;

fname = strcat('Global-Attention-Pos-Increase-GN','.nii');
descrip = 'Global-Attention-Pos-Increase';
input_data = attention_contrast_pos_increase; 
real_save_image;

fname = strcat('Global-Attention-Neg-Increase-GN','.nii');
descrip = 'Global-Attention-Neg-Increase';
input_data = attention_contrast_neg_increase; 
real_save_image;

fname = strcat('Global-Stimulus-Pos-Increase-GN','.nii');
descrip = 'Global-Stimulus-Pos-Increase';
input_data = stimulus_contrast_pos_increase; 
real_save_image;

fname = strcat('Global-Stimulus-Neg-Increase-GN','.nii');
descrip = 'Global-Stimulus-Neg-Increase';
input_data = stimulus_contrast_neg_increase; 
real_save_image;

fname = strcat('Global-Attention-Pos-Decrease-GN','.nii');
descrip = 'Global-Attention-Pos-Decrease';
input_data = attention_contrast_pos_decrease; 
real_save_image;

fname = strcat('Global-Attention-Neg-Decrease-GN','.nii');
descrip = 'Global-Attention-Neg-Decrease';
input_data = attention_contrast_neg_decrease; 
real_save_image;

fname = strcat('Global-Stimulus-Pos-Decrease-GN','.nii');
descrip = 'Global-Stimulus-Pos-Decrease';
input_data = stimulus_contrast_pos_decrease; 
real_save_image;

fname = strcat('Global-Stimulus-Neg-Decrease-GN','.nii');
descrip = 'Global-Stimulus-Neg-Decrease';
input_data = stimulus_contrast_neg_decrease; 
real_save_image;

end

function doContrastsGlobalNormalizationPerRun

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nROI = 90;
idx_ROI = 1:nROI;

load('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\GLOBAL\densities\outside-volumes\2nd\allVoxelsDensities.mat');

allVoxelsDensities.track_neg = allVoxelsDensities.track_neg(:,1:32);
allVoxelsDensities.track_pos = allVoxelsDensities.track_pos(:,1:32);

allVoxelsDensities.passive_neg = allVoxelsDensities.passive_neg(:,1:32);
allVoxelsDensities.passive_pos = allVoxelsDensities.passive_pos(:,1:32);

allVoxelsDensities.rest_neg = allVoxelsDensities.rest_neg(:,1:32);
allVoxelsDensities.rest_pos = allVoxelsDensities.rest_pos(:,1:32);

nSubjects = 8;
nRunsPerSubject = 4;
nRuns = 32;

allVoxelsDensities.track_neg_norm = [];
allVoxelsDensities.track_pos_norm = [];
allVoxelsDensities.passive_neg_norm = [];
allVoxelsDensities.passive_pos_norm = [];
allVoxelsDensities.rest_neg_norm = [];
allVoxelsDensities.rest_pos_norm = [];

for iRun=1:nRuns
    
    all_track_neg = allVoxelsDensities.track_neg(:,iRun);
    mean_all_track_neg = mean(all_track_neg(:));
    all_track_neg_norm = all_track_neg ./ mean_all_track_neg;
    
    allVoxelsDensities.track_neg_norm = [allVoxelsDensities.track_neg_norm,all_track_neg_norm];
    
    all_track_pos = allVoxelsDensities.track_pos(:,iRun);
    mean_all_track_pos = mean(all_track_pos(:));
    all_track_pos_norm = all_track_pos ./ mean_all_track_pos;
    
    allVoxelsDensities.track_pos_norm = [allVoxelsDensities.track_pos_norm,all_track_pos_norm];
    
    all_passive_neg = allVoxelsDensities.passive_neg(:,iRun);
    mean_all_passive_neg = mean(all_passive_neg(:));
    all_passive_neg_norm = all_passive_neg ./ mean_all_passive_neg;
    
    allVoxelsDensities.passive_neg_norm = [allVoxelsDensities.passive_neg_norm,all_passive_neg_norm];
    
    all_passive_pos = allVoxelsDensities.passive_pos(:,iRun);
    mean_all_passive_pos = mean(all_passive_pos(:));
    all_passive_pos_norm = all_passive_pos ./ mean_all_passive_pos;
    
    allVoxelsDensities.passive_pos_norm = [allVoxelsDensities.passive_pos_norm,all_passive_pos_norm];
    
    all_rest_neg = allVoxelsDensities.rest_neg(:,iRun);
    mean_all_rest_neg = mean(all_rest_neg(:));
    all_rest_neg_norm = all_rest_neg ./ mean_all_rest_neg;
    
    allVoxelsDensities.rest_neg_norm = [allVoxelsDensities.rest_neg_norm,all_rest_neg_norm];
    
    all_rest_pos = allVoxelsDensities.rest_pos(:,iRun);
    mean_all_rest_pos = mean(all_rest_pos(:));
    all_rest_pos_norm = all_rest_pos ./ mean_all_rest_pos;
    
    allVoxelsDensities.rest_pos_norm = [allVoxelsDensities.rest_pos_norm,all_rest_pos_norm];
    
end

attention_contrast_pos_decrease = zeros(size(AAL_img));
attention_contrast_neg_decrease = zeros(size(AAL_img));

attention_contrast_pos_increase = zeros(size(AAL_img));
attention_contrast_neg_increase = zeros(size(AAL_img));

stimulus_contrast_pos_decrease = zeros(size(AAL_img));
stimulus_contrast_neg_decrease = zeros(size(AAL_img));

stimulus_contrast_pos_increase = zeros(size(AAL_img));
stimulus_contrast_neg_increase = zeros(size(AAL_img));

nTotalVoxels = 160990;

H_attention_neg = zeros(1,nTotalVoxels);
H_attention_pos = zeros(1,nTotalVoxels);

H_stimulus_neg = zeros(1,nTotalVoxels);
H_stimulus_pos = zeros(1,nTotalVoxels);

p_attention_neg = zeros(1,nTotalVoxels);
p_attention_pos = zeros(1,nTotalVoxels);

p_stimulus_neg = zeros(1,nTotalVoxels);
p_stimulus_pos = zeros(1,nTotalVoxels);

iiVoxel = 0;

for iROI=1:nROI

    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    
    area_label = strrep(label_ROI{iROI},'_','-');       
    
    disp(area_label);
    
    idx_voxels = find(AAL_img==AAL_ROI(idx_ROI(iROI)).ID);
    
    nVoxels = length(idx_voxels);
    
    disp(strcat('nVoxes=',int2str(nVoxels)));
    
    for iVoxel=1:nVoxels
        
        iiVoxel = iiVoxel + 1;
        
        if mod(iVoxel,1000) == 0; disp(strcat('...voxels so far:',int2str(iVoxel))); end
       
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
    
        [H_attention_neg(iiVoxel),p_attention_neg(iiVoxel)] = ttest(squeeze(allVoxelsDensities.track_neg_norm(iiVoxel,:)),squeeze(allVoxelsDensities.passive_neg_norm(iiVoxel,:)));
        
        if nanmean(allVoxelsDensities.track_neg_norm(iiVoxel,:)) > nanmean(allVoxelsDensities.passive_neg_norm(iiVoxel,:))
           attention_contrast_neg_increase(idxx,idxy,idxz) = p_attention_neg(iiVoxel);
        elseif nanmean(allVoxelsDensities.track_neg_norm(iiVoxel,:)) < nanmean(allVoxelsDensities.passive_neg_norm(iiVoxel,:))
           attention_contrast_neg_decrease(idxx,idxy,idxz) = p_attention_neg(iiVoxel);
        end
        
        [H_attention_pos(iiVoxel),p_attention_pos(iiVoxel)] = ttest(squeeze(allVoxelsDensities.track_pos_norm(iiVoxel,:)),squeeze(allVoxelsDensities.passive_pos_norm(iiVoxel,:))); 
   
        if nanmean(allVoxelsDensities.track_pos_norm(iiVoxel,:)) > nanmean(allVoxelsDensities.passive_pos_norm(iiVoxel,:))
           attention_contrast_pos_increase(idxx,idxy,idxz) = p_attention_pos(iiVoxel);
        elseif nanmean(allVoxelsDensities.track_pos_norm(iiVoxel,:)) < nanmean(allVoxelsDensities.passive_pos_norm(iiVoxel,:))
           attention_contrast_pos_decrease(idxx,idxy,idxz) = p_attention_pos(iiVoxel);
        end
        
        [H_stimulus_neg(iiVoxel),p_stimulus_neg(iiVoxel)] = ttest(squeeze(allVoxelsDensities.passive_neg_norm(iiVoxel,:)),squeeze(allVoxelsDensities.rest_neg_norm(iiVoxel,:)));
        
        if nanmean(allVoxelsDensities.passive_neg_norm(iiVoxel,:)) > nanmean(allVoxelsDensities.rest_neg_norm(iiVoxel,:))
           stimulus_contrast_neg_increase(idxx,idxy,idxz) = p_stimulus_neg(iiVoxel);
        elseif nanmean(allVoxelsDensities.passive_neg_norm(iiVoxel,:)) < nanmean(allVoxelsDensities.rest_neg_norm(iiVoxel,:))
           stimulus_contrast_neg_decrease(idxx,idxy,idxz) = p_stimulus_neg(iiVoxel);
        end
        
        [H_stimulus_pos(iiVoxel),p_stimulus_pos(iiVoxel)] = ttest(squeeze(allVoxelsDensities.passive_pos_norm(iiVoxel,:)),squeeze(allVoxelsDensities.rest_pos_norm(iiVoxel,:))); 

        if nanmean(allVoxelsDensities.passive_pos_norm(iiVoxel,:)) > nanmean(allVoxelsDensities.rest_pos_norm(iiVoxel,:))
           stimulus_contrast_pos_increase(idxx,idxy,idxz) = p_stimulus_pos(iiVoxel);
        elseif nanmean(allVoxelsDensities.passive_pos_norm(iiVoxel,:)) < nanmean(allVoxelsDensities.rest_pos_norm(iiVoxel,:))
           stimulus_contrast_pos_decrease(idxx,idxy,idxz) = p_stimulus_pos(iiVoxel);
        end
        
    end
    
end

save('allVoxelsDensitiesGlobalNormalizationPerRun.mat','allVoxelsDensities','H_attention_neg','H_attention_pos','H_stimulus_neg','H_stimulus_pos','p_attention_neg','p_attention_pos','p_stimulus_neg','p_stimulus_pos');

scl_slope = 1;
scl_inter = 0; 
dim = [size(AAL_img,1),size(AAL_img,2),size(AAL_img,3)];
dtype = 'FLOAT32';
offset = 0;
nifti_file.mat = load_aal.mat;
nifti_file.mat_intent = load_aal.mat_intent;
nifti_file.mat0 = load_aal.mat0;
nifti_file.mat0_intent = load_aal.mat0_intent;

fname = strcat('Global-Attention-Pos-Increase-GN-PerRun','.nii');
descrip = 'Global-Attention-Pos-Increase';
input_data = attention_contrast_pos_increase; 
real_save_image;

fname = strcat('Global-Attention-Neg-Increase-GN-PerRun','.nii');
descrip = 'Global-Attention-Neg-Increase';
input_data = attention_contrast_neg_increase; 
real_save_image;

fname = strcat('Global-Stimulus-Pos-Increase-GN-PerRun','.nii');
descrip = 'Global-Stimulus-Pos-Increase';
input_data = stimulus_contrast_pos_increase; 
real_save_image;

fname = strcat('Global-Stimulus-Neg-Increase-GN-PerRun','.nii');
descrip = 'Global-Stimulus-Neg-Increase';
input_data = stimulus_contrast_neg_increase; 
real_save_image;

fname = strcat('Global-Attention-Pos-Decrease-GN-PerRun','.nii');
descrip = 'Global-Attention-Pos-Decrease';
input_data = attention_contrast_pos_decrease; 
real_save_image;

fname = strcat('Global-Attention-Neg-Decrease-GN-PerRun','.nii');
descrip = 'Global-Attention-Neg-Decrease';
input_data = attention_contrast_neg_decrease; 
real_save_image;

fname = strcat('Global-Stimulus-Pos-Decrease-GN-PerRun','.nii');
descrip = 'Global-Stimulus-Pos-Decrease';
input_data = stimulus_contrast_pos_decrease; 
real_save_image;

fname = strcat('Global-Stimulus-Neg-Decrease-GN-PerRun','.nii');
descrip = 'Global-Stimulus-Neg-Decrease';
input_data = stimulus_contrast_neg_decrease; 
real_save_image;

end


