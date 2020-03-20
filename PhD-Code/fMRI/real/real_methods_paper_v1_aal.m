function real_methods_paper_v1_aal

%%% AAL REGIONS

function getMeanTimeSeriesFrom90AALRunByRun

nROI = 90;
nTotalRuns = 32;
nRuns = 4;
nTR = 150;

all_settings = getAllSettings;

nSubjects = length(all_settings);

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

mean_area_RestingState_voxel = zeros(nROI,nTR);
mean_area_Passive_voxel = zeros(nROI,nTR);
mean_area_Track_voxel = zeros(nROI,nTR);

iiRun = 0;

for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
    
    kind = 'RestingState';
    for irun=1:nRuns;
        [RestingState(irun).run, RestingState(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
    kind = 'Passive';
    for irun=1:nRuns;
        [Passive(irun).run, Passive(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
    kind = 'Track';
    for irun=1:nRuns;
        [Track(irun).run, Track(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
    end
    
    for irun=1:nRuns
        
        iiRun = iiRun + 1;
        
        disp(strcat('Run:',int2str(iiRun)));
       
        for iROI=1:nROI

            idx_ROI = AAL_ROI(iROI).ID;
            idx_ROI_Label = AAL_ROI(iROI).Nom_L;
            idx_ROI_Label = strrep(idx_ROI_Label,'_','-');

            idx_voxels = find(AAL_img == idx_ROI);

            nVoxels = length(idx_voxels);

            area_RestingState = zeros(nVoxels,nTR);
            area_Passive = zeros(nVoxels,nTR);
            area_Track = zeros(nVoxels,nTR);

            izr = 0;
            izp = 0;
            izt = 0;

            zeros_RestingState = [];
            zeros_Passive = [];
            zeros_Track = [];

            for iVoxel=1:nVoxels

                [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

                area_RestingState(iVoxel,:) = RestingState(irun).run(idxx,idxy,idxz,:);
                area_Passive(iVoxel,:) = Passive(irun).run(idxx,idxy,idxz,:);
                area_Track(iVoxel,:) = Track(irun).run(idxx,idxy,idxz,:);

                if sum(area_RestingState(iVoxel,:)) == 0

                    izr = izr + 1;
                    zeros_RestingState(izr) = iVoxel;

                end
                
                if sum(area_Passive(iVoxel,:)) == 0

                    izp = izp + 1;
                    zeros_Passive(izp) = iVoxel;

                end
                
                if sum(area_Track(iVoxel,:)) == 0

                    izt = izt + 1;
                    zeros_Track(izt) = iVoxel;

                end


            end

            area_RestingState(zeros_RestingState,:) = [];
            area_Passive(zeros_Passive,:) = [];
            area_Track(zeros_Track,:) = [];

            mn_RestingState = mean(area_RestingState,1);
            mn_Passive = mean(area_Passive,1);
            mn_Track = mean(area_Track,1);

            mean_area_RestingState_voxel(iROI,:) = mn_RestingState(:);
            mean_area_Passive_voxel(iROI,:) = mn_Passive(:);
            mean_area_Track_voxel(iROI,:) = mn_Track(:);

        end
        
        mean_run(iiRun).rest = mean_area_RestingState_voxel;
        mean_run(iiRun).passive = mean_area_Passive_voxel;
        mean_run(iiRun).track = mean_area_Track_voxel;
        
        clear mean_area_RestingState_voxel
        clear mean_area_Passive_voxel
        clear mean_area_Track_voxel
        
    end 

end

save('AAL_ROI_mean_run.mat','mean_run');

end

function getMeanTimeSeriesFrom90AALRunByRunPETRA

nROI = 90;
nTotalRuns = 27;
nSegments = 3;
nTR = 150;

all_settings = getAllSettingsPetra;
nSubjects = length(all_settings);

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

run = 1;
all_settings = getAllSettingsPetra;

iRun = 0;
for iSet=1:length(all_settings)
    
    settings = all_settings(iSet).settings;

    get_at_this_preprocessed_step = settings.FSL.folders.melodic;
    
    file = settings.FSL.files.functional.residual;

    [residual, settings] = real_get_data_FSL_Petra(settings,run,file,get_at_this_preprocessed_step);
    
    for iSeg=1:nSegments
        
        iRun = iRun + 1;
        
        Subject(iRun).residual(:,:,:,:) = residual(:,:,:,(1+(iSeg-1)*nTR):(nTR+(iSeg-1)*nTR));
        
    end

end
    
for iRun=1:nTotalRuns

    disp(strcat('Run:',int2str(iRun)));

    for iROI=1:nROI

        idx_ROI = AAL_ROI(iROI).ID;
        idx_ROI_Label = AAL_ROI(iROI).Nom_L;
        idx_ROI_Label = strrep(idx_ROI_Label,'_','-');

        idx_voxels = find(AAL_img == idx_ROI);

        nVoxels = length(idx_voxels);

        area_RestingState = zeros(nVoxels,nTR);

        izr = 0;

        zeros_RestingState = [];

        for iVoxel=1:nVoxels

            [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

            area_RestingState(iVoxel,:) = Subject(iRun).residual(idxx,idxy,idxz,:);

            if sum(area_RestingState(iVoxel,:)) == 0

                izr = izr + 1;
                zeros_RestingState(izr) = iVoxel;

            end

        end

        area_RestingState(zeros_RestingState,:) = [];

        mn_RestingState = mean(area_RestingState,1);

        mean_area_RestingState_voxel(iROI,:) = mn_RestingState(:);

    end

    mean_run(iRun).rest = mean_area_RestingState_voxel;

    clear mean_area_RestingState_voxel

end 

save('AAL_ROI_mean_run-Petra.mat','mean_run');

end

function getMeanTimeSeriesFrom90AALRunByRunHCP

nROI = 90;
nTotalRuns = 38*4;
nTotalRunsQuarter = 38*4/4;
nRuns = 4;
nTR = 1200;

all_settings = getAllSettingsHCP;
nSubjects = length(all_settings);

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

iiRun = 0;
for iSet=1:5
    
    settings = all_settings(iSet).settings;
    
    disp(strcat('Subj:',int2str(iSet)));

    for iRun=1:nRuns
        
        iiRun = iiRun + 1;
        
        [all_runs(iiRun).run, settings] = real_get_data_HCP(settings,iRun);
    
    end
    
end

i_global_Run = 0;

for iRun=1:20
    
    i_global_Run = i_global_Run + 1;

    disp(strcat('Run:',int2str(iRun)));

    for iROI=1:nROI

        idx_ROI = AAL_ROI(iROI).ID;
        idx_ROI_Label = AAL_ROI(iROI).Nom_L;
        idx_ROI_Label = strrep(idx_ROI_Label,'_','-');

        idx_voxels = find(AAL_img == idx_ROI);

        nVoxels = length(idx_voxels);

        area_RestingState = zeros(nVoxels,nTR);

        izr = 0;

        zeros_RestingState = [];

        for iVoxel=1:nVoxels

            [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

            area_RestingState(iVoxel,:) = all_runs(iRun).run(idxx,idxy,idxz,:);

            if sum(area_RestingState(iVoxel,:)) == 0

                izr = izr + 1;
                zeros_RestingState(izr) = iVoxel;

            end

        end

        area_RestingState(zeros_RestingState,:) = [];

        mn_RestingState = mean(area_RestingState,1);

        mean_area_RestingState_voxel(iROI,:) = mn_RestingState(:);

    end

    mean_run(i_global_Run).rest = mean_area_RestingState_voxel;

    clear mean_area_RestingState_voxel

end 

iiRun = 0;
for iSet=6:8
    
    settings = all_settings(iSet).settings;
    
    disp(strcat('Subj:',int2str(iSet)));

    for iRun=1:nRuns
        
        iiRun = iiRun + 1;
        
        [all_runs(iiRun).run, settings] = real_get_data_HCP(settings,iRun);
    
    end
    
end

for iRun=1:12
    
    i_global_Run = i_global_Run + 1;

    disp(strcat('Run:',int2str(iRun)));

    for iROI=1:nROI

        idx_ROI = AAL_ROI(iROI).ID;
        idx_ROI_Label = AAL_ROI(iROI).Nom_L;
        idx_ROI_Label = strrep(idx_ROI_Label,'_','-');

        idx_voxels = find(AAL_img == idx_ROI);

        nVoxels = length(idx_voxels);

        area_RestingState = zeros(nVoxels,nTR);

        izr = 0;

        zeros_RestingState = [];

        for iVoxel=1:nVoxels

            [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));

            area_RestingState(iVoxel,:) = all_runs(iRun).run(idxx,idxy,idxz,:);

            if sum(area_RestingState(iVoxel,:)) == 0

                izr = izr + 1;
                zeros_RestingState(izr) = iVoxel;

            end

        end

        area_RestingState(zeros_RestingState,:) = [];

        mn_RestingState = mean(area_RestingState,1);

        mean_area_RestingState_voxel(iROI,:) = mn_RestingState(:);

    end

    mean_run(i_global_Run).rest = mean_area_RestingState_voxel;

    clear mean_area_RestingState_voxel

end 

save('AAL_ROI_mean_run-HCP.mat','mean_run');

end

function getCorrFromMeanTimeSeriesFrom90AALRunByRun

load('AAL_ROI_mean_run.mat');

nTotalRuns = 32;

for iRun=1:nTotalRuns
    
   run = mean_run(iRun).rest;
   run = run';
   rho(iRun).rest_corr = corr(run);
   
   run = mean_run(iRun).passive;
   run = run';
   rho(iRun).passive_corr = corr(run);
   
   run = mean_run(iRun).track;
   run = run';
   rho(iRun).track_corr = corr(run);
    
end

save('AAL_ROI_mean_run_corr.mat','rho');

end

function getCorrFromMeanTimeSeriesFrom90AALRunByRunPETRA

load('AAL_ROI_mean_run-Petra.mat');

nTotalRuns = 27;

for iRun=1:nTotalRuns
    
   run = mean_run(iRun).rest;
   run = run';
   rho(iRun).rest_corr = corr(run);
    
end

save('AAL_ROI_mean_run_corr-Petra.mat','rho');

end

function getCorrFromMeanTimeSeriesFrom90AALRunByRunHCP

load('AAL_ROI_mean_run-HCP.mat');

nTotalRuns = 32;

for iRun=1:nTotalRuns
    
   run = mean_run(iRun).rest;
   run = run';
   rho(iRun).rest_corr = corr(run);
    
end

save('AAL_ROI_mean_run_corr-HCP.mat','rho');

end

function getAndPlotVarianceDistributionFromFC90AAL

load('AAL_ROI_mean_run_corr.mat');

nROI = 90;
nTotalRuns = 32;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            end
            
        end
        
    end
    
end

% iPair = 0;
% 
% for iROI=1:nROI
%     
%     for iiROI=1:nROI
%         
%         if iROI ~= iiROI
%             
%             iPair = iPair + 1;
%         
%             for iRun=1:nTotalRuns
%             
%                 single_rho = rho(iRun).passive_corr(iROI,iiROI);
%             
%                 pairs_passive(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
%                 
%             end
%             
%         end
%         
%     end
%     
% end

% iPair = 0;
% 
% for iROI=1:nROI
%     
%     for iiROI=1:nROI
%         
%         if iROI ~= iiROI
%             
%             iPair = iPair + 1;
%         
%             for iRun=1:nTotalRuns
%             
%                 single_rho = rho(iRun).track_corr(iROI,iiROI);
%             
%                 pairs_track(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
%             
%             end
%             
%         end
%         
%     end
%     
% end

my_var_aal = var(pairs_rest);
my_cv_aal = std(pairs_rest) ./ abs(mean(pairs_rest,1));
my_mean_aal = mean(pairs_rest,1);
my_std_aal = std(pairs_rest);

vector_one = my_std_aal;
vector_two = my_mean_aal;
vector_three = my_cv_aal;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'AAL-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

fraction = FractionSignificantConsistent( vector_two, vector_one, z_threshold, cv_threshold )

% my_var_aal = var(pairs_passive);
% my_cv_aal = std(pairs_passive) ./ abs(mean(pairs_passive,1));
% my_mean_aal = mean(pairs_passive,1);
% my_std_aal = std(pairs_passive);
% 
% vector_one = my_std_aal;
% vector_two = my_mean_aal;
% vector_three = my_cv_aal;
% label_one = 'STD';
% label_two = 'Mean';
% label_three = 'CV';
% title_label = 'AAL-Passive';
% 
% MyContourCorrelation_v5(vector_two,vector_one,title_label,1);
% 
% my_var_aal = var(pairs_track);
% my_cv_aal = std(pairs_track) ./ abs(mean(pairs_track,1));
% my_mean_aal = mean(pairs_track,1);
% my_std_aal = std(pairs_track);
% 
% vector_one = my_std_aal;
% vector_two = my_mean_aal;
% vector_three = my_cv_aal;
% label_one = 'STD';
% label_two = 'Mean';
% label_three = 'CV';
% title_label = 'AAL-Track';
% 
% MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function getAndPlotVarianceDistributionFromFC90AALONESubject(nStartRun,nEndRun)

load('AAL_ROI_mean_run_corr.mat');

nROI = 90;
nTotalRuns = 4;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            iiRun = 0;
            for iRun=nStartRun:nEndRun
                
                iiRun = iiRun + 1;
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest(iiRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            end
            
        end
        
    end
    
end

my_var_aal = var(pairs_rest);
my_cv_aal = std(pairs_rest) ./ abs(mean(pairs_rest,1));
my_mean_aal = mean(pairs_rest,1);
my_std_aal = std(pairs_rest);

vector_one = my_std_aal;
vector_two = my_mean_aal;
vector_three = my_cv_aal;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = strcat('AAL-Rest','-Runs-',int2str(nStartRun),'-',int2str(nEndRun));

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

r_threshold = 0.186722;
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1.0;

fraction = FractionSignificantConsistent( vector_two, vector_one, z_threshold, cv_threshold )

end

function getAndPlotVarianceDistributionFromFC90AALPETRA

load('AAL_ROI_mean_run_corr-Petra.mat');

nROI = 90;
nTotalRuns = 27;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            end
            
        end
        
    end
    
end

my_var_aal = var(pairs_rest);
my_cv_aal = std(pairs_rest) ./ abs(mean(pairs_rest,1));
my_mean_aal = mean(pairs_rest,1);
my_std_aal = std(pairs_rest);

vector_one = my_std_aal;
vector_two = my_mean_aal;
vector_three = my_cv_aal;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'AAL-Petra-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function getAndPlotVarianceDistributionFromFC90AALHCP

load('AAL_ROI_mean_run_corr-HCP.mat');

nROI = 90;
nTotalRuns = 32;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
            
            end
            
        end
        
    end
    
end

my_var_aal = var(pairs_rest);
my_cv_aal = std(pairs_rest) ./ abs(mean(pairs_rest,1));
my_mean_aal = mean(pairs_rest,1);
my_std_aal = std(pairs_rest);

vector_one = my_std_aal;
vector_two = my_mean_aal;
vector_three = my_cv_aal;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'AAL-HCP-Rest';

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function getAndPlotVarianceDistributionFromFC90AALInsideVoxelLevel

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
nTotalRuns = 32;

all_settings = getAllSettings;

nSubjects = length(all_settings);

for iiRun=1:nTotalRuns
    
    aal_rhos(iiRun).rho = [];
    
end

for iROI=1:nROI

    label_ROI{iROI} = AAL_ROI(idx_ROI(iROI)).Nom_L;
    area_label = strrep(label_ROI{iROI},'_','-');       
    disp(area_label);
        
    iiRun = 0;
    for iSubject=1:length(all_settings)
        
        settings = all_settings(iSubject).settings;
  
        RestingStateRHOs = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI-corr','-','RestingState','-',area_label,'.mat'));
        
        for irun=1:nRuns
            
            iiRun = iiRun + 1;
       
            rest_rho = RestingStateRHOs.FC_Voxels.run(irun).rho_RestingState;
       
            aal_rhos(iiRun).rho = [aal_rhos(iiRun).rho(:)',rest_rho(:)'];
            
        end
     
    end

end

nRhos = length(aal_rhos(1).rho);

AAL_rhos = zeros(nTotalRuns,nRhos);

for iRun=1:nTotalRuns
    
    AAL_rhos(iRun,:) = aal_rhos(iRun).rho(:);

end

clear aal_rhos
            
AAL_rhos_z = (1/2).*log((1 + AAL_rhos)./(1 - AAL_rhos));

my_mean_aal = mean(AAL_rhos_z,1);
my_std_aal = std(AAL_rhos_z,0,1);

vector_one = my_std_aal;
vector_two = my_mean_aal;
label_one = 'STD';
label_two = 'Mean';
label_three = 'CV';
title_label = 'AAL-Inside-Rest';

r_threshold = 0.186728;  % corresponds to p-value = 0.01 for N=150
z_threshold = 0.5 * log( (1+r_threshold) ./ (1-r_threshold) );
cv_threshold = 1;

IM_total = TotalInformationContent( vector_two, vector_one, z_threshold, cv_threshold );

MyContourCorrelation_v5(vector_two,vector_one,title_label,1);

end

function plotMeanFCFrom90AAL

load('AAL_ROI_mean_run_corr.mat');

nROI = 90;
nTotalRuns = 32;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                pairs_rest_rho(iRun,iPair) = single_rho;
            
            end
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z,1);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ mean(pairs_rest_z,1);



mean_corr = zeros(nROI,nROI);
s_corr = zeros(nROI,nROI);

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            %if my_cv_aal(iPair) < cv_criterion & my_cv_aal(iPair) > 0 & p_pairs(iPair) < p_criterion
                
                mean_corr(iROI,iiROI) =  m_pairs_z(iPair);
                s_corr(iROI,iiROI) =  s_pairs_z(iPair);
                
            %end
            
        end
        
    end
    
end

label = 'FC-AAL';
plotFCgeneral(mean_corr,s_corr,label);

% f = figure;
% 
% clmap = colormap('jet');
% clmap(65,:) = clmap(64,:);
% clmap(33,:) = [1 1 1];
% 
% min_C = -1;
% max_C = 1;
% 
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% 
% print(f,'-depsc','FC-AAL.eps');
% export_fig('FC-AAL','-pdf');

end

function plotMeanFCFrom90AALONESubject(nStartRun,nEndRun)

load('AAL_ROI_mean_run_corr.mat');

nROI = 90;
nTotalRuns = 4;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            iiRun = 0;
            for iRun=nStartRun:nEndRun
                
                iiRun = iiRun + 1;
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest_z(iiRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                pairs_rest_rho(iiRun,iPair) = single_rho;
            
            end
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z,1);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ mean(pairs_rest_z,1);



mean_corr = zeros(nROI,nROI);
s_corr = zeros(nROI,nROI);

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            %if my_cv_aal(iPair) < cv_criterion & my_cv_aal(iPair) > 0 & p_pairs(iPair) < p_criterion
                
                mean_corr(iROI,iiROI) =  m_pairs_z(iPair);
                s_corr(iROI,iiROI) =  s_pairs_z(iPair);
                
            %end
            
        end
        
    end
    
end

label = strcat('FC-AAL','-Runs-',int2str(nStartRun),'-',int2str(nEndRun));
plotFCgeneral(mean_corr,s_corr,label);

% f = figure;
% 
% clmap = colormap('jet');
% clmap(65,:) = clmap(64,:);
% clmap(33,:) = [1 1 1];
% 
% min_C = -1;
% max_C = 1;
% 
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% 
% print(f,'-depsc','FC-AAL.eps');
% export_fig('FC-AAL','-pdf');

end

function plotMeanFCFrom90AALPETRA

load('AAL_ROI_mean_run_corr-Petra.mat');

nROI = 90;
nTotalRuns = 27;
nTR = 150;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                pairs_rest_rho(iRun,iPair) = single_rho;
            
            end
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z,1);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ mean(pairs_rest_z,1);



mean_corr = zeros(nROI,nROI);
s_corr = zeros(nROI,nROI);

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            %if my_cv_aal(iPair) < cv_criterion & my_cv_aal(iPair) > 0 & p_pairs(iPair) < p_criterion
                
                mean_corr(iROI,iiROI) =  m_pairs_z(iPair);
                s_corr(iROI,iiROI) =  s_pairs_z(iPair);
                
            %end
            
        end
        
    end
    
end

label = 'FC-AAL-Petra';
plotFCgeneral(mean_corr,s_corr,label);

% f = figure;
% 
% clmap = colormap('jet');
% clmap(65,:) = clmap(64,:);
% clmap(33,:) = [1 1 1];
% 
% min_C = -1;
% max_C = 1;
% 
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% 
% print(f,'-depsc','FC-AAL.eps');
% export_fig('FC-AAL','-pdf');

end

function plotMeanFCFrom90AALHCP

load('AAL_ROI_mean_run_corr-HCP.mat');

nROI = 90;
nTotalRuns = 32;
nTR = 1200;
cv_criterion = 0.75;
p_criterion = 0.01;

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            for iRun=1:nTotalRuns
                
                single_rho = rho(iRun).rest_corr(iROI,iiROI);
            
                pairs_rest_z(iRun,iPair) = (1/2)*log((1 + single_rho)/(1 - single_rho));
                
                pairs_rest_rho(iRun,iPair) = single_rho;
            
            end
            
        end
        
    end
    
end

m_pairs = mean(pairs_rest_rho,1);
m_pairs_z = mean(pairs_rest_z,1);

s_pairs_z = std(pairs_rest_z,1);

t_pairs = m_pairs ./ sqrt((1-m_pairs.^2)/(nTR-2));

p_pairs = 1-tcdf(t_pairs,nTR-1);

my_cv_aal = std(pairs_rest_z) ./ mean(pairs_rest_z,1);



mean_corr = zeros(nROI,nROI);
s_corr = zeros(nROI,nROI);

iPair = 0;

for iROI=1:nROI
    
    for iiROI=1:nROI
        
        if iROI ~= iiROI
            
            iPair = iPair + 1;
        
            %if my_cv_aal(iPair) < cv_criterion & my_cv_aal(iPair) > 0 & p_pairs(iPair) < p_criterion
                
                mean_corr(iROI,iiROI) =  m_pairs_z(iPair);
                s_corr(iROI,iiROI) =  s_pairs_z(iPair);
                
            %end
            
        end
        
    end
    
end

label = 'FC-AAL-HCP';
plotFCgeneral(mean_corr,s_corr,label);

% f = figure;
% 
% clmap = colormap('jet');
% clmap(65,:) = clmap(64,:);
% clmap(33,:) = [1 1 1];
% 
% min_C = -1;
% max_C = 1;
% 
% 
% imagesc(mean_corr);
% colorbar;
% caxis([min_C max_C]);
% colormap(clmap);
% 
% 
% print(f,'-depsc','FC-AAL.eps');
% export_fig('FC-AAL','-pdf');

end

end

