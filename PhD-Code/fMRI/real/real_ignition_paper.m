function real_ignition_paper

% doIgnition;

% level = 'cluster';
% pcriterion = 0.05;
% doIgnitionContrast(level,pcriterion);

% level = 'cluster';
% doIgnition(level,'All');

% level = 'voxel';
% Condition_Label = 'Track';
% doIgnition(level,Condition_Label);
% 
% level = 'voxel';
% Condition_Label = 'Passive';
% doIgnition(level,Condition_Label);
% 
% level = 'voxel';
% Condition_Label = 'RestingState';
% doIgnition(level,Condition_Label);

% saveIgnitionContrast3DVolume('cluster');

% saveIgnitionMEvokedInteg23DVolume('cluster');

% plotIgnitionStrength;

% level = 'cluster';
% pcriterion = 0.01;
% doIgnitionContrastSubjectLevel(level,pcriterion);

% level = 'cluster';
% pcriterion = 0.01;
% doIgnitionContrastSubjectLevelOnlyMax(level,pcriterion);

% level = 'cluster';
% which_kind = 'TS';
% getReportIgnitionContrastSubjectLevel(level,which_kind);

% level = 'cluster';
% Condition_Label = 'All';
% freq = [0.01 0.25];
% freq = [0.04 0.14];
% freq = [0.0001 0.002];
% freq = [0.04 0.07];
% pcriterion = 0.01;
% doIgnition(level,Condition_Label,freq);
% doIgnitionContrast(level,Condition_Label,pcriterion,freq);

% level = 'cluster';
% Condition_Label = 'All';
% doSpectralAnalysis(level,Condition_Label);

% level = 'cluster';
% freq = [0.04 0.07];
% pcriterion = 0.05;
% Condition_Label = 'All';
% doIgnitionPhases(level,'All',freq);
% doIgnitionContrast(level,Condition_Label,pcriterion,freq);

% pcriterion = 0.05;
% saveIgnitionContrast3DVolume('cluster',pcriterion);

% ROI = getMD758Clusters;
% level = 'cluster';
% freq = [0.04 0.07];
% pcriterion = 0.05;
% Condition_Label = 'All';
% selection_label = 'CORTEX';
% doIgnitionPhasesSelection(level,selection_label,Condition_Label,freq,ROI);
% doIgnitionContrastSelection(level,selection_label,Condition_Label,pcriterion,freq);
% selection_label = 'SUBCORTEX';
% doIgnitionPhasesSelection(level,selection_label,Condition_Label,freq,ROI);
% doIgnitionContrastSelection(level,selection_label,Condition_Label,pcriterion,freq);

% idx_cluster = 699;
% checkOverlapWithFunctionalNetworks(idx_cluster);

% level = 'cluster';
% freq = [0.04 0.07];
% pcriterion = 0.01;
% Condition_Label = 'All';
% doIgnitionContrastSTD(level,Condition_Label,pcriterion,freq);

% level = 'cluster';
% pcriterion = 0.01;
% saveIgnitionContrast3DVolume(level,pcriterion);

% level = 'cluster';
% freq = [0.04 0.07];
% pcriterion = 0.05;
% Condition_Label = 'All';
% doIgnition(level,Condition_Label,freq);

% saveIgnitionSTDEvokedInteg23DVolume('cluster');

% freq = [0.04 0.07];
% pcriterion = 0.01;
% Condition_Label = 'All';
% level = 'cluster';
% getPostTriggered3DVolume(level);
% doPostTriggeredContrast(level,Condition_Label,pcriterion,freq);
% getPostTriggeredContrast3DVolume(level,pcriterion);

% level = 'cluster';
% checkBCTparam(level);

% level = 'cluster';
% checkTRANSITIVITY_BU(level);
% checkMODULARITY_UND(level);
% checkSUBGRAPH_CENTRALITY(level);
% checkMODULE_DEGREE_ZSCORE(level);
% checkCLUSTERING_COEF_BU(level);
% checkSUBGRAPH_PEREVENT_CORBETTA(level);

% getSUBGRAPH_PHASE;

% plotGraphReliability;

plotGraphConcat;

% checkSUBGRAPH_PEREVENT_CORBETTA_cluster('cluster');

% doFreqMembershipMedianContrastPerSeed;

% doBinomialStatOnFreqMembershipContrastPerSeed;

% plotBinomialStatOnFreqMembershipContrastPerSeed;

% plotBinomialStatOnFreqMembershipContrastPerSeedGreenRedBlueGray;

% plotBinomialStatOnFreqMembershipContrastPerSeedGreenRedBlueGray;

% plotBinomialStatOnFreqMembershipContrast_v2;

% plot3DclustersFromContrast;

end

function Condition = loadAllRunsLHR(Condition_Label)

nRuns = 4;

all_settings = getAllSettings;
nSubjects = length(all_settings);

iiRun_track = 0;
iiRun_passive = 0;
iiRun_rest = 0;
for iSubject=1:nSubjects
    
    settings = all_settings(iSubject).settings;
    
    %% LOAD DATA
    get_at_this_preprocessed_step = settings.FSL.folders.custom;
    file = settings.FSL.files.functional.custom.residual_voxel;
    mask = settings.FSL.files.mask.custom;
    
    switch Condition_Label
        
        case 'Track'
            
            kind = 'Track';
            for irun=1:nRuns;
                iiRun_track = iiRun_track + 1;
                [Condition(iiRun_track).run, Condition(iiRun_track).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
            end
    
        case 'Passive'
            
            kind = 'Passive';
            for irun=1:nRuns;
                iiRun_passive = iiRun_passive + 1;
                [Condition(iiRun_passive).run, Condition(iiRun_passive).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
            end
    
        case 'RestingState'
            
            kind = 'RestingState';
            for irun=1:nRuns;
                iiRun_rest = iiRun_rest + 1;
                [Condition(iiRun_rest).run, Condition(iiRun_rest).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
            end
    
    end
    
end


end

function ROI = getMD758Clusters

load('Z:\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat');

end

function ROI_info = getROIInfoOnClusters

load('Z:\_DATA\Parcellation\758-Cluster\ROI_info.mat');

end

function [Events, Phases] = getEventsAndPhases(all_data,freq)

[nRun,nCond,nSeed,nTR] = size(all_data);

disp('...reshape');
all_data_rs = reshape(all_data,[nRun*nCond*nSeed,nTR]);
clear all_data
all_data_rs = all_data_rs';

TR = 2;

%%%%%%%%%%%%%
% flp = .04;                        % lowpass frequency of filter
% fhi = .07;                        % highpass
flp = freq(1);
fhi = freq(2);
delt = TR;                        % sampling interval
k = 2;                            % 2nd order butterworth filter
fnq = 1/(2*delt);                 % Nyquist frequency
min_freq = flp/fnq;
Wn = [flp/fnq (fhi/fnq)-min_freq];           % butterworth bandpass non-dimensional frequency
[bfilt2,afilt2] = butter(k,Wn);   % construct the filter

integration_cutoff = 10;
integration_cutoff2 = 9;

Tmax = nTR; 
T = integration_cutoff:Tmax-integration_cutoff;  

disp('...demean, detrend');
x = demean(detrend(all_data_rs),1);

clear all_data all_data_rs

disp('...zero-phase filter');
timeseriedata = filtfilt(bfilt2,afilt2,x);    % zero phase filter the data

clear x

disp('...Hilbert');
Xanalytic = hilbert(demean(timeseriedata));

disp('...detrend, demean');
tise = detrend(demean(timeseriedata,1));

clear timeseriedata

disp('...threshold events');
s_tise = std(tise,0,1);
m_tise = mean(tise,1);
ms_tise = s_tise + m_tise;
for iE=1:(nRun*nCond*nSeed)
    ev1(:,iE) = squeeze(tise(:,iE)) > ms_tise(iE);
end
ev2 = [zeros(1,nRun*nCond*nSeed); ev1(1:nTR-1,:)];

clear tise

event_contrast = (ev1-ev2)>0; 

event_contrast_tp = event_contrast';
z = num2cell(event_contrast_tp,2);

Xanalytic_tp = Xanalytic';
w = num2cell(Xanalytic_tp,2);

clear s_tise m_tise ms_tise ev1 ev2 event_contrast event_contrast_tp Xanalytic_tp

disp('...reshape again');
% transform back 
Events = zeros(nRun,nCond,nSeed,nTR);
iiSeed = 0;
for iRun=1:nRun
    for iCond=1:nCond
        for iSeed=1:nSeed
            iiSeed = iiSeed  + 1;
            seed = z{iiSeed};
            analytic = w{iiSeed};
            Events(iRun,iCond,iSeed,:) = seed(:);
            Phases(iRun,iCond,iSeed,:) = angle(analytic);
        end
    end
end

end

function all_data = getFormattedData(level,Condition_Label)

nTR = 150;
nTotalClusters = 758;
nSubjects = 8;
nROI = 90;
nRuns = 4;
nTotalRuns = 32;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking
TR = 2; %sampling interval(TR)  
nTotalVoxels = 160990;
MNI_size = [91 109 91];

disp(strcat('doing:',level));

switch level
    
    case 'voxel'
        
        disp('Load AAL');
        load_aal = nifti('ROI_MNI_V4.nii');
        load_aal.dat.fname = strcat('Z:\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
        AAL_img = load_aal.dat(:,:,:);
        load_roi = load('ROI_MNI_V4_List.mat');
        AAL_ROI = load_roi.ROI;
        
        Condition = loadAllRunsLHR(Condition_Label);
        
        N = nTotalVoxels; %758;%90;regions
        
        all_data = zeros(nTotalRuns,nConditions,N,nTR);
        
        iiRun = 0;
        for iSubject=1:nSubjects

            for iRun=1:nRuns

               iiRun = iiRun + 1;

               iiVoxel = 0;

               for iROI=1:nROI

                   label_ROI{iROI} = AAL_ROI(iROI).Nom_L;
                   area_label = strrep(label_ROI{iROI},'_','-');       
                   disp(strcat('Run:',int2str(iiRun),':',area_label));
                   
                   idx_ROI = AAL_ROI(iROI).ID;
                   
                   idx_voxels = find(AAL_img==idx_ROI);
                   
                   nVoxels = length(idx_voxels);

                   for iVoxel=1:nVoxels

                       iiVoxel = iiVoxel + 1;
                       
                       [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

                       all_data(iiRun,1,iiVoxel,1:nTR) = squeeze(Condition(iiRun).run(idxx,idxy,idxz,1:nTR));

                   end

               end

            end

        end
        
    case 'cluster' 
        
        ROI = getMD758Clusters;
        N = nTotalClusters; %758;%90;regions
        
        all_data = zeros(nTotalRuns,nConditions,N,nTR);

        iiRun = 0;
        for iSubject=1:nSubjects

            for iRun=1:nRuns

               iiRun = iiRun + 1;

               iiCluster = 0;

               for iROI=1:nROI

                   nClusters = ROI(iROI).nClusters;

                   for iCluster=1:nClusters

                       iiCluster = iiCluster + 1;

                       all_data(iiRun,1,iiCluster,1:nTR) = squeeze(ROI(iROI).clusters(iCluster).rest_mean(iiRun,1:nTR));
                       all_data(iiRun,2,iiCluster,1:nTR) = squeeze(ROI(iROI).clusters(iCluster).passive_mean(iiRun,1:nTR));
                       all_data(iiRun,3,iiCluster,1:nTR) = squeeze(ROI(iROI).clusters(iCluster).track_mean(iiRun,1:nTR));

                   end

               end

            end

        end
        
end

end

function all_data = getFormattedDataSelection(level,selection_label,Condition_Label,ROI)

nTR = 150;
nTotalClusters = 758;
nSubjects = 8;
nROI = 90;
nRuns = 4;
nTotalRuns = 32;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking
TR = 2; %sampling interval(TR)  
nTotalVoxels = 160990;
MNI_size = [91 109 91];

switch selection_label
         
    case 'CORTEX' 
        
        SELECTED_ROIS = 1:90;
        SELECTED_ROIS([31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78]) = [];
        
        N = 0;
        for iROI=SELECTED_ROIS
            
            N = N + ROI(iROI).nClusters;
            
        end
        
    case 'SUBCORTEX'
        
        SELECTED_ROIS = [31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];
        
        N = 0;
        for iROI=SELECTED_ROIS
            
            N = N + ROI(iROI).nClusters;
            
        end
        
end

disp(strcat('doing:',level));

switch level
         
    case 'cluster' 
        
        ROI = getMD758Clusters;
        % N = nTotalClusters;
        % %758;%90;regions;%NumberOfClustersInsideROISelection
        
        all_data = zeros(nTotalRuns,nConditions,N,nTR);

        iiRun = 0;
        for iSubject=1:nSubjects

            for iRun=1:nRuns

               iiRun = iiRun + 1;

               iiCluster = 0;

               for iROI=SELECTED_ROIS

                   nClusters = ROI(iROI).nClusters;

                   for iCluster=1:nClusters

                       iiCluster = iiCluster + 1;

                       all_data(iiRun,1,iiCluster,1:nTR) = squeeze(ROI(iROI).clusters(iCluster).rest_mean(iiRun,1:nTR));
                       all_data(iiRun,2,iiCluster,1:nTR) = squeeze(ROI(iROI).clusters(iCluster).passive_mean(iiRun,1:nTR));
                       all_data(iiRun,3,iiCluster,1:nTR) = squeeze(ROI(iROI).clusters(iCluster).track_mean(iiRun,1:nTR));

                   end

               end

            end

        end
        
end

end

function doIgnition(level,Condition_Label,freq)

load('Z:\_PAPERS\ignition\Events\analysis\Ignition-v1-cluster-Contrast.mat');

idx_att = find(contrast_attention_z);
idx_stim = find(contrast_stimulus_z);

all_idx_contrast = [idx_att,idx_stim];

nTR = 150;
nTotalClusters = 758;
nSubjects = 8;
nROI = 90;
nRuns = 4;
nTotalRuns = 32;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking
TR = 2; %sampling interval(TR)  
nTotalVoxels = 160990;
MNI_size = [91 109 91];

integration_cutoff = 10;
integration_cutoff2 = 9;

disp('getting Data...');

all_data = getFormattedData(level,Condition_Label);

N = size(all_data,3);

%%% EVENTS

disp('doing Events...');

[Events, ~] = getEventsAndPhases(all_data,freq);

clear all_data 

%%% INTEGRATION

disp('doing Integration...');
 
%  T=1:1:size(Phases,2);
T=integration_cutoff:1:(size(Events,4)-integration_cutoff);

for iRun=1:nTotalRuns
    
    for iCondition=1:nConditions

        for t = T
            
            if mod(t,10) == 0; disp(strcat('Integration:iRun:',int2str(iRun),':iCondition:',int2str(iCondition),':t:',int2str(t))); end

           %  phasematrix(i,j)=exp(-3*adif(Phases(i,t),Phases(j,t)));
           
           line_seeds = squeeze(Events(iRun,iCondition,:,t-integration_cutoff2))';
           column_seeds = squeeze(Events(iRun,iCondition,:,t-integration_cutoff2));
           
           line_seeds_matrix = repmat(line_seeds,[N,1]);
           column_seeds_matrix = repmat(column_seeds,[1,N]);
           
           Run(iRun).Cond(iCondition).T(t-integration_cutoff2).phasematrix(:,:) = line_seeds_matrix .* column_seeds_matrix;
           
           %%% Run(iRun).Cond(iCondition).T(t-integration_cutoff2).column_seeds_matrix = column_seeds_matrix;

           cc = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).phasematrix; %*Cbin;
           cc = cc - eye(N);
           
           if any(cc(:))

                [comps csize] = get_components(cc);
           
           else
               
               comps = [];
               csize = [];
               
           end
           
           Run(iRun).Cond(iCondition).T(t-integration_cutoff2).integ = max(csize)/N;
           Run(iRun).Cond(iCondition).T(t-integration_cutoff2).comps = comps;
           Run(iRun).Cond(iCondition).T(t-integration_cutoff2).csize = csize;

        end

    end
    
end


 %%%% EVENT TRIGGER
 
 disp('doing Event Trigger...');
 
 for iRun=1:nTotalRuns
     
     for iCondition=1:nConditions

         Run(iRun).Cond(iCondition).nevents = zeros(1,N);
         
         nevents = zeros(1,N);
         nevents2 = zeros(1,N);
         
         local_events = squeeze(Events(iRun,iCondition,:,:));

         for seed=1:N
             
           if mod(seed,10) == 0; disp(strcat('Trigger:iRun:',int2str(iRun),':iCondition:',int2str(iCondition),':seed:',int2str(seed))); end

           flag = 0;

           for t=T

               % if events(seed,t-9)==1 && flag==0
               if local_events(seed,t-9)==1 && flag==0
                   
                    flag = 1;
                    
                    Run(iRun).Cond(iCondition).nevents(seed) = Run(iRun).Cond(iCondition).nevents(seed) + 1;
                    
                    nevents(seed)= nevents(seed)+1;
                    nevents2(seed)= nevents2(seed)+1;
                    
               end
                
               if flag>0
                    
                   IntegStim(iRun,iCondition,seed,flag,nevents(seed)) = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).integ;
                   IntegStim2(iRun,iCondition,seed,flag,nevents2(seed)) = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).integ;
                   
                   IntegStim2Comps(iRun,iCondition,seed,flag,nevents(seed)).comps = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).comps;
                   IntegStim2Comps(iRun,iCondition,seed,flag,nevents(seed)).csize = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).csize;
                   %%% IntegStim2Comps(iRun,iCondition,seed,flag,nevents(seed)).phasematrix = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).phasematrix;
                   %%% IntegStim2Comps(iRun,iCondition,seed,flag,nevents(seed)).column_seeds_matrix = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).column_seeds_matrix;
                   
                   cc = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).phasematrix;
                   cc = cc - eye(N);
                   uam = cc;
                   comps = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).comps;

                   
                   
%                    if ~isempty(find(all_idx_contrast==seed))
% 
%                         STATUS_LABEL = strcat(int2str(iRun),':',int2str(iCondition),':',int2str(seed),':',int2str(t));
%  
%                         BCT(iRun,iCondition,seed,flag,nevents(seed)).allBCTparam = getBCTparam_v2( uam, comps, STATUS_LABEL );
% 
%                    end
%                 
                   flag  = flag + 1;
                
               end
                
               if flag==5
                
                   flag=0;
                   
               end
             
           end
           
         end
         
         for seed=1:N

            mevokedinteg2(iRun,iCondition,seed) = max(mean(squeeze(IntegStim2(iRun,iCondition,seed,:,1:nevents2(seed))),2));
             
         end
         
        mignitionon2(iRun,iCondition) = mean(squeeze(mevokedinteg2(iRun,iCondition,:)));
        stdignitionon2(iRun,iCondition) = std(squeeze(mevokedinteg2(iRun,iCondition,:)));

     end

 end

 disp('...doing allComps');
 
allComps = zeros(nTotalRuns,nConditions,N,N);
 
for iRun=1:nTotalRuns
    
    for iCondition=1:nConditions
   
        for seed=1:N
             
             comps = zeros(1,N);
             
             nEvent = squeeze(Run(iRun).Cond(iCondition).nevents(seed));
             nFlag = 4;
             
             iiFlagEvent = 0;
             
             tmp_graph = zeros(nEvent*nFlag,N);
             
             all_sizes = zeros(1,iiFlagEvent);
             
             for iFlag=1:nFlag
                 
                 for iEvent=1:nEvent
                     
                    iiFlagEvent = iiFlagEvent + 1;
                     
                    this_event_comps = IntegStim2Comps(iRun,iCondition,seed,iFlag,iEvent).comps;
                    
                    %%% this_event_phasematrix = IntegStim2Comps(iRun,iCondition,seed,iFlag,iEvent).phasematrix;
                    
                    %%% this_event_column_seeds_matrix = IntegStim2Comps(iRun,iCondition,seed,iFlag,iEvent).column_seeds_matrix;
                    
                    csize = IntegStim2Comps(iRun,iCondition,seed,iFlag,iEvent).csize;
                    
                    all_sizes(iiFlagEvent) = max(csize(:));
                    
                    if ~isempty(this_event_comps)
                     
                        [uComp, cComp] = count_unique(this_event_comps(:)); 
                    
                    else
                        
                        uComp = [];
                        
                    end
                    
                    if ~isempty(uComp)
                        
                        idx_biggest = find(cComp == max(cComp));
                        
                        if length(idx_biggest) == nTotalClusters
                            
                            idx_nodes = seed;
                            
                        else
                            
                            idx_component_biggest = uComp(idx_biggest);
                            
                            idx_nodes = find(IntegStim2Comps(iRun,iCondition,seed,iFlag,iEvent).comps == idx_component_biggest);
                    
                        end
                     
                        tmp_comp = zeros(1,N);
             
%                         for iIDX=1:length(idx_component_biggest)
% 
%                             idx_this_nodes = find(IntegStim2Comps(iRun,iCondition,seed,iFlag,iEvent).comps == idx_component_biggest(iIDX));
%                             
%                             idx_nodes = [idx_nodes,idx_this_nodes];
%                             
%                         end
                        
                        tmp_comp(idx_nodes) = 1;
                        
                        %%% biggest_subgraphs(iRun,iCondition,seed).column_seeds_matrix{iiFlagEvent} = this_event_column_seeds_matrix;
                        
                        %%% biggest_subgraphs(iRun,iCondition,seed).phasematrix{iiFlagEvent} = this_event_phasematrix;
                        
                        biggest_subgraphs(iRun,iCondition,seed).comps{iiFlagEvent} = this_event_comps;
                        
                        biggest_subgraphs(iRun,iCondition,seed).perEvent{iEvent,iFlag}.comps = this_event_comps;
                        
                        biggest_subgraphs(iRun,iCondition,seed).subgraph{iiFlagEvent} = idx_nodes;
                        
                        biggest_subgraphs(iRun,iCondition,seed).perEvent{iEvent,iFlag}.subgraph = idx_nodes;
                        
                        biggest_subgraphs(iRun,iCondition,seed).perEvent{iEvent,iFlag}.membership = tmp_comp;
                        
                        tmp_graph(iiFlagEvent,:) = tmp_comp;

                        comps = comps + tmp_comp;
                    
                    end
                    
                 end
                 
             end
             
             biggest_subgraphs(iRun,iCondition,seed).subsize = all_sizes;
             
             biggest_subgraphs(iRun,iCondition,seed).membership = tmp_graph;
             
             biggest_subgraphs(iRun,iCondition,seed).reliability = getReliability(tmp_graph);
             
             clear tmp_size tmp_graph
             
             norm_comps = comps ./ ( nFlag*nEvent );
             allComps(iRun,iCondition,seed,:) =  norm_comps;
          
        end
        
    end
    
end

mignitionon = mean(mignitionon2,1);
stdignitionon = mean(stdignitionon2,1);

%minteg=mean(integt);

for iRun=1:nTotalRuns
    
    for iCondition=1:nConditions
        
        for seed=1:N

            mevokedinteg(iRun,iCondition,seed,:) = mean(squeeze(IntegStim(iRun,iCondition,seed,:,1:nevents(seed))),2);
            mev(iRun,iCondition,seed) = mean(squeeze(max(squeeze(IntegStim(iRun,iCondition,seed,:,1:nevents(seed))),[],1)));
            stdev(iRun,iCondition,seed) = std(squeeze(max(squeeze(IntegStim(iRun,iCondition,seed,:,1:nevents(seed))),[],1)));
            % mevokedQ(iRun,iCondition,seed,:) = mean(squeeze(IntegStimQ(iRun,iCondition,seed,:,1:nevents(seed))),2);

        end

    end
    
end

disp('...saving');

save(strcat('Ignition-v6','-',level,'-',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'-allOthers.mat'), 'mev', 'stdev', 'mevokedinteg','mignitionon','stdignitionon','mignitionon2','stdignitionon2','mevokedinteg2','Events','IntegStim2','-v7.3');
save(strcat('Ignition-v6','-',level,'-',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'-Run.mat'), 'Run', '-v7.3');
save(strcat('Ignition-v6','-',level,'-',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'-IntegStim2Comps-allComps.mat'),'IntegStim2Comps','allComps','-v7.3');
save(strcat('Ignition-v6','-',level,'-',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'-biggest_subgraphs.mat'), 'biggest_subgraphs', '-v7.3');
% save(strcat('Ignition-v6','-',level,'-',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'-BCT.mat'), 'BCT', '-v7.3');

end

function doIgnitionContrast(level,Condition_Label,pcriterion,freq)

disp('Loading Ignition...');
load(strcat('Ignition-v1','-',level,'-',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'-Phases.mat'));

nTotalClusters = 758;
% pcriterion = 0.01;

N = nTotalClusters;

idx_RestingState = 1;
idx_PassiveViewing = 2;
idx_Attention = 3;

p_attention = zeros(1,N);
h_attention = zeros(1,N);
p_stimulus = zeros(1,N);
h_stimulus = zeros(1,N);

for iIg=1:N
    
   disp(strcat('Ig:',int2str(iIg)));
      
   p_attention(iIg) = ranksum(squeeze(mevokedinteg2(:,idx_Attention,iIg)),squeeze(mevokedinteg2(:,idx_PassiveViewing,iIg)));
   
   if median(squeeze(mevokedinteg2(:,idx_Attention,iIg))) > median(squeeze(mevokedinteg2(:,idx_PassiveViewing,iIg))); h_attention(iIg) = -1; end
   if median(squeeze(mevokedinteg2(:,idx_Attention,iIg))) < median(squeeze(mevokedinteg2(:,idx_PassiveViewing,iIg))); h_attention(iIg) = 1; end
  
   p_stimulus(iIg) = ranksum(squeeze(mevokedinteg2(:,idx_PassiveViewing,iIg)),squeeze(mevokedinteg2(:,idx_RestingState,iIg)));
   
   if median(squeeze(mevokedinteg2(:,idx_PassiveViewing,iIg))) > median(squeeze(mevokedinteg2(:,idx_RestingState,iIg))); h_stimulus(iIg) = -1; end
   if median(squeeze(mevokedinteg2(:,idx_PassiveViewing,iIg))) < median(squeeze(mevokedinteg2(:,idx_RestingState,iIg))); h_stimulus(iIg) = 1; end

end

contrast_attention = zeros(1,N);
contrast_stimulus = zeros(1,N);

contrast_attention_z = zeros(1,N);
contrast_stimulus_z = zeros(1,N);

for iIg=1:N
      
    contrast_attention(iIg) = p_attention(iIg);
    contrast_stimulus(iIg) = p_stimulus(iIg);
    
end

contrast_attention(contrast_attention>pcriterion) = 0;
contrast_stimulus(contrast_stimulus>pcriterion) = 0;

for iIg=1:N
     
    if contrast_attention(iIg) ~= 0; contrast_attention_z(iIg) = norminv(contrast_attention(iIg)); end
    if contrast_stimulus(iIg) ~= 0; contrast_stimulus_z(iIg) = norminv(contrast_stimulus(iIg)); end

    contrast_attention_z(iIg) = h_attention(iIg) * contrast_attention_z(iIg);
    contrast_stimulus_z(iIg) = h_stimulus(iIg) * contrast_stimulus_z(iIg);
    
end

save(strcat('Ignition-v1','-',level,'-','Contrast','-',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'-Phases.mat'),'contrast_stimulus_z','contrast_attention_z');

end

function doIgnitionContrastSTD(level,Condition_Label,pcriterion,freq)

disp('Loading Ignition...');
% load(strcat('Ignition-v1','-',level,'-',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'-Phases.mat'));
load(strcat('Ignition-v1','-',level,'-',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'.mat'));

nTotalClusters = 758;
% pcriterion = 0.01;
nTotalRuns = 32;
nConditions = 3;

N = nTotalClusters;

idx_RestingState = 1;
idx_PassiveViewing = 2;
idx_Attention = 3;

p_attention = zeros(1,N);
h_attention = zeros(1,N);
p_stimulus = zeros(1,N);
h_stimulus = zeros(1,N);

integration_cutoff = 10;
integration_cutoff2 = 9;
 
%  T=1:1:size(Phases,2);
T=integration_cutoff:1:(size(Events,4)-integration_cutoff);

sevokedinteg2 = zeros(nTotalRuns,nConditions,nTotalClusters);

%%%% EVENT TRIGGER
 
 disp('doing Event Trigger...');
 
 for iRun=1:nTotalRuns
     
     for iCondition=1:nConditions
 
         nevents = zeros(1,N);
         nevents2 = zeros(1,N);
         
         local_events = squeeze(Events(iRun,iCondition,:,:));

         for seed=1:N
             
           if mod(seed,1000) == 0; disp(strcat('Trigger:iRun:',int2str(iRun),':iCondition:',int2str(iCondition),':seed:',int2str(seed))); end

           flag = 0;

           for t=T

               % if events(seed,t-9)==1 && flag==0
               if local_events(seed,t-9)==1 && flag==0
                   
                    flag = 1;
                    nevents(seed)= nevents(seed)+1;
                    nevents2(seed)= nevents2(seed)+1;
                    
               end
                
               if flag>0
                    
                   IntegStim(iRun,iCondition,seed,flag,nevents(seed)) = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).integ;
                   IntegStim2(iRun,iCondition,seed,flag,nevents2(seed)) = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).integ;
                
                   flag  = flag + 1;
                
               end
                
               if flag==5
                
                   flag=0;
                   
               end
             
           end
           
         end

         for seed=1:N
            
             % mevokedinteg2(iRun,iCondition,seed) = max(mean(squeeze(IntegStim2(iRun,iCondition,seed,:,1:nevents2(seed))),2));
             sevokedinteg2(iRun,iCondition,seed) = max(std(squeeze(IntegStim2(iRun,iCondition,seed,:,1:nevents2(seed))),0,2));

         end
         
     end

end

for iIg=1:N
    
   disp(strcat('Ig:',int2str(iIg)));
      
   p_attention(iIg) = ranksum(squeeze(sevokedinteg2(:,idx_Attention,iIg)),squeeze(sevokedinteg2(:,idx_PassiveViewing,iIg)));
   
   if median(squeeze(sevokedinteg2(:,idx_Attention,iIg))) > median(squeeze(sevokedinteg2(:,idx_PassiveViewing,iIg))); h_attention(iIg) = -1; end
   if median(squeeze(sevokedinteg2(:,idx_Attention,iIg))) < median(squeeze(sevokedinteg2(:,idx_PassiveViewing,iIg))); h_attention(iIg) = 1; end
  
   p_stimulus(iIg) = ranksum(squeeze(sevokedinteg2(:,idx_PassiveViewing,iIg)),squeeze(sevokedinteg2(:,idx_RestingState,iIg)));
   
   if median(squeeze(sevokedinteg2(:,idx_PassiveViewing,iIg))) > median(squeeze(sevokedinteg2(:,idx_RestingState,iIg))); h_stimulus(iIg) = -1; end
   if median(squeeze(sevokedinteg2(:,idx_PassiveViewing,iIg))) < median(squeeze(sevokedinteg2(:,idx_RestingState,iIg))); h_stimulus(iIg) = 1; end

end

contrast_attention = zeros(1,N);
contrast_stimulus = zeros(1,N);

contrast_attention_z = zeros(1,N);
contrast_stimulus_z = zeros(1,N);

for iIg=1:N
      
    contrast_attention(iIg) = p_attention(iIg);
    contrast_stimulus(iIg) = p_stimulus(iIg);
    
end

contrast_attention(contrast_attention>pcriterion) = 0;
contrast_stimulus(contrast_stimulus>pcriterion) = 0;

for iIg=1:N
     
    if contrast_attention(iIg) ~= 0; contrast_attention_z(iIg) = norminv(contrast_attention(iIg)); end
    if contrast_stimulus(iIg) ~= 0; contrast_stimulus_z(iIg) = norminv(contrast_stimulus(iIg)); end

    contrast_attention_z(iIg) = h_attention(iIg) * contrast_attention_z(iIg);
    contrast_stimulus_z(iIg) = h_stimulus(iIg) * contrast_stimulus_z(iIg);
    
end

save(strcat('Ignition-v1','-',level,'-','Contrast-STD','-',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'.mat'),'contrast_stimulus_z','contrast_attention_z','sevokedinteg2');

end

function doIgnitionContrastSubjectLevel(level,pcriterion)

disp('Loading Ignition...');
load(strcat('Ignition-v1-',level,'-All.mat'));

nTotalClusters = 758;
% pcriterion = 0.01;
nRuns = 4;
nTotalRuns = 32;

N = nTotalClusters;
nSubjects = 8;

idx_RestingState = 1;
idx_PassiveViewing = 2;
idx_Attention = 3;

p_attention = zeros(1,N);
h_attention = zeros(1,N);
p_stimulus = zeros(1,N);
h_stimulus = zeros(1,N);

[nTotalRuns,nConditions,nSeeds,nFlags,nEvents] = size(IntegStim2);

r_IntegStim = reshape(IntegStim2,[nTotalRuns nConditions nSeeds nFlags*nEvents]);

iSubject = 0;
for iRun=1:nRuns:nTotalRuns
    
    iSubject = iSubject + 1;

    for iSeed=1:N

       disp(strcat('Subject:',int2str(iSubject),':Seed:',int2str(iSeed)));

       samples_attention = squeeze(mean(r_IntegStim(iRun:iRun+nRuns-1,idx_Attention,iSeed,:),1));
       samples_passive = squeeze(mean(r_IntegStim(iRun:iRun+nRuns-1,idx_PassiveViewing,iSeed,:),1));
       samples_rest = squeeze(mean(r_IntegStim(iRun:iRun+nRuns-1,idx_RestingState,iSeed,:),1));
       
       [s,I] = sort(samples_attention,'descend');
       samples_attention(I(1)) = [];
       samples_attention(samples_attention==0) = [];

       [s,I] = sort(samples_passive,'descend');
       samples_passive(I(1)) = [];
       samples_passive(samples_passive==0) = [];
       
       [s,I] = sort(samples_rest,'descend');
       samples_rest(I(1)) = [];
       samples_rest(samples_rest==0) = [];
       
       p_attention(iSubject,iSeed) = ranksum(samples_attention(:),samples_passive(:));

       if median(samples_attention(:)) > median(samples_passive(:)); h_attention(iSubject,iSeed) = -1; end
       if median(samples_attention(:)) < median(samples_passive(:)); h_attention(iSubject,iSeed) = 1; end

       p_stimulus(iSubject,iSeed) = ranksum(samples_passive(:),samples_rest(:));

       if median(samples_passive(:)) > median(samples_rest(:)); h_stimulus(iSubject,iSeed) = -1; end
       if median(samples_passive(:)) < median(samples_rest(:)); h_stimulus(iSubject,iSeed) = 1; end

    end

end

contrast_attention = zeros(iSubject,N);
contrast_stimulus = zeros(iSubject,N);

contrast_attention_z = zeros(iSubject,N);
contrast_stimulus_z = zeros(iSubject,N);

for iSubject=1:nSubjects
    
    for iSeed=1:N

        contrast_attention(iSubject,iSeed) = p_attention(iSubject,iSeed);
        contrast_stimulus(iSubject,iSeed) = p_stimulus(iSubject,iSeed);

    end

end

contrast_attention(contrast_attention>pcriterion) = 0;
contrast_stimulus(contrast_stimulus>pcriterion) = 0;

for iSubject=1:nSubjects
    
    for iSeed=1:N

        if contrast_attention(iSubject,iSeed) ~= 0; contrast_attention_z(iSubject,iSeed) = norminv(contrast_attention(iSubject,iSeed)); end
        if contrast_stimulus(iSubject,iSeed) ~= 0; contrast_stimulus_z(iSubject,iSeed) = norminv(contrast_stimulus(iSubject,iSeed)); end

        contrast_attention_z(iSubject,iSeed) = h_attention(iSubject,iSeed) * contrast_attention_z(iSubject,iSeed);
        contrast_stimulus_z(iSubject,iSeed) = h_stimulus(iSubject,iSeed) * contrast_stimulus_z(iSubject,iSeed);

    end

end

save(strcat('Ignition-v1-',level,'-Contrast-',num2str(pcriterion),'-Subjects','.mat'),'contrast_stimulus_z','contrast_attention_z');

end

function doIgnitionContrastSubjectLevelOnlyMax(level,pcriterion)

disp('Loading Ignition...');
load(strcat('Ignition-v1-',level,'-All.mat'));

nTotalClusters = 758;
% pcriterion = 0.01;
nRuns = 4;
nTotalRuns = 32;

N = nTotalClusters;
nSubjects = 8;

idx_RestingState = 1;
idx_PassiveViewing = 2;
idx_Attention = 3;

p_attention = zeros(1,N);
h_attention = zeros(1,N);
p_stimulus = zeros(1,N);
h_stimulus = zeros(1,N);

[nTotalRuns,nConditions,nSeeds,nFlags,nEvents] = size(IntegStim2);

r_IntegStim = reshape(IntegStim2,[nTotalRuns nConditions nSeeds nFlags*nEvents]);

iSubject = 0;
for iRun=1:nRuns:nTotalRuns
    
    iSubject = iSubject + 1;

    for iSeed=1:N

       disp(strcat('Subject:',int2str(iSubject),':Seed:',int2str(iSeed)));

       samples_attention = squeeze(max(r_IntegStim(iRun:iRun+nRuns-1,idx_Attention,iSeed,:),[],4));
       samples_passive = squeeze(max(r_IntegStim(iRun:iRun+nRuns-1,idx_PassiveViewing,iSeed,:),[],4));
       samples_rest = squeeze(max(r_IntegStim(iRun:iRun+nRuns-1,idx_RestingState,iSeed,:),[],4));

       p_attention(iSubject,iSeed) = ranksum(samples_attention(:),samples_passive(:));

       if median(samples_attention(:)) > median(samples_passive(:)); h_attention(iSubject,iSeed) = -1; end
       if median(samples_attention(:)) < median(samples_passive(:)); h_attention(iSubject,iSeed) = 1; end

       p_stimulus(iSubject,iSeed) = ranksum(samples_passive(:),samples_rest(:));

       if median(samples_passive(:)) > median(samples_rest(:)); h_stimulus(iSubject,iSeed) = -1; end
       if median(samples_passive(:)) < median(samples_rest(:)); h_stimulus(iSubject,iSeed) = 1; end

    end

end

contrast_attention = zeros(iSubject,N);
contrast_stimulus = zeros(iSubject,N);

contrast_attention_z = zeros(iSubject,N);
contrast_stimulus_z = zeros(iSubject,N);

for iSubject=1:nSubjects
    
    for iSeed=1:N

        contrast_attention(iSubject,iSeed) = p_attention(iSubject,iSeed);
        contrast_stimulus(iSubject,iSeed) = p_stimulus(iSubject,iSeed);

    end

end

contrast_attention(contrast_attention>pcriterion) = 0;
contrast_stimulus(contrast_stimulus>pcriterion) = 0;

for iSubject=1:nSubjects
    
    for iSeed=1:N

        if contrast_attention(iSubject,iSeed) ~= 0; contrast_attention_z(iSubject,iSeed) = norminv(contrast_attention(iSubject,iSeed)); end
        if contrast_stimulus(iSubject,iSeed) ~= 0; contrast_stimulus_z(iSubject,iSeed) = norminv(contrast_stimulus(iSubject,iSeed)); end

        contrast_attention_z(iSubject,iSeed) = h_attention(iSubject,iSeed) * contrast_attention_z(iSubject,iSeed);
        contrast_stimulus_z(iSubject,iSeed) = h_stimulus(iSubject,iSeed) * contrast_stimulus_z(iSubject,iSeed);

    end

end

save(strcat('Ignition-v1-',level,'-Contrast-',num2str(pcriterion),'-Subjects-OnlyMax','.mat'),'contrast_stimulus_z','contrast_attention_z');

end

function saveIgnitionContrast3DVolume(level,pcriterion)


% load(strcat('Ignition-v1-',level,'-Contrast-All-0.04-0.07-Phases-p',num2str(pcriterion),'.mat'));
load(strcat('Ignition-v1-',level,'-Contrast-STD-All-0.04-0.07-p',num2str(pcriterion),'.mat'));

func_vol = nifti('LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
func_vol.dat.fname = 'Z:\_DATA\Parcellation\758-Cluster\LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii';

func_clu = func_vol.dat(:,:,:);

z_threshold = abs ( -sqrt(2) * erfcinv(pcriterion*2) );
MNI_size = [91 109 91];

%% Attention

idx_increase_att = find(contrast_attention_z>z_threshold);
idx_decrease_att = find(contrast_attention_z<((-1)*z_threshold));

idx_att_all = [idx_increase_att(:);idx_decrease_att(:)];

att_increase_vol = zeros(size(func_clu));
att_decrease_vol = zeros(size(func_clu));

att_all_vol = zeros(size(func_clu));
att_all_vol_idx = zeros(size(func_clu));

for idx=idx_increase_att
    
    idx_voxels = find(func_clu==idx);
    
    nVoxels = length(idx_voxels);
    
    for iVoxel=1:nVoxels
        
        [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
        
        att_increase_vol(idxx,idxy,idxz) = contrast_attention_z(idx);
        
        att_all_vol_idx(idxx,idxy,idxz) = idx;
        
    end
    
end

for idx=idx_decrease_att
    
    idx_voxels = find(func_clu==idx);
    
    nVoxels = length(idx_voxels);
    
    for iVoxel=1:nVoxels
        
        [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
        
        att_decrease_vol(idxx,idxy,idxz) = contrast_attention_z(idx);
        
        att_all_vol_idx(idxx,idxy,idxz) = idx;
        
    end
    
end

for idx=idx_att_all'
    
    idx_voxels = find(func_clu==idx);
    
    nVoxels = length(idx_voxels);
    
    for iVoxel=1:nVoxels
        
        [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
        
        att_all_vol(idxx,idxy,idxz) = contrast_attention_z(idx);
        
    end
    
end

%% Stimulus

idx_increase_stim = find(contrast_stimulus_z>z_threshold);
idx_decrease_stim = find(contrast_stimulus_z<((-1)*z_threshold));

idx_stim_all = [idx_increase_stim(:);idx_decrease_stim(:)];

stim_increase_vol = zeros(size(func_clu));
stim_decrease_vol = zeros(size(func_clu));

stim_all_vol = zeros(size(func_clu));
stim_all_vol_idx = zeros(size(func_clu));

for idx=idx_increase_stim
    
    idx_voxels = find(func_clu==idx);
    
    nVoxels = length(idx_voxels);
    
    for iVoxel=1:nVoxels
        
        [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
        
        stim_increase_vol(idxx,idxy,idxz) = contrast_stimulus_z(idx);
        
        stim_all_vol_idx(idxx,idxy,idxz) = idx;
        
    end
    
end

for idx=idx_decrease_stim
    
    idx_voxels = find(func_clu==idx);
    
    nVoxels = length(idx_voxels);
    
    for iVoxel=1:nVoxels
        
        [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
        
        stim_decrease_vol(idxx,idxy,idxz) = contrast_stimulus_z(idx);
        
        stim_all_vol_idx(idxx,idxy,idxz) = idx;
        
    end
    
end

for idx=idx_stim_all'
    
    idx_voxels = find(func_clu==idx);
    
    nVoxels = length(idx_voxels);
    
    for iVoxel=1:nVoxels
        
        [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
        
        stim_all_vol(idxx,idxy,idxz) = contrast_stimulus_z(idx);
        
    end
    
end

pcriterion = round(pcriterion * 100);

nifti_file = func_vol;
offset = func_vol.dat.offset;
scl_slope = func_vol.dat.scl_slope;
scl_inter = func_vol.dat.scl_inter;
dtype = 'FLOAT32';
offset = 0;
dim = func_vol.dat.dim;

descrip = 'Attention-Increase';
fname = strcat('Ignition-v1-',level,'-Contrast-','Attention-Increase-',num2str(pcriterion),'.nii');
input_data = att_increase_vol; 
real_save_image;

descrip = 'Attention-Decrease';
fname = strcat('Ignition-v1-',level,'-Contrast-','Attention-Decrease-',num2str(pcriterion),'.nii');
input_data = att_decrease_vol; 
real_save_image;

descrip = 'Attention-All';
fname = strcat('Ignition-v1-',level,'-Contrast-','Attention-All-',num2str(pcriterion),'.nii');
input_data = att_all_vol; 
real_save_image;

descrip = 'Attention-All-idx';
fname = strcat('Ignition-v1-',level,'-Contrast-','Attention-All-idx-',num2str(pcriterion),'.nii');
input_data = att_all_vol_idx; 
real_save_image;

descrip = 'Stimulus-Increase';
fname = strcat('Ignition-v1-',level,'-Contrast-','Stimulus-Increase-',num2str(pcriterion),'.nii');
input_data = stim_increase_vol; 
real_save_image;

descrip = 'Stimulus-Decrease';
fname = strcat('Ignition-v1-',level,'-Contrast-','Stimulus-Decrease-',num2str(pcriterion),'.nii');
input_data = stim_decrease_vol; 
real_save_image;

descrip = 'Stimulus-All';
fname = strcat('Ignition-v1-',level,'-Contrast-','Stimulus-All-',num2str(pcriterion),'.nii');
input_data = stim_all_vol; 
real_save_image;

descrip = 'Stimulus-All-idx';
fname = strcat('Ignition-v1-',level,'-Contrast-','Stimulus-All-idx-',num2str(pcriterion),'.nii');
input_data = stim_all_vol_idx; 
real_save_image;

end

function saveIgnitionMEvokedInteg23DVolume(level)

disp('Loading Ignition...');
load(strcat('Ignition-v1-',level,'-All-0.04-0.07-Phases.mat'));

MNI_size = [91 109 91];

func_vol = nifti('LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
func_vol.dat.fname = 'Z:\_DATA\Parcellation\758-Cluster\LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii';

func_clu = func_vol.dat(:,:,:);

nTotalClusters = 758;

idx_RestingState = 1;
idx_PassiveViewing = 2;
idx_Attention = 3;

m_mevokedinteg2 = squeeze(mean(mevokedinteg2,1));

m_attention = zeros(MNI_size);
m_passive = zeros(MNI_size);
m_resting = zeros(MNI_size);

for iCluster=1:nTotalClusters
    
   idx_voxels = find(func_clu==iCluster);
   
   nVoxels = length(idx_voxels);
   
   for iVoxel=1:nVoxels
      
       [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
       
       m_attention(idxx,idxy,idxz) = m_mevokedinteg2(idx_Attention,iCluster); 
       m_passive(idxx,idxy,idxz) = m_mevokedinteg2(idx_PassiveViewing,iCluster);
       m_resting(idxx,idxy,idxz) = m_mevokedinteg2(idx_RestingState,iCluster);
       
   end
    
end

nifti_file = func_vol;
offset = func_vol.dat.offset;
scl_slope = func_vol.dat.scl_slope;
scl_inter = func_vol.dat.scl_inter;
dtype = 'FLOAT32';
offset = 0;
dim = func_vol.dat.dim;

descrip = 'Attention';
fname = strcat('Ignition-v1-',level,'-Mean-Attention','-All-0.04-0.07-Phases.nii');
input_data = m_attention; 
real_save_image;

descrip = 'Passive';
fname = strcat('Ignition-v1-',level,'-Mean-Passive','-All-0.04-0.07-Phases.nii');
input_data = m_passive; 
real_save_image;

descrip = 'RestingState';
fname = strcat('Ignition-v1-',level,'-Mean-RestingState','-All-0.04-0.07-Phases.nii');
input_data = m_resting; 
real_save_image;


end

function saveIgnitionSTDEvokedInteg23DVolume(level)

disp('Loading Ignition...');
% load(strcat('Ignition-v1-',level,'-All-0.04-0.07.mat'));
load(strcat('Ignition-v1-',level,'-Contrast-STD-All-0.04-0.07.mat'));

MNI_size = [91 109 91];

func_vol = nifti('LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
func_vol.dat.fname = 'Z:\_DATA\Parcellation\758-Cluster\LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii';

func_clu = func_vol.dat(:,:,:);

nTotalClusters = 758;

idx_RestingState = 1;
idx_PassiveViewing = 2;
idx_Attention = 3;

s_sevokedinteg2 = squeeze(mean(sevokedinteg2,1));

s_attention = zeros(MNI_size);
s_passive = zeros(MNI_size);
s_resting = zeros(MNI_size);

for iCluster=1:nTotalClusters
    
   idx_voxels = find(func_clu==iCluster);
   
   nVoxels = length(idx_voxels);
   
   for iVoxel=1:nVoxels
      
       [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
       
       s_attention(idxx,idxy,idxz) = s_sevokedinteg2(idx_Attention,iCluster); 
       s_passive(idxx,idxy,idxz) = s_sevokedinteg2(idx_PassiveViewing,iCluster);
       s_resting(idxx,idxy,idxz) = s_sevokedinteg2(idx_RestingState,iCluster);
       
   end
    
end

nifti_file = func_vol;
offset = func_vol.dat.offset;
scl_slope = func_vol.dat.scl_slope;
scl_inter = func_vol.dat.scl_inter;
dtype = 'FLOAT32';
offset = 0;
dim = func_vol.dat.dim;

descrip = 'Attention';
fname = strcat('Ignition-v1-',level,'-STD-Attention','-All-0.04-0.07.nii');
input_data = s_attention; 
real_save_image;

descrip = 'Passive';
fname = strcat('Ignition-v1-',level,'-STD-Passive','-All-0.04-0.07.nii');
input_data = s_passive; 
real_save_image;

descrip = 'RestingState';
fname = strcat('Ignition-v1-',level,'-STD-RestingState','-All-0.04-0.07.nii');
input_data = s_resting; 
real_save_image;


end

function getReportIgnitionContrastSubjectLevel(level,which_kind)

load(strcat('Ignition-v1-',level,'-Contrast-0.01-Subjects-',which_kind,'.mat'));

ROI_info = getROIInfoOnClusters;

nSubjects = 8;
z_threshold = 2.3;
nContrasts = 4;
nROIs = 90;

contrast_label = {'Attention-Pos' 'Attention-Neg' 'Stimulus-Pos' 'Stimulus-Neg'};

for iSubject=1:nSubjects
    
   attention_contrast = squeeze(contrast_attention_z(iSubject,:));
   stimulus_contrast = squeeze(contrast_stimulus_z(iSubject,:));
    
   idx_attention_pos = find(attention_contrast>z_threshold);
   idx_attention_neg = find(attention_contrast<(-1)*z_threshold);
   
   idx_stimulus_pos = find(stimulus_contrast>z_threshold);
   idx_stimulus_neg = find(stimulus_contrast<(-1)*z_threshold);

    for iContrast=1:nContrasts

        switch contrast_label{iContrast}

            case 'Attention-Pos'

                idx_contrast = idx_attention_pos;

            case 'Attention-Neg'

                idx_contrast = idx_attention_neg;

            case 'Stimulus-Pos'

                idx_contrast = idx_stimulus_pos;

            case 'Stimulus-Neg'

                idx_contrast = idx_stimulus_neg;

        end
        
        if strfind(contrast_label{iContrast},'Att')
        
            % report_attention{iSubject+1,1} = strcat('SUBJECT:',int2str(iSubject));
            
            for idx=1:length(idx_contrast)
                
                idx_cluster = idx_contrast(idx);
                
                report_attention(iSubject,idx_cluster) = attention_contrast(idx_cluster);

            end
            
        end     

        if strfind(contrast_label{iContrast},'Stim')
        
            % report_stimulus{iSubject+1,1} = strcat('SUBJECT:',int2str(iSubject));
            
            for idx=1:length(idx_contrast)
                
                idx_cluster = idx_contrast(idx);
                
                report_stimulus(iSubject,idx_cluster) = attention_contrast(idx_cluster);

            end
            
        end     

    end
   
 end

for iROI=1:nROIs

   idx_start = ROI_info{iROI,4};
   idx_end = ROI_info{iROI,5}; 
   iROI_label = ROI_info{iROI,2};
   
   for idx=idx_start:idx_end
      
       report_attention_labels{1,idx+1} = iROI_label;
       report_stimulus_labels{1,idx+1} = iROI_label;
       
   end

end          

save(strcat('Ignition-v1-',level,'-Contrast-Subjects-',which_kind,'-Report','.mat'),'report_attention','report_stimulus','report_attention_labels','report_stimulus_labels');

end

function plotIgnitionStrength

load('Ignition-v1-cluster-All.mat');
load('Ignition-v1-cluster-Contrast.mat');

ROI_info = getROIInfoOnClusters;

z_threshold = 2.3;
nRuns = 4;

[nTotalRuns,nConditions,nSeeds,nFlags,nEvents] = size(IntegStim2);

r_IntegStim = reshape(IntegStim2,[nTotalRuns nConditions nSeeds nFlags*nEvents]);

mr_IntegStim =squeeze(mean(r_IntegStim(:,:,:,:),1));

max_Integ = max(mr_IntegStim(:));

nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking
nROIs = 90;
nContrasts = 4;

Condition_label = {'RestingState' 'PassiveViewing' 'AttentiveTracking'};

color_cond{1} = 'b-';
color_cond{2} = 'r-';
color_cond{3} = 'k-';

idx_attention_pos = find(contrast_attention_z>z_threshold);
idx_attention_neg = find(contrast_attention_z<(-1)*z_threshold);

idx_stimulus_pos = find(contrast_stimulus_z>z_threshold);
idx_stimulus_neg = find(contrast_stimulus_z<(-1)*z_threshold);

contrast_label = {'Attention-Pos' 'Attention-Neg' 'Stimulus-Pos' 'Stimulus-Neg'};

for iContrast=1:nContrasts
    
    switch contrast_label{iContrast}
        
        case 'Attention-Pos'
            
            idx_contrast = idx_attention_pos;
            
        case 'Attention-Neg'
            
            idx_contrast = idx_attention_neg;
            
        case 'Stimulus-Pos'
            
            idx_contrast = idx_stimulus_pos;
            
        case 'Stimulus-Neg'
            
            idx_contrast = idx_stimulus_neg;
            
    end
    
    if ~isempty(idx_contrast)

        for iSeed=idx_contrast

            for iROI=1:nROIs

                if iSeed >= ROI_info{iROI,4} && iSeed <= ROI_info{iROI,5}; seed_AAL_label = ROI_info{iROI,2}; end

            end

            %%% MEAN OF ALL RUNS
            
            f = figure;

            for iCond=1:nConditions

                plot(squeeze(mr_IntegStim(iCond,iSeed,:)),color_cond{iCond});
                hold on

            end

            title(strcat(contrast_label{iContrast},':',seed_AAL_label,':',int2str(iSeed)));
            legend(Condition_label);
            xlabel('Sequence Of Events');
            ylabel('Max Comp Size of Phase Matrix');
            xlim([0 (nEvents*nFlags)+1]);
            ylim([0 max_Integ]);

            print(f,'-depsc',strcat('Ignition','-',contrast_label{iContrast},'-',seed_AAL_label,'-',int2str(iSeed),'.eps'));
            
            %%% ALL RUNS, AT SUBJECT LEVEL
            
            f = figure;
            
            for iRun=1:4:nTotalRuns
                
                for iCond=1:nConditions

                    plot(squeeze(mean(r_IntegStim(iRun:iRun+nRuns-1,iCond,iSeed,:),1)),color_cond{iCond});
                    hold on

                end
            
            end

            title(strcat('Subjects',':',contrast_label{iContrast},':',seed_AAL_label,':',int2str(iSeed)));
            xlabel('Sequence Of Events');
            ylabel('Max Comp Size of Phase Matrix');
            xlim([0 (nEvents*nFlags)+1]);
            ylim([0 max_Integ]);

            print(f,'-depsc',strcat('Ignition','-',contrast_label{iContrast},'-',seed_AAL_label,'-',int2str(iSeed),'-Subjects','.eps'));
            
            %%% ALL RUNS, ONE BY ONE
            
            f = figure;
            
            for iRun=1:nTotalRuns
                
                for iCond=1:nConditions

                    plot(squeeze(r_IntegStim(iRun,iCond,iSeed,:)),color_cond{iCond});
                    hold on

                end
            
            end

            title(strcat('AllRuns',':',contrast_label{iContrast},':',seed_AAL_label,':',int2str(iSeed)));
            xlabel('Sequence Of Events');
            ylabel('Max Comp Size of Phase Matrix');
            xlim([0 (nEvents*nFlags)+1]);
            ylim([0 max_Integ]);

            print(f,'-depsc',strcat('Ignition','-',contrast_label{iContrast},'-',seed_AAL_label,'-',int2str(iSeed),'-AllRuns','.eps'));

        end

    end

end

end

function doSpectralAnalysis(level,Condition_Label)

disp('...loading data');
all_data = getFormattedData(level,Condition_Label);

[nRun,nCond,nSeed,nTR] = size(all_data);

disp('...reshape');
all_data_rs = reshape(all_data,[nRun*nCond*nSeed,nTR]);
clear all_data
all_data_rs = all_data_rs';

N = nTR;
all_data_FFT = fft(all_data_rs,[],2);
all_data_xdft = all_data_FFT(1:N/2+1,:);

Fs = 1/2;
all_data_psdx = (1/(Fs*N)) * abs(all_data_xdft).^2;
all_data_psdx(2:end-1,:) = 2*all_data_psdx(2:end-1,:);
freq = 0:Fs/N:Fs/2;

db_all_data_psdx = 10*log10(all_data_psdx);
s_db_all_data_psdx = sum(db_all_data_psdx,2);

for iSeed=1:nRun*nCond*nSeed
    
    [pxx(iSeed,:), w(iSeed,:)] = periodogram(squeeze(all_data_rs(:,iSeed)),[],[],Fs);
    
end

s_pxx = sum(pxx,1);
w = squeeze(w(1,:));

med_pxx = max(s_pxx)/2;

idx_mean = find(s_pxx>med_pxx);

mean_pxx = s_pxx(idx_mean);
mean_w = w(idx_mean);

f = figure;
plot(mean_w,mean_pxx);
xlabel('Frequencies');
ylabel('Power');
title('Periodogram');



end

function doIgnitionPhases(level,Condition_Label,freq)

nTR = 150;
nTotalClusters = 758;
nSubjects = 8;
nROI = 90;
nRuns = 4;
nTotalRuns = 32;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking
TR = 2; %sampling interval(TR)  
nTotalVoxels = 160990;
MNI_size = [91 109 91];

integration_cutoff = 10;
integration_cutoff2 = 9;

disp('getting Data...');

all_data = getFormattedData(level,Condition_Label);

N = size(all_data,3);

%%% EVENTS

disp('doing Events...');

[Events, Phases] = getEventsAndPhases(all_data,freq);

clear all_data 

%%% INTEGRATION

disp('doing Integration...');
 
%  T=1:1:size(Phases,2);
T=integration_cutoff:1:(size(Phases,4)-integration_cutoff);

ttotal(1:nConditions) = 1;

for iRun=1:nTotalRuns
    
    for iCondition=1:nConditions

        for t = T
            
            if mod(t,10) == 0; disp(strcat('Integration:iRun:',int2str(iRun),':iCondition:',int2str(iCondition),':t:',int2str(t))); end

%            line_seeds = squeeze(Events(iRun,iCondition,:,t-integration_cutoff2))';
%            column_seeds = squeeze(Events(iRun,iCondition,:,t-integration_cutoff2));
%            
%            line_seeds_matrix = repmat(line_seeds,[N,1]);
%            column_seeds_matrix = repmat(column_seeds,[1,N]);
%            
%            Run(iRun).Cond(iCondition).T(t-integration_cutoff2).phasematrix(:,:) = line_seeds_matrix .* column_seeds_matrix;

           for iSeed=1:N
               for jSeed=1:N
                    Run(iRun).Cond(iCondition).T(t-integration_cutoff2).phasematrix(iSeed,jSeed) = exp(-3*adif(Phases(iRun,iCondition,iSeed,t),Phases(iRun,iCondition,jSeed,t)));
               end
           end             

           cc = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).phasematrix; %*Cbin;
           cc = cc - eye(N);

           % [comps csize] = get_components(cc);
           
           % Run(iRun).Cond(iCondition).T(t-integration_cutoff2).integ = max(csize)/N; 
           
          pp=1;
          PR=0:0.01:0.99;
          for p=PR
            A = abs(cc) > p;
            [comps, csize] = get_components(A);
            cs(pp) = max(csize);
            pp = pp+1;
          end
          Run(iRun).Cond(iCondition).T(t-integration_cutoff2).integ = sum(cs)*0.01/N;

          Cond(iCondition).integt(ttotal(iCondition)) = sum(cs)*0.01/N;
          
          [M modularity] = community_louvain(Run(iRun).Cond(iCondition).T(t-integration_cutoff2).phasematrix);
          
          Run(iRun).Cond(iCondition).T(t-integration_cutoff2).M = M;
          Run(iRun).Cond(iCondition).T(t-integration_cutoff2).modularity = modularity;
          
          ttotal(iCondition) = ttotal(iCondition) + 1;

        end

    end
    
end


 %%%% EVENT TRIGGER
 
 disp('doing Event Trigger...');
 
 for iRun=1:nTotalRuns
     
     for iCondition=1:nConditions
 
         nevents = zeros(1,N);
         nevents2 = zeros(1,N);
         
         local_events = squeeze(Events(iRun,iCondition,:,:));

         for seed=1:N
             
           if mod(seed,1000) == 0; disp(strcat('Trigger:iRun:',int2str(iRun),':iCondition:',int2str(iCondition),':seed:',int2str(seed))); end

           flag = 0;

           for t=T

               % if events(seed,t-9)==1 && flag==0
               if local_events(seed,t-9)==1 && flag==0
                   
                    flag = 1;
                    nevents(seed)= nevents(seed)+1;
                    nevents2(seed)= nevents2(seed)+1;
                    
               end
                
               if flag>0
                    
                   IntegStim(iRun,iCondition,seed,flag,nevents(seed)) = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).integ;
                   IntegStim2(iRun,iCondition,seed,flag,nevents2(seed)) = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).integ;
                   
                   IntegStimQ(iRun,iCondition,seed,flag,nevents(seed)) = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).modularity;
                
                   flag  = flag + 1;
                
               end
                
               if flag==5
                
                   flag=0;
                   
               end
             
           end
           
         end

         for seed=1:N
            
             mevokedinteg2(iRun,iCondition,seed) = max(mean(squeeze(IntegStim2(iRun,iCondition,seed,:,1:nevents2(seed))),2));
          
         end
         
         mignitionon2(iRun,iCondition) = mean(squeeze(mevokedinteg2(iRun,iCondition,:)));
         stdignitionon2(iRun,iCondition) = std(squeeze(mevokedinteg2(iRun,iCondition,:)));
         
     end

end

mignitionon = mean(mignitionon2,1);
stdignitionon = mean(stdignitionon2,1);

%minteg=mean(integt);

for iRun=1:nTotalRuns
    
    for iCondition=1:nConditions
        
        for seed=1:N

            mevokedinteg(iRun,iCondition,seed,:) = mean(squeeze(IntegStim(iRun,iCondition,seed,:,1:nevents(seed))),2);
            mev(iRun,iCondition,seed) = mean(squeeze(max(squeeze(IntegStim(iRun,iCondition,seed,:,1:nevents(seed))),[],1)));
            stdev(iRun,iCondition,seed) = std(squeeze(max(squeeze(IntegStim(iRun,iCondition,seed,:,1:nevents(seed))),[],1)));
            mevokedQ(iRun,iCondition,seed,:) = mean(squeeze(IntegStimQ(iRun,iCondition,seed,:,1:nevents(seed))),2);

        end

    end
    
end

save (strcat('Ignition-v1','-',level,'-',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'-Phases.mat'), 'mev', 'stdev', 'mevokedinteg','mevokedQ','mignitionon','stdignitionon','mignitionon2','stdignitionon2','mevokedinteg2','Run','Cond','Events','Phases','IntegStim2','-v7.3');

end

function doIgnitionPhasesSelection(level,selection_label,Condition_Label,freq,ROI)

nTR = 150;
nTotalClusters = 758;
nSubjects = 8;
nROI = 90;
nRuns = 4;
nTotalRuns = 32;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking
TR = 2; %sampling interval(TR)  
nTotalVoxels = 160990;
MNI_size = [91 109 91];

integration_cutoff = 10;
integration_cutoff2 = 9;

disp('getting Data...');

all_data = getFormattedDataSelection(level,selection_label,Condition_Label,ROI);

N = size(all_data,3);

%%% EVENTS

disp('doing Events...');

[Events, Phases] = getEventsAndPhases(all_data,freq);

clear all_data 

%%% INTEGRATION

disp('doing Integration...');
 
%  T=1:1:size(Phases,2);
T=integration_cutoff:1:(size(Phases,4)-integration_cutoff);

ttotal(1:nConditions) = 1;

for iRun=1:nTotalRuns
    
    for iCondition=1:nConditions

        for t = T
            
            if mod(t,10) == 0; disp(strcat('Integration:iRun:',int2str(iRun),':iCondition:',int2str(iCondition),':t:',int2str(t))); end

%            line_seeds = squeeze(Events(iRun,iCondition,:,t-integration_cutoff2))';
%            column_seeds = squeeze(Events(iRun,iCondition,:,t-integration_cutoff2));
%            
%            line_seeds_matrix = repmat(line_seeds,[N,1]);
%            column_seeds_matrix = repmat(column_seeds,[1,N]);
%            
%            Run(iRun).Cond(iCondition).T(t-integration_cutoff2).phasematrix(:,:) = line_seeds_matrix .* column_seeds_matrix;

           for iSeed=1:N
               for jSeed=1:N
                    Run(iRun).Cond(iCondition).T(t-integration_cutoff2).phasematrix(iSeed,jSeed) = exp(-3*adif(Phases(iRun,iCondition,iSeed,t),Phases(iRun,iCondition,jSeed,t)));
               end
           end             

           cc = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).phasematrix; %*Cbin;
           cc = cc - eye(N);

           % [comps csize] = get_components(cc);
           
           % Run(iRun).Cond(iCondition).T(t-integration_cutoff2).integ = max(csize)/N; 
           
          pp=1;
          PR=0:0.01:0.99;
          for p=PR
            A = abs(cc) > p;
            [comps, csize] = get_components(A);
            cs(pp) = max(csize);
            pp = pp+1;
          end
          Run(iRun).Cond(iCondition).T(t-integration_cutoff2).integ = sum(cs)*0.01/N;

          Cond(iCondition).integt(ttotal(iCondition)) = sum(cs)*0.01/N;
          
          [M modularity] = community_louvain(Run(iRun).Cond(iCondition).T(t-integration_cutoff2).phasematrix);
          
          Run(iRun).Cond(iCondition).T(t-integration_cutoff2).M = M;
          Run(iRun).Cond(iCondition).T(t-integration_cutoff2).modularity = modularity;
          
          ttotal(iCondition) = ttotal(iCondition) + 1;

        end

    end
    
end


 %%%% EVENT TRIGGER
 
 disp('doing Event Trigger...');
 
 for iRun=1:nTotalRuns
     
     for iCondition=1:nConditions
 
         nevents = zeros(1,N);
         nevents2 = zeros(1,N);
         
         local_events = squeeze(Events(iRun,iCondition,:,:));

         for seed=1:N
             
           if mod(seed,1000) == 0; disp(strcat('Trigger:iRun:',int2str(iRun),':iCondition:',int2str(iCondition),':seed:',int2str(seed))); end

           flag = 0;

           for t=T

               % if events(seed,t-9)==1 && flag==0
               if local_events(seed,t-9)==1 && flag==0
                   
                    flag = 1;
                    nevents(seed)= nevents(seed)+1;
                    nevents2(seed)= nevents2(seed)+1;
                    
               end
                
               if flag>0
                    
                   IntegStim(iRun,iCondition,seed,flag,nevents(seed)) = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).integ;
                   IntegStim2(iRun,iCondition,seed,flag,nevents2(seed)) = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).integ;
                   
                   IntegStimQ(iRun,iCondition,seed,flag,nevents(seed)) = Run(iRun).Cond(iCondition).T(t-integration_cutoff2).modularity;
                
                   flag  = flag + 1;
                
               end
                
               if flag==5
                
                   flag=0;
                   
               end
             
           end
           
         end

         for seed=1:N
            
             mevokedinteg2(iRun,iCondition,seed) = max(mean(squeeze(IntegStim2(iRun,iCondition,seed,:,1:nevents2(seed))),2));
          
         end
         
         mignitionon2(iRun,iCondition) = mean(squeeze(mevokedinteg2(iRun,iCondition,:)));
         stdignitionon2(iRun,iCondition) = std(squeeze(mevokedinteg2(iRun,iCondition,:)));
         
     end

end

mignitionon = mean(mignitionon2,1);
stdignitionon = mean(stdignitionon2,1);

%minteg=mean(integt);

for iRun=1:nTotalRuns
    
    for iCondition=1:nConditions
        
        for seed=1:N

            mevokedinteg(iRun,iCondition,seed,:) = mean(squeeze(IntegStim(iRun,iCondition,seed,:,1:nevents(seed))),2);
            mev(iRun,iCondition,seed) = mean(squeeze(max(squeeze(IntegStim(iRun,iCondition,seed,:,1:nevents(seed))),[],1)));
            stdev(iRun,iCondition,seed) = std(squeeze(max(squeeze(IntegStim(iRun,iCondition,seed,:,1:nevents(seed))),[],1)));
            mevokedQ(iRun,iCondition,seed,:) = mean(squeeze(IntegStimQ(iRun,iCondition,seed,:,1:nevents(seed))),2);

        end

    end
    
end

save (strcat('Ignition-v1','-',level,'-',selection_label,'-','-',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'-Phases.mat'), 'mev', 'stdev', 'mevokedinteg','mevokedQ','mignitionon','stdignitionon','mignitionon2','stdignitionon2','mevokedinteg2','Run','Cond','Events','Phases','IntegStim2','-v7.3');

end

function doIgnitionContrastSelection(level,selection_label,Condition_Label,pcriterion,freq)

disp('Loading Ignition...');
load(strcat('Ignition-v1','-',level,'-',selection_label,'--',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'-Phases.mat'));

nTotalClusters = 758;
% pcriterion = 0.01;

% N = nTotalClusters;
N = size(mevokedinteg2,3);

idx_RestingState = 1;
idx_PassiveViewing = 2;
idx_Attention = 3;

p_attention = zeros(1,N);
h_attention = zeros(1,N);
p_stimulus = zeros(1,N);
h_stimulus = zeros(1,N);

for iIg=1:N
    
   disp(strcat('Ig:',int2str(iIg)));
      
   p_attention(iIg) = ranksum(squeeze(mevokedinteg2(:,idx_Attention,iIg)),squeeze(mevokedinteg2(:,idx_PassiveViewing,iIg)));
   
   if median(squeeze(mevokedinteg2(:,idx_Attention,iIg))) > median(squeeze(mevokedinteg2(:,idx_PassiveViewing,iIg))); h_attention(iIg) = -1; end
   if median(squeeze(mevokedinteg2(:,idx_Attention,iIg))) < median(squeeze(mevokedinteg2(:,idx_PassiveViewing,iIg))); h_attention(iIg) = 1; end
  
   p_stimulus(iIg) = ranksum(squeeze(mevokedinteg2(:,idx_PassiveViewing,iIg)),squeeze(mevokedinteg2(:,idx_RestingState,iIg)));
   
   if median(squeeze(mevokedinteg2(:,idx_PassiveViewing,iIg))) > median(squeeze(mevokedinteg2(:,idx_RestingState,iIg))); h_stimulus(iIg) = -1; end
   if median(squeeze(mevokedinteg2(:,idx_PassiveViewing,iIg))) < median(squeeze(mevokedinteg2(:,idx_RestingState,iIg))); h_stimulus(iIg) = 1; end

end

contrast_attention = zeros(1,N);
contrast_stimulus = zeros(1,N);

contrast_attention_z = zeros(1,N);
contrast_stimulus_z = zeros(1,N);

for iIg=1:N
      
    contrast_attention(iIg) = p_attention(iIg);
    contrast_stimulus(iIg) = p_stimulus(iIg);
    
end

contrast_attention(contrast_attention>pcriterion) = 0;
contrast_stimulus(contrast_stimulus>pcriterion) = 0;

for iIg=1:N
     
    if contrast_attention(iIg) ~= 0; contrast_attention_z(iIg) = norminv(contrast_attention(iIg)); end
    if contrast_stimulus(iIg) ~= 0; contrast_stimulus_z(iIg) = norminv(contrast_stimulus(iIg)); end

    contrast_attention_z(iIg) = h_attention(iIg) * contrast_attention_z(iIg);
    contrast_stimulus_z(iIg) = h_stimulus(iIg) * contrast_stimulus_z(iIg);
    
end

save(strcat('Ignition-v1','-',level,'-',selection_label,'-','Contrast','-',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'-Phases.mat'),'contrast_stimulus_z','contrast_attention_z');

end

function checkOverlapWithFunctionalNetworks(idx_cluster)

MD758 = nifti('LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
MD758.dat.fname = ('Z:\_DATA\Parcellation\758-Cluster\LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
MD758 = MD758.dat(:,:,:);

idx_MD = find(MD758==idx_cluster);

folder = 'Z:\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually';

networks = {'DAN' 'VAN' 'VIS' 'AUD' 'FPC' 'SMN' 'DMN' 'LAN'};
    
DAN = nifti('DAN-bin.nii');
DAN.dat.fname = 'Z:\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually\DAN-bin.nii';
DAN = DAN.dat(:,:,:);
idx_DAN = find(DAN);
disp(strcat('DAN:',int2str(length(find(ismember(idx_MD,idx_DAN)))),':voxels'));

VAN = nifti('VAN-bin.nii');
VAN.dat.fname = 'Z:\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually\VAN-bin.nii';
VAN = VAN.dat(:,:,:);
idx_VAN = find(VAN);
disp(strcat('VAN:',int2str(length(find(ismember(idx_MD,idx_VAN)))),':voxels'));

VIS = nifti('VIS-bin.nii');
VIS.dat.fname = 'Z:\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually\VIS-bin.nii';
VIS = VIS.dat(:,:,:);
idx_VIS = find(VIS);
disp(strcat('VIS:',int2str(length(find(ismember(idx_MD,idx_VIS)))),':voxels'));

AUD = nifti('AUD-bin.nii');
AUD.dat.fname = 'Z:\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually\AUD-bin.nii';
AUD = AUD.dat(:,:,:);
idx_AUD = find(AUD);
disp(strcat('AUD:',int2str(length(find(ismember(idx_MD,idx_AUD)))),':voxels'));

LAN = nifti('LAN-bin.nii');
LAN.dat.fname = 'Z:\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually\LAN-bin.nii';
LAN = LAN.dat(:,:,:);
idx_LAN = find(LAN);
disp(strcat('LAN:',int2str(length(find(ismember(idx_MD,idx_LAN)))),':voxels'));

DMN = nifti('DMN-bin.nii');
DMN.dat.fname = 'Z:\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually\DMN-bin.nii';
DMN = DMN.dat(:,:,:);
idx_DMN = find(DMN);
disp(strcat('DMN:',int2str(length(find(ismember(idx_MD,idx_DMN)))),':voxels'));

SMN = nifti('SMN-bin.nii');
SMN.dat.fname = 'Z:\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually\SMN-bin.nii';
SMN = SMN.dat(:,:,:);
idx_SMN = find(SMN);
disp(strcat('SMN:',int2str(length(find(ismember(idx_MD,idx_SMN)))),':voxels'));

FPC = nifti('FPC-bin.nii');
FPC.dat.fname = 'Z:\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually\FPC-bin.nii';
FPC = FPC.dat(:,:,:);
idx_FPC = find(FPC);
disp(strcat('FPC:',int2str(length(find(ismember(idx_MD,idx_FPC)))),':voxels'));
    
    
end

function getPostTriggered3DVolume(level)

func_vol = nifti('LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
func_vol.dat.fname = 'Z:\_DATA\Parcellation\758-Cluster\LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii';
func_clu = func_vol.dat(:,:,:);

load(strcat('Ignition-v3-',level,'-All-0.04-0.07-IntegStim2Comps-allComps.mat'));

load('Z:\_PAPERS\ignition\Events\analysis\Ignition-v1-cluster-Contrast.mat');

nTR = 150;
nTotalClusters = 758;
nSubjects = 8;
nROI = 90;
nRuns = 4;
nTotalRuns = 32;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking
TR = 2; %sampling interval(TR)  
nTotalVoxels = 160990;
MNI_size = [91 109 91];
nFlags = 4;

idx_attention = find(contrast_attention_z);
% idx_stimulus = find(contrast_stimulus_z);

for idx=idx_attention
    
    % rs_comps = squeeze(allComps(:,1,idx,:));
    pv_comps = squeeze(allComps(:,2,idx,:));
    at_comps = squeeze(allComps(:,3,idx,:));
    
    % m_rs_comps = mean(rs_comps(:,:),1);
    m_pv_comps = mean(pv_comps(:,:),1);
    m_at_comps = mean(at_comps(:,:),1);
    
    % m_rs_comps_vol = zeros(MNI_size);
    m_pv_comps_vol = zeros(MNI_size);
    m_at_comps_vol = zeros(MNI_size);
    
    for iClu=1:nTotalClusters
        
        idx_voxels = find(func_clu == iClu);
        
        nVoxels = length(idx_voxels);
    
        for iVoxel=1:nVoxels

            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

            %m_rs_comps_vol(idxx,idxy,idxz) = m_rs_comps(iClu);
            m_pv_comps_vol(idxx,idxy,idxz) = m_pv_comps(iClu);
            m_at_comps_vol(idxx,idxy,idxz) = m_at_comps(iClu);

        end

    end
    
    nifti_file = func_vol;
    offset = func_vol.dat.offset;
    scl_slope = func_vol.dat.scl_slope;
    scl_inter = func_vol.dat.scl_inter;
    dtype = 'FLOAT32';
    offset = 0;
    dim = func_vol.dat.dim;

    descrip = 'Attention-Triggered';
    fname = strcat('Ignition-v3-',level,'-TriggeredBy-',int2str(idx),'-onAttContrast-Att','.nii');
    input_data = m_at_comps_vol; 
    real_save_image;
    
    descrip = 'Attention-Triggered';
    fname = strcat('Ignition-v3-',level,'-TriggeredBy-',int2str(idx),'-onAttContrast-Pass','.nii');
    input_data = m_pv_comps_vol; 
    real_save_image;
    
%     descrip = 'Attention-Triggered';
%     fname = strcat('Ignition-v3-',level,'-TriggeredBy-',int2str(idx),'-onAttContrast-Rest','.nii');
%     input_data = m_rs_comps_vol; 
%     real_save_image;
    
end
    
%idx_attention = find(contrast_attention_z);
idx_stimulus = find(contrast_stimulus_z);

for idx=idx_stimulus
    
    rs_comps = squeeze(allComps(:,1,idx,:));
    pv_comps = squeeze(allComps(:,2,idx,:));
    %at_comps = squeeze(allComps(:,3,idx,:));
    
    m_rs_comps = mean(rs_comps(:,:),1);
    m_pv_comps = mean(pv_comps(:,:),1);
    %m_at_comps = mean(at_comps(:,:),1);
    
    m_rs_comps_vol = zeros(MNI_size);
    m_pv_comps_vol = zeros(MNI_size);
    %m_at_comps_vol = zeros(MNI_size);
    
    for iClu=1:nTotalClusters
        
        idx_voxels = find(func_clu == iClu);
        
        nVoxels = length(idx_voxels);
    
        for iVoxel=1:nVoxels

            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

            m_rs_comps_vol(idxx,idxy,idxz) = m_rs_comps(iClu);
            m_pv_comps_vol(idxx,idxy,idxz) = m_pv_comps(iClu);
            %m_at_comps_vol(idxx,idxy,idxz) = m_at_comps(iClu);

        end

    end
    
    nifti_file = func_vol;
    offset = func_vol.dat.offset;
    scl_slope = func_vol.dat.scl_slope;
    scl_inter = func_vol.dat.scl_inter;
    dtype = 'FLOAT32';
    offset = 0;
    dim = func_vol.dat.dim;

%     descrip = 'Stimulus-Triggered';
%     fname = strcat('Ignition-v3-',level,'-TriggeredBy-',int2str(idx),'-onStimContrast-Att','.nii');
%     input_data = m_at_comps_vol; 
%     real_save_image;
    
    descrip = 'Stimulus-Triggered';
    fname = strcat('Ignition-v3-',level,'-TriggeredBy-',int2str(idx),'-onStimContrast-Pass','.nii');
    input_data = m_pv_comps_vol; 
    real_save_image;
    
    descrip = 'Stimulus-Triggered';
    fname = strcat('Ignition-v3-',level,'-TriggeredBy-',int2str(idx),'-onStimContrast-Rest','.nii');
    input_data = m_rs_comps_vol; 
    real_save_image;
    
end
    
end

function doPostTriggeredContrast(level,Condition_Label,pcriterion,freq)

disp('Loading Ignition...');
load(strcat('Ignition-v3-',level,'-All-0.04-0.07-IntegStim2Comps-allComps.mat'));
load('Z:\_PAPERS\ignition\Events\analysis\Ignition-v1-cluster-Contrast.mat');

nTotalClusters = 758;
% pcriterion = 0.01;

N = nTotalClusters;

idx_RestingState = 1;
idx_PassiveViewing = 2;
idx_Attention = 3;

p_attention = zeros(N,N);
h_attention = zeros(N,N);
p_stimulus = zeros(N,N);
h_stimulus = zeros(N,N);

idx_attention = find(contrast_attention_z);
idx_stimulus = find(contrast_stimulus_z);

for idx=idx_attention
    
    for iIg=1:N

       disp(strcat('Ig:',int2str(iIg)));

       p_attention(idx,iIg) = ranksum(squeeze(allComps(:,idx_Attention,idx,iIg)),squeeze(allComps(:,idx_PassiveViewing,idx,iIg)));

       if median(squeeze(allComps(:,idx_Attention,idx,iIg))) > median(squeeze(allComps(:,idx_PassiveViewing,idx,iIg))); h_attention(idx,iIg) = -1; end
       if median(squeeze(allComps(:,idx_Attention,idx,iIg))) < median(squeeze(allComps(:,idx_PassiveViewing,idx,iIg))); h_attention(idx,iIg) = 1; end

    end

end

for idx=idx_attention
    
    trig(idx).contrast_attention = zeros(1,N);

    trig(idx).contrast_attention_z = zeros(1,N);

    for iIg=1:N

        trig(idx).contrast_attention(iIg) = p_attention(idx,iIg);

    end

    trig(idx).contrast_attention(trig(idx).contrast_attention>pcriterion) = 0;

    for iIg=1:N

        if trig(idx).contrast_attention(iIg) ~= 0; trig(idx).contrast_attention_z(iIg) = norminv(trig(idx).contrast_attention(iIg)); end
 
        trig(idx).contrast_attention_z(iIg) = h_attention(idx,iIg) * trig(idx).contrast_attention_z(iIg);
  
    end

end

trig_att = trig;

clear trig 

for idx=idx_stimulus
    
    for iIg=1:N

       disp(strcat('Ig:',int2str(iIg)));

       p_stimulus(idx,iIg) = ranksum(squeeze(allComps(:,idx_PassiveViewing,idx,iIg)),squeeze(allComps(:,idx_RestingState,idx,iIg)));

       if median(squeeze(allComps(:,idx_PassiveViewing,idx,iIg))) > median(squeeze(allComps(:,idx_RestingState,idx,iIg))); h_stimulus(idx,iIg) = -1; end
       if median(squeeze(allComps(:,idx_PassiveViewing,idx,iIg))) < median(squeeze(allComps(:,idx_RestingState,idx,iIg))); h_stimulus(idx,iIg) = 1; end

    end

end

for idx=idx_stimulus
    
    trig(idx).contrast_stimulus = zeros(1,N);

    trig(idx).contrast_stimulus_z = zeros(1,N);

    for iIg=1:N

        trig(idx).contrast_stimulus(iIg) = p_stimulus(idx,iIg);

    end

    trig(idx).contrast_stimulus(trig(idx).contrast_stimulus>pcriterion) = 0;

    for iIg=1:N

        if trig(idx).contrast_stimulus(iIg) ~= 0; trig(idx).contrast_stimulus_z(iIg) = norminv(trig(idx).contrast_stimulus(iIg)); end
 
        trig(idx).contrast_stimulus_z(iIg) = h_stimulus(idx,iIg) * trig(idx).contrast_stimulus_z(iIg);
  
    end

end

trig_stim = trig;

clear trig

save(strcat('Ignition-v3','-',level,'-','PostTrigContrast','-',Condition_Label,'-',num2str(freq(1)),'-',num2str(freq(2)),'.mat'),'trig_stim','trig_att');

end

function getPostTriggeredContrast3DVolume(level,pcriterion)

load(strcat('Ignition-v3-',level,'-PostTrigContrast-All-0.04-0.07.mat'));

pcriterion_label = round(pcriterion * 100);

nTR = 150;
nTotalClusters = 758;
nSubjects = 8;
nROI = 90;
nRuns = 4;
nTotalRuns = 32;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking
TR = 2; %sampling interval(TR)  
nTotalVoxels = 160990;
MNI_size = [91 109 91];
nFlags = 4;

contrast_idx_Trig_Att = [];
contrast_idx_Trig_Stim = [];

iidx = 0;
for idx=1:length(trig_att)
    if ~isempty(find(trig_att(idx).contrast_attention_z))
        iidx = iidx + 1;
        contrast_idx_Trig_Att(iidx) = idx;
    end
end

iidx = 0;
for idx=1:length(trig_stim)
    if ~isempty(find(trig_stim(idx).contrast_stimulus_z))
        iidx = iidx + 1;
        contrast_idx_Trig_Stim(iidx) = idx;
    end
end

func_vol = nifti('LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
func_vol.dat.fname = 'Z:\_DATA\Parcellation\758-Cluster\LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii';
func_clu = func_vol.dat(:,:,:);

z_threshold = abs ( -sqrt(2) * erfcinv(pcriterion*2) );
MNI_size = [91 109 91];

%% Attention

for idx_Trig_Att=contrast_idx_Trig_Att

    contrast_attention_z = trig_att(idx_Trig_Att).contrast_attention_z;

    idx_increase_att = find(contrast_attention_z>z_threshold);
    idx_decrease_att = find(contrast_attention_z<((-1)*z_threshold));

    idx_att_all = [idx_increase_att(:);idx_decrease_att(:)];

    att_increase_vol = zeros(size(func_clu));
    att_decrease_vol = zeros(size(func_clu));

    att_all_vol = zeros(size(func_clu));
    att_all_vol_idx = zeros(size(func_clu));

    for idx=idx_increase_att

        idx_voxels = find(func_clu==idx);

        nVoxels = length(idx_voxels);

        for iVoxel=1:nVoxels

            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

            att_increase_vol(idxx,idxy,idxz) = contrast_attention_z(idx);

            att_all_vol_idx(idxx,idxy,idxz) = idx;

        end

    end

    for idx=idx_decrease_att

        idx_voxels = find(func_clu==idx);

        nVoxels = length(idx_voxels);

        for iVoxel=1:nVoxels

            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

            att_decrease_vol(idxx,idxy,idxz) = contrast_attention_z(idx);

            att_all_vol_idx(idxx,idxy,idxz) = idx;

        end

    end

    for idx=idx_att_all'

        idx_voxels = find(func_clu==idx);

        nVoxels = length(idx_voxels);

        for iVoxel=1:nVoxels

            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

            att_all_vol(idxx,idxy,idxz) = contrast_attention_z(idx);

        end

    end

    nifti_file = func_vol;
    offset = func_vol.dat.offset;
    scl_slope = func_vol.dat.scl_slope;
    scl_inter = func_vol.dat.scl_inter;
    dtype = 'FLOAT32';
    offset = 0;
    dim = func_vol.dat.dim;

    descrip = 'Attention-Increase';
    fname = strcat('Ignition-v3-',level,'-PostTrigContrast-',int2str(idx_Trig_Att),'-Attention-Increase-',num2str(pcriterion_label),'.nii');
    input_data = att_increase_vol; 
    real_save_image;

    descrip = 'Attention-Decrease';
    fname = strcat('Ignition-v3-',level,'-PostTrigContrast-',int2str(idx_Trig_Att),'-Attention-Decrease-',num2str(pcriterion_label),'.nii');
    input_data = att_decrease_vol; 
    real_save_image;

    descrip = 'Attention-All';
    fname = strcat('Ignition-v3-',level,'-PostTrigContrast-',int2str(idx_Trig_Att),'-Attention-All-',num2str(pcriterion_label),'.nii');
    input_data = att_all_vol; 
    real_save_image;

    descrip = 'Attention-All-idx';
    fname = strcat('Ignition-v3-',level,'-PostTrigContrast-',int2str(idx_Trig_Att),'-Attention-All-idx-',num2str(pcriterion_label),'.nii');
    input_data = att_all_vol_idx; 
    real_save_image;

end

%% Stimulus

for idx_Trig_Stim=contrast_idx_Trig_Stim
    
    contrast_stimulus_z = trig_stim(idx_Trig_Stim).contrast_stimulus_z;

    idx_increase_stim = find(contrast_stimulus_z>z_threshold);
    idx_decrease_stim = find(contrast_stimulus_z<((-1)*z_threshold));

    idx_stim_all = [idx_increase_stim(:);idx_decrease_stim(:)];

    stim_increase_vol = zeros(size(func_clu));
    stim_decrease_vol = zeros(size(func_clu));

    stim_all_vol = zeros(size(func_clu));
    stim_all_vol_idx = zeros(size(func_clu));

    for idx=idx_increase_stim

        idx_voxels = find(func_clu==idx);

        nVoxels = length(idx_voxels);

        for iVoxel=1:nVoxels

            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

            stim_increase_vol(idxx,idxy,idxz) = contrast_stimulus_z(idx);

            stim_all_vol_idx(idxx,idxy,idxz) = idx;

        end

    end

    for idx=idx_decrease_stim

        idx_voxels = find(func_clu==idx);

        nVoxels = length(idx_voxels);

        for iVoxel=1:nVoxels

            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

            stim_decrease_vol(idxx,idxy,idxz) = contrast_stimulus_z(idx);

            stim_all_vol_idx(idxx,idxy,idxz) = idx;

        end

    end

    for idx=idx_stim_all'

        idx_voxels = find(func_clu==idx);

        nVoxels = length(idx_voxels);

        for iVoxel=1:nVoxels

            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));

            stim_all_vol(idxx,idxy,idxz) = contrast_stimulus_z(idx);

        end

    end

    nifti_file = func_vol;
    offset = func_vol.dat.offset;
    scl_slope = func_vol.dat.scl_slope;
    scl_inter = func_vol.dat.scl_inter;
    dtype = 'FLOAT32';
    offset = 0;
    dim = func_vol.dat.dim;

    descrip = 'Stimulus-Increase';
    fname = strcat('Ignition-v3-',level,'-PostTrigContrast-',int2str(idx_Trig_Stim),'-Stimulus-Increase-',num2str(pcriterion_label),'.nii');
    input_data = stim_increase_vol; 
    real_save_image;

    descrip = 'Stimulus-Decrease';
    fname = strcat('Ignition-v3-',level,'-PostTrigContrast-',int2str(idx_Trig_Stim),'-Stimulus-Decrease-',num2str(pcriterion_label),'.nii');
    input_data = stim_decrease_vol; 
    real_save_image;

    descrip = 'Stimulus-All';
    fname = strcat('Ignition-v3-',level,'-PostTrigContrast-',int2str(idx_Trig_Stim),'-Stimulus-All-',num2str(pcriterion_label),'.nii');
    input_data = stim_all_vol; 
    real_save_image;

    descrip = 'Stimulus-All-idx';
    fname = strcat('Ignition-v3-',level,'-PostTrigContrast-',int2str(idx_Trig_Stim),'-Stimulus-All-idx-',num2str(pcriterion_label),'.nii');
    input_data = stim_all_vol_idx; 
    real_save_image;

end
    
end

function checkSubGraph(level)

load('Ignition-v3-cluster-All-0.04-0.07-biggest_subgraphs.mat');
load('Ignition-v3-cluster-All-0.04-0.07-allOthers.mat');

nTR = 150;
nTotalClusters = 758;
nSubjects = 8;
nROI = 90;
nRuns = 4;
nTotalRuns = 32;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking
TR = 2; %sampling interval(TR)  
nTotalVoxels = 160990;
MNI_size = [91 109 91];

for iCondition=1:nConditions
    
    for iSeed=1:nTotalClusters
        
        all_tmp_subsize = [];
        for iRun=1:nTotalRuns
            
            tmp_subsize = biggest_subgraphs(iRun,iCondition,iSeed).subsize;
            
            all_tmp_subsize = [all_tmp_subsize, tmp_subsize];
            
            tmp_reliability(iRun) = biggest_subgraphs(iRun,iCondition,iSeed).reliability;
            
        end
        
        med_subsize(iCondition,iSeed) = median(all_tmp_subsize);
        m_subsize(iCondition,iSeed) = mean(all_tmp_subsize);
        max_subsize(iCondition,iSeed) = max(all_tmp_subsize);
        std_subsize(iCondition,iSeed) = std(all_tmp_subsize);
        m_reliability(iCondition,iSeed) = mean(tmp_reliability);
        
    end
    
end

for iCondition=1:nConditions
    
    for iSeed=1:nTotalClusters
        
        tmp_integ = [];
        
        for iRun=1:nTotalRuns
            
            tmp_integ(iRun) = mevokedinteg2(iRun,iCondition,iSeed);
            
        end
        
        med_integ(iCondition,iSeed) = median(tmp_integ);
        
    end
    
end

for iCondition=1:nConditions
    
    for iSeed=1:nTotalClusters
        
        mem = [];
        
        for iRun=1:nTotalRuns

            tmp_mem = biggest_subgraphs(iRun,iCondition,iSeed).membership;
            
            mem = [mem;tmp_mem];

        end
        
        allReal(iCondition,iSeed) = getReliability(mem);
    
    end
    
end

networks = {'DAN' 'VAN' 'VIS' 'AUD' 'LAN' 'FPC' 'SMN' 'DMN'};

for iNet=1:length(networks)
    
    network_label = networks{iNet};
    
    netPerMD = getNetPerMD(network_label);
    
    for iCondition=1:nConditions
    
        for iSeed=1:nTotalClusters
            
            disp(strcat(int2str(iNet),':',int2str(iCondition),':',int2str(iSeed)));

            mem = [];

            for iRun=1:nTotalRuns

                tmp_mem = biggest_subgraphs(iRun,iCondition,iSeed).membership;

                mem = [mem;tmp_mem];

            end
            
            mem = mem(:,find(netPerMD));

            allReal_perNet(iNet,iCondition,iSeed) = getReliability(mem);

        end
    
    end
    
end

save('Ignition-v3-cluster-All-0.04-0.07-biggest_subgraphs-CHECK.mat','med_subsize','m_subsize','max_subsize','std_subsize','m_reliability','med_integ','allReal','allReal_perNet');

end

function checkBCTparam(level)

disp(strcat('...starting loading:',datestr(now)));

load(strcat('Ignition-v4-',level,'-All-0.04-0.07-BCT.mat'));

load('Z:\_PAPERS\ignition\Events\analysis\Ignition-v1-cluster-Contrast.mat');

load('Ignition-v4-cluster-All-0.04-0.07-Run.mat');

disp(strcat('...ending loading:',datestr(now)));

idx_att = find(contrast_attention_z);
idx_stim = find(contrast_stimulus_z);

all_idx_contrast = [idx_att,idx_stim];

nTR = 150;
nTotalClusters = 758;
nSubjects = 8;
nROI = 90;
nRuns = 4;
nTotalRuns = 32;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking
TR = 2; %sampling interval(TR)  
nTotalVoxels = 160990;
MNI_size = [91 109 91];
nFlags = 4;

%%% CLUSTERING_COEF_BU

for idx_clu=1:length(all_idx_contrast)

    seed = all_idx_contrast(idx_clu);

    for iCondition=1:nConditions

        iiSample = 0;

        for iRun=1:nTotalRuns

            nEvents = Run(iRun).Cond(iCondition).nevents(seed);

            for iEvent=1:nEvents

                for iFlag=1:nFlags

                    iiSample = iiSample + 1;

                    CLUSTERING_COEF_BU(seed).C(iCondition,iiSample,:) = BCT(iRun,iCondition,seed,iFlag,iEvent).allBCTparam.CLUSTERING_COEF_BU.C(:);

                end

            end

        end

    end

end

%%% MODULARITY_UND

for idx_clu=1:length(all_idx_contrast)

    seed = all_idx_contrast(idx_clu);

    for iCondition=1:nConditions

        iiSample = 0;

        for iRun=1:nTotalRuns

            nEvents = Run(iRun).Cond(iCondition).nevents(seed);

            for iEvent=1:nEvents

                for iFlag=1:nFlags

                    iiSample = iiSample + 1;

                    MODULARITY_UND(seed).Ci(iCondition,iiSample,:) = BCT(iRun,iCondition,seed,iFlag,iEvent).allBCTparam.MODULARITY_UND.Ci(:);
                    MODULARITY_UND(seed).Q(iCondition,iiSample) = BCT(iRun,iCondition,seed,iFlag,iEvent).allBCTparam.MODULARITY_UND.Q;

                end

            end

        end

    end

end

%%% MODULE_DEGREE_ZSCORE

for idx_clu=1:length(all_idx_contrast)

    seed = all_idx_contrast(idx_clu);

    for iCondition=1:nConditions

        iiSample = 0;

        for iRun=1:nTotalRuns

            nEvents = Run(iRun).Cond(iCondition).nevents(seed);

            for iEvent=1:nEvents

                for iFlag=1:nFlags

                    iiSample = iiSample + 1;

                    MODULE_DEGREE_ZSCORE(seed).Z(iCondition,iiSample,:) = BCT(iRun,iCondition,seed,iFlag,iEvent).allBCTparam.MODULE_DEGREE_ZSCORE.Z(:);
                    
                end

            end

        end

    end

end

%%% SUBGRAPH_CENTRALITY

for idx_clu=1:length(all_idx_contrast)

    seed = all_idx_contrast(idx_clu);

    for iCondition=1:nConditions

        iiSample = 0;

        for iRun=1:nTotalRuns

            nEvents = Run(iRun).Cond(iCondition).nevents(seed);

            for iEvent=1:nEvents

                for iFlag=1:nFlags

                    iiSample = iiSample + 1;

                    SUBGRAPH_CENTRALITY(seed).Cs(iCondition,iiSample,:) = BCT(iRun,iCondition,seed,iFlag,iEvent).allBCTparam.SUBGRAPH_CENTRALITY.Cs(:);
                    
                end

            end

        end

    end

end

%%% TRANSITIVITY_BU

for idx_clu=1:length(all_idx_contrast)

    seed = all_idx_contrast(idx_clu);

    for iCondition=1:nConditions

        iiSample = 0;

        for iRun=1:nTotalRuns

            nEvents = Run(iRun).Cond(iCondition).nevents(seed);

            for iEvent=1:nEvents

                for iFlag=1:nFlags

                    iiSample = iiSample + 1;

                    TRANSITIVITY_BU(seed).C_tri(iCondition,iiSample) = BCT(iRun,iCondition,seed,iFlag,iEvent).allBCTparam.TRANSITIVITY_BU.C_tri;
                    
                end

            end

        end

    end

end

for iRun=1:nTotalRuns

    for iCondition=1:nConditions

        for seed=1:nTotalClusters

            nAllEvents(iCondition,iRun,seed) = Run(iRun).Cond(iCondition).nevents(seed);

        end

    end

end

save(strcat('Ignition-v4-',level,'-All-0.04-0.07-BCT-selected.mat'),'CLUSTERING_COEF_BU','MODULARITY_UND','MODULE_DEGREE_ZSCORE','SUBGRAPH_CENTRALITY','TRANSITIVITY_BU','nAllEvents');

end

function checkTRANSITIVITY_BU(level)

label = 'TRANSITIVITY_BU';

load(strcat('Ignition-v4-',level,'-All-0.04-0.07-BCT-selected.mat'));
load('Z:\_PAPERS\ignition\Events\analysis\Ignition-v1-cluster-Contrast.mat');
load('Ignition-v4-cluster-All-0.04-0.07-allPhasesMatricesBCT.mat');

nTotalRuns = 32;
nTotalClusters = 758;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking

nEffectiveTRs = 131;
nFlags = 4;

idx_resting = 1;
idx_passive = 2;
idx_attention = 3;

idx_att = find(contrast_attention_z);
idx_stim = find(contrast_stimulus_z);

all_idx_contrast = [idx_att,idx_stim];

for idx_clu=1:length(all_idx_contrast)

    seed = all_idx_contrast(idx_clu);
   
    trig(idx_clu).seed = seed;

    my_samples = TRANSITIVITY_BU(seed).C_tri;

    my_samples = my_samples';

    nSamples = size(my_samples,1);

    [trig(idx_clu).contrast_attention_z, trig(idx_clu).contrast_stimulus_z, trig(idx_clu).p_attention, trig(idx_clu).p_stimulus, trig(idx_clu).h_attention, trig(idx_clu).h_stimulus] = doAnyMedianContrastPerParam(my_samples,level,label);

    disp('...doing BOOTstrap');
   
    iiSample = 0;
    for iRun=1:nTotalRuns
        for iTR=1:nEffectiveTRs
            iiSample = iiSample + 1;
            y_att_pass(iiSample) = BCT(idx_attention,iRun,iTR).allBCTparam.TRANSITIVITY_BU.C_tri - BCT(idx_passive,iRun,iTR).allBCTparam.TRANSITIVITY_BU.C_tri;
        end
    end

    y_att_pass(isinf(abs(y_att_pass))) = 0;
    y_att_pass(isnan(abs(y_att_pass))) = 0;

    [trig(idx_clu).ci_att_pass,trig(idx_clu).bootstat_att_pass] = bootci(nSamples, @(x) [mean(x) median(x) std(x)], y_att_pass);

    disp('...doing BOOTstrap');

    iiSample = 0;
    for iRun=1:nTotalRuns
        for iTR=1:nEffectiveTRs
            iiSample = iiSample + 1;
            y_pass_rest(iiSample) = BCT(idx_passive,iRun,iTR).allBCTparam.TRANSITIVITY_BU.C_tri - BCT(idx_resting,iRun,iTR).allBCTparam.TRANSITIVITY_BU.C_tri;
        end
    end

    y_pass_rest(isinf(abs(y_pass_rest))) = 0;
    y_pass_rest(isnan(abs(y_pass_rest))) = 0;

    [trig(idx_clu).ci_pass_rest,trig(idx_clu).bootstat_pass_rest] = bootci(nSamples, @(x) [mean(x) median(x) std(x)], y_pass_rest);

   clear y_att_pass y_pass_rest

end

save(strcat('Ignition-v4','-',level,'-','Contrast-',label,'.mat'),'trig');

end

function checkMODULARITY_UND(level)

label = 'MODULARITY_UND';

load(strcat('Ignition-v4-',level,'-All-0.04-0.07-BCT-selected.mat'));
load('Z:\_PAPERS\ignition\Events\analysis\Ignition-v1-cluster-Contrast.mat');
load('Ignition-v4-cluster-All-0.04-0.07-allPhasesMatricesBCT.mat');

nTotalRuns = 32;
nTotalClusters = 758;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking

nEffectiveTRs = 131;
nFlags = 4;

idx_resting = 1;
idx_passive = 2;
idx_attention = 3;

idx_att = find(contrast_attention_z);
idx_stim = find(contrast_stimulus_z);

all_idx_contrast = [idx_att,idx_stim];

for idx_clu=1:length(all_idx_contrast)

    seed = all_idx_contrast(idx_clu);
   
    trig(idx_clu).seed = seed;

    my_samples = MODULARITY_UND(seed).Q;

    my_samples = my_samples';

    nSamples = size(my_samples,1);

    [trig(idx_clu).contrast_attention_z, trig(idx_clu).contrast_stimulus_z, trig(idx_clu).p_attention, trig(idx_clu).p_stimulus, trig(idx_clu).h_attention, trig(idx_clu).h_stimulus] = doAnyMedianContrastPerParam(my_samples,level,label);

    disp('...doing BOOTstrap');
   
    iiSample = 0;
    for iRun=1:nTotalRuns
        for iTR=1:nEffectiveTRs
            iiSample = iiSample + 1;
            y_att_pass(iiSample) = BCT(idx_attention,iRun,iTR).allBCTparam.MODULARITY_UND.Q - BCT(idx_passive,iRun,iTR).allBCTparam.MODULARITY_UND.Q;
        end
    end

    y_att_pass(isinf(abs(y_att_pass))) = 0;
    y_att_pass(isnan(abs(y_att_pass))) = 0;

    [trig(idx_clu).ci_att_pass,trig(idx_clu).bootstat_att_pass] = bootci(nSamples, @(x) [mean(x) median(x) std(x)], y_att_pass);

    disp('...doing BOOTstrap');

    iiSample = 0;
    for iRun=1:nTotalRuns
        for iTR=1:nEffectiveTRs
            iiSample = iiSample + 1;
            y_pass_rest(iiSample) = BCT(idx_passive,iRun,iTR).allBCTparam.MODULARITY_UND.Q - BCT(idx_resting,iRun,iTR).allBCTparam.MODULARITY_UND.Q;
        end
    end

    y_pass_rest(isinf(abs(y_pass_rest))) = 0;
    y_pass_rest(isnan(abs(y_pass_rest))) = 0;

    [trig(idx_clu).ci_pass_rest,trig(idx_clu).bootstat_pass_rest] = bootci(nSamples, @(x) [mean(x) median(x) std(x)], y_pass_rest);

   clear y_att_pass y_pass_rest

end

save(strcat('Ignition-v4','-',level,'-','Contrast-',label,'.mat'),'trig');


end

function checkCLUSTERING_COEF_BU(level)

label = 'CLUSTERING_COEF_BU';

load(strcat('Ignition-v4-',level,'-All-0.04-0.07-BCT-selected.mat'));
load('Z:\_PAPERS\ignition\Events\analysis\Ignition-v1-cluster-Contrast.mat');
load('Ignition-v4-cluster-All-0.04-0.07-allPhasesMatricesBCT.mat');

nTotalRuns = 32;
nTotalClusters = 758;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking

nEffectiveTRs = 131;
nFlags = 4;

idx_resting = 1;
idx_passive = 2;
idx_attention = 3;

idx_att = find(contrast_attention_z);
idx_stim = find(contrast_stimulus_z);

all_idx_contrast = [idx_att,idx_stim];

for idx_clu=1:length(all_idx_contrast)

    seed = all_idx_contrast(idx_clu);
   
    trig(idx_clu).seed = seed;

    my_samples = CLUSTERING_COEF_BU(seed).C;

    nSamples = nAllEvents(:,:,idx_clu);
    nSamples = max(nSamples(:)) * nFlags * nTotalRuns;

    my_samples= permute(my_samples,[2 1 3]);

    [trig(idx_clu).contrast_attention_z, trig(idx_clu).contrast_stimulus_z, trig(idx_clu).p_attention, trig(idx_clu).p_stimulus, trig(idx_clu).h_attention, trig(idx_clu).h_stimulus] = doAnyMedianContrastPerSeed(my_samples,level,label);

    for target=1:nTotalClusters

        disp(strcat(label,':target:',int2str(target),'...doing BOOTstrap'));

        iiSample = 0;
        for iRun=1:nTotalRuns
            for iTR=1:nEffectiveTRs
                iiSample = iiSample + 1;
                y_att_pass(iiSample) = BCT(idx_attention,iRun,iTR).allBCTparam.CLUSTERING_COEF_BU.C(target) - BCT(idx_passive,iRun,iTR).allBCTparam.CLUSTERING_COEF_BU.C(target);
            end
        end

        y_att_pass(isinf(abs(y_att_pass))) = 0;
        y_att_pass(isnan(abs(y_att_pass))) = 0;

        [trig(idx_clu).cluster(target).ci_att_pass,trig(idx_clu).cluster(target).bootstat_att_pass] = bootci(nSamples, @(x) [mean(x) median(x) std(x)], y_att_pass);

        disp(strcat(label,':target:',int2str(target),'...doing BOOTstrap'));

        iiSample = 0;
        for iRun=1:nTotalRuns
            for iTR=1:nEffectiveTRs
                iiSample = iiSample + 1;
                y_pass_rest(iiSample) = BCT(idx_passive,iRun,iTR).allBCTparam.CLUSTERING_COEF_BU.C(target) - BCT(idx_resting,iRun,iTR).allBCTparam.CLUSTERING_COEF_BU.C(target);
            end
        end

        y_pass_rest(isinf(abs(y_pass_rest))) = 0;
        y_pass_rest(isnan(abs(y_pass_rest))) = 0;

        [trig(idx_clu).cluster(target).ci_pass_rest,trig(idx_clu).cluster(target).bootstat_pass_rest] = bootci(nSamples, @(x) [mean(x) median(x) std(x)], y_pass_rest);

       clear y_att_pass y_pass_rest

   end

end

save(strcat('Ignition-v4','-',level,'-','Contrast-',label,'.mat'),'trig');

end

function checkSUBGRAPH_CENTRALITY(level)

label = 'SUBGRAPH_CENTRALITY';

load(strcat('Ignition-v4-',level,'-All-0.04-0.07-BCT-selected.mat'));
load('Z:\_PAPERS\ignition\Events\analysis\Ignition-v1-cluster-Contrast.mat');
load('Ignition-v4-cluster-All-0.04-0.07-allPhasesMatricesBCT.mat');

nTotalRuns = 32;
nTotalClusters = 758;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking

nEffectiveTRs = 131;
nFlags = 4;

idx_resting = 1;
idx_passive = 2;
idx_attention = 3;

idx_att = find(contrast_attention_z);
idx_stim = find(contrast_stimulus_z);

all_idx_contrast = [idx_att,idx_stim];

for idx_clu=1:length(all_idx_contrast)

    seed = all_idx_contrast(idx_clu);
   
    trig(idx_clu).seed = seed;

    my_samples = SUBGRAPH_CENTRALITY(seed).Cs;

    nSamples = nAllEvents(:,:,idx_clu);
    nSamples = max(nSamples(:)) * nFlags * nTotalRuns;

    my_samples = permute(my_samples,[2 1 3]);

    [trig(idx_clu).contrast_attention_z, trig(idx_clu).contrast_stimulus_z, trig(idx_clu).p_attention, trig(idx_clu).p_stimulus, trig(idx_clu).h_attention, trig(idx_clu).h_stimulus] = doAnyMedianContrastPerSeed(my_samples,level,label);

    for target=1:nTotalClusters

        disp(strcat(label,':target:',int2str(target),'...doing BOOTstrap'));

        iiSample = 0;
        for iRun=1:nTotalRuns
            for iTR=1:nEffectiveTRs
                iiSample = iiSample + 1;
                y_att_pass(iiSample) = BCT(idx_attention,iRun,iTR).allBCTparam.SUBGRAPH_CENTRALITY.Cs(target) - BCT(idx_passive,iRun,iTR).allBCTparam.SUBGRAPH_CENTRALITY.Cs(target);
            end
        end

        y_att_pass(isinf(abs(y_att_pass))) = 0;
        y_att_pass(isnan(abs(y_att_pass))) = 0;

        [trig(idx_clu).cluster(target).ci_att_pass,trig(idx_clu).cluster(target).bootstat_att_pass] = bootci(nSamples, @(x) [mean(x) median(x) std(x)], y_att_pass);

        disp(strcat(label,':target:',int2str(target),'...doing BOOTstrap'));

        iiSample = 0;
        for iRun=1:nTotalRuns
            for iTR=1:nEffectiveTRs
                iiSample = iiSample + 1;
                y_pass_rest(iiSample) = BCT(idx_passive,iRun,iTR).allBCTparam.SUBGRAPH_CENTRALITY.Cs(target) - BCT(idx_resting,iRun,iTR).allBCTparam.SUBGRAPH_CENTRALITY.Cs(target);
            end
        end

        y_pass_rest(isinf(abs(y_pass_rest))) = 0;
        y_pass_rest(isnan(abs(y_pass_rest))) = 0;

        [trig(idx_clu).cluster(target).ci_pass_rest,trig(idx_clu).cluster(target).bootstat_pass_rest] = bootci(nSamples, @(x) [mean(x) median(x) std(x)], y_pass_rest);

       clear y_att_pass y_pass_rest

   end

end

save(strcat('Ignition-v4','-',level,'-','Contrast-',label,'.mat'),'trig');

end

function checkMODULE_DEGREE_ZSCORE(level)

label = 'MODULE_DEGREE_ZSCORE';

load(strcat('Ignition-v4-',level,'-All-0.04-0.07-BCT-selected.mat'));
load('Z:\_PAPERS\ignition\Events\analysis\Ignition-v1-cluster-Contrast.mat');
load('Ignition-v4-cluster-All-0.04-0.07-allPhasesMatricesBCT.mat');

nTotalRuns = 32;
nTotalClusters = 758;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking

nEffectiveTRs = 131;
nFlags = 4;

idx_resting = 1;
idx_passive = 2;
idx_attention = 3;

idx_att = find(contrast_attention_z);
idx_stim = find(contrast_stimulus_z);

all_idx_contrast = [idx_att,idx_stim];

for idx_clu=1:length(all_idx_contrast)

    seed = all_idx_contrast(idx_clu);
   
    trig(idx_clu).seed = seed;

    my_samples = MODULE_DEGREE_ZSCORE(seed).Z;

    nSamples = nAllEvents(:,:,idx_clu);
    nSamples = max(nSamples(:)) * nFlags * nTotalRuns;

    my_samples= permute(my_samples,[2 1 3]);

    [trig(idx_clu).contrast_attention_z, trig(idx_clu).contrast_stimulus_z, trig(idx_clu).p_attention, trig(idx_clu).p_stimulus, trig(idx_clu).h_attention, trig(idx_clu).h_stimulus] = doAnyMedianContrastPerSeed(my_samples,level,label);

    for target=1:nTotalClusters

        disp(strcat(label,':target:',int2str(target),'...doing BOOTstrap'));

        iiSample = 0;
        for iRun=1:nTotalRuns
            for iTR=1:nEffectiveTRs
                iiSample = iiSample + 1;
                y_att_pass(iiSample) = BCT(idx_attention,iRun,iTR).allBCTparam.MODULE_DEGREE_ZSCORE.Z(target) - BCT(idx_passive,iRun,iTR).allBCTparam.MODULE_DEGREE_ZSCORE.Z(target);
            end
        end

        y_att_pass(isinf(abs(y_att_pass))) = 0;
        y_att_pass(isnan(abs(y_att_pass))) = 0;

        [trig(idx_clu).cluster(target).ci_att_pass,trig(idx_clu).cluster(target).bootstat_att_pass] = bootci(nSamples, @(x) [mean(x) median(x) std(x)], y_att_pass);

        disp(strcat(label,':target:',int2str(target),'...doing BOOTstrap'));

        iiSample = 0;
        for iRun=1:nTotalRuns
            for iTR=1:nEffectiveTRs
                iiSample = iiSample + 1;
                y_pass_rest(iiSample) = BCT(idx_passive,iRun,iTR).allBCTparam.MODULE_DEGREE_ZSCORE.Z(target) - BCT(idx_resting,iRun,iTR).allBCTparam.MODULE_DEGREE_ZSCORE.Z(target);
            end
        end

        y_pass_rest(isinf(abs(y_pass_rest))) = 0;
        y_pass_rest(isnan(abs(y_pass_rest))) = 0;

        [trig(idx_clu).cluster(target).ci_pass_rest,trig(idx_clu).cluster(target).bootstat_pass_rest] = bootci(nSamples, @(x) [mean(x) median(x) std(x)], y_pass_rest);

       clear y_att_pass y_pass_rest

   end

end

save(strcat('Ignition-v4','-',level,'-','Contrast-',label,'.mat'),'trig');

end

function [contrast_attention_z, contrast_stimulus_z, p_attention, p_stimulus, h_attention, h_stimulus] = doAnyMedianContrastPerSeed(my_samples,level,label)

nTotalClusters = 758;
pcriterion = 0.01;

N = nTotalClusters;

idx_RestingState = 1;
idx_PassiveViewing = 2;
idx_Attention = 3;

p_attention = zeros(1,N);
h_attention = zeros(1,N);
p_stimulus = zeros(1,N);
h_stimulus = zeros(1,N);

for iIg=1:N
    
   disp(strcat('Ig:',int2str(iIg)));
      
   p_attention(iIg) = ranksum(squeeze(my_samples(:,idx_Attention,iIg)),squeeze(my_samples(:,idx_PassiveViewing,iIg)));
   
   if median(squeeze(my_samples(:,idx_Attention,iIg))) > median(squeeze(my_samples(:,idx_PassiveViewing,iIg))); h_attention(iIg) = -1; end
   if median(squeeze(my_samples(:,idx_Attention,iIg))) < median(squeeze(my_samples(:,idx_PassiveViewing,iIg))); h_attention(iIg) = 1; end
  
   p_stimulus(iIg) = ranksum(squeeze(my_samples(:,idx_PassiveViewing,iIg)),squeeze(my_samples(:,idx_RestingState,iIg)));
   
   if median(squeeze(my_samples(:,idx_PassiveViewing,iIg))) > median(squeeze(my_samples(:,idx_RestingState,iIg))); h_stimulus(iIg) = -1; end
   if median(squeeze(my_samples(:,idx_PassiveViewing,iIg))) < median(squeeze(my_samples(:,idx_RestingState,iIg))); h_stimulus(iIg) = 1; end

end

contrast_attention = zeros(1,N);
contrast_stimulus = zeros(1,N);

contrast_attention_z = zeros(1,N);
contrast_stimulus_z = zeros(1,N);

for iIg=1:N
      
    contrast_attention(iIg) = p_attention(iIg);
    contrast_stimulus(iIg) = p_stimulus(iIg);
    
end

contrast_attention(contrast_attention>pcriterion) = 0;
contrast_stimulus(contrast_stimulus>pcriterion) = 0;

for iIg=1:N
     
    if contrast_attention(iIg) ~= 0; contrast_attention_z(iIg) = norminv(contrast_attention(iIg)); end
    if contrast_stimulus(iIg) ~= 0; contrast_stimulus_z(iIg) = norminv(contrast_stimulus(iIg)); end

    contrast_attention_z(iIg) = h_attention(iIg) * contrast_attention_z(iIg);
    contrast_stimulus_z(iIg) = h_stimulus(iIg) * contrast_stimulus_z(iIg);
    
end

% save(strcat('Ignition-v4','-',level,'-','Contrast-',label,'.mat'),'contrast_stimulus_z','contrast_attention_z');

end

function [contrast_attention_z, contrast_stimulus_z, p_attention, p_stimulus, h_attention, h_stimulus] = doAnyMedianContrastPerParam(my_samples,level,label)

disp(strcat('isINF:',num2str(length(find(isinf(abs(my_samples)))) / length(my_samples(:)))));
my_samples(find(isinf(abs(my_samples)))) = 0;

nTotalClusters = 758;
pcriterion = 0.01;

N = nTotalClusters;

idx_RestingState = 1;
idx_PassiveViewing = 2;
idx_Attention = 3;

h_attention = 0;
h_stimulus = 0;

p_attention = ranksum(squeeze(my_samples(:,idx_Attention)),squeeze(my_samples(:,idx_PassiveViewing)));

if median(squeeze(my_samples(:,idx_Attention))) > median(squeeze(my_samples(:,idx_PassiveViewing))); h_attention = -1; end
if median(squeeze(my_samples(:,idx_Attention))) < median(squeeze(my_samples(:,idx_PassiveViewing))); h_attention = 1; end

p_stimulus = ranksum(squeeze(my_samples(:,idx_PassiveViewing)),squeeze(my_samples(:,idx_RestingState)));

if median(squeeze(my_samples(:,idx_PassiveViewing))) > median(squeeze(my_samples(:,idx_RestingState))); h_stimulus = -1; end
if median(squeeze(my_samples(:,idx_PassiveViewing))) < median(squeeze(my_samples(:,idx_RestingState))); h_stimulus = 1; end

contrast_attention = p_attention;
contrast_stimulus = p_stimulus;

if contrast_attention ~= 0; contrast_attention_z = norminv(contrast_attention); end
if contrast_stimulus ~= 0; contrast_stimulus_z = norminv(contrast_stimulus); end

contrast_attention_z = h_attention * contrast_attention_z;
contrast_stimulus_z = h_stimulus * contrast_stimulus_z;

% save(strcat('Ignition-v4','-',level,'-','Contrast-',label,'.mat'),'contrast_stimulus_z','contrast_attention_z');

end

function doBinomialStatOnFreqMembershipContrastPerSeed

level = 'cluster';

disp('...loading biggest subgraph');
load(strcat('Ignition-v5-',level,'-All-0.04-0.07-biggest_subgraphs.mat'));

disp('...loading contrast');
load('Z:\_PAPERS\ignition\Events\analysis\Ignition-v1-cluster-Contrast.mat');

nTotalRuns = 32;
nConditions = 3;
nTotalClusters = 758;
max_Flag = 1;

pcriterion = 0.01;

idx_rs = 1;
idx_pv = 2;
idx_at = 3;

idx_att = find(contrast_attention_z);
idx_stim = find(contrast_stimulus_z);

all_idx_contrast = [idx_att,idx_stim];

%%% for idx_clu=1:length(all_idx_contrast)
    %%% thisSeed = all_idx_contrast(idx_clu);
for idx_clu=1:nTotalClusters    
    thisSeed = idx_clu;

    total_membership = zeros(nConditions,nTotalClusters);
    for iCondition=1:nConditions

        iiEvent(iCondition) = 0;
        for iRun=1:nTotalRuns
            
            nEvents(iCondition) = size(biggest_subgraphs(iRun,iCondition,thisSeed).perEvent,1);
            
            for iEvent=1:nEvents(iCondition)
                
                iiEvent(iCondition) = iiEvent(iCondition) + 1;

                membership = biggest_subgraphs(iRun,iCondition,thisSeed).perEvent{iEvent,max_Flag}.membership;
                
                total_membership(iCondition,:) = total_membership(iCondition,:) + membership;
                
            end
    
        end

    end
    
    Contrast(idx_clu).thisSeed = thisSeed;
    Contrast(idx_clu).total_membership = total_membership;
    
    for iCondition=1:nConditions
        Contrast(idx_clu).total_membership_norm(iCondition,:) = total_membership(iCondition,:) ./ iiEvent(iCondition);
    end
    
    for iCluster=1:nTotalClusters
        
        NA = iiEvent(idx_rs);
        NB = iiEvent(idx_pv);
        
        OnA = total_membership(idx_rs,iCluster);
        OnB = total_membership(idx_pv,iCluster);
        
        [Contrast(idx_clu).stim_pvalue(iCluster),Contrast(idx_clu).stim_rs_pvalue(iCluster),Contrast(idx_clu).stim_pv_pvalue(iCluster)] = Binostats( NA, OnA, NB, OnB );
        
        NA = iiEvent(idx_pv);
        NB = iiEvent(idx_at);
        
        OnA = total_membership(idx_pv,iCluster);
        OnB = total_membership(idx_at,iCluster);
        
        [Contrast(idx_clu).att_pvalue(iCluster),Contrast(idx_clu).att_pv_pvalue(iCluster),Contrast(idx_clu).att_at_pvalue(iCluster)] = Binostats( NA, OnA, NB, OnB );
        
    end

end

for iTrig=1:length(Contrast)
   
    cluster = Contrast(iTrig).thisSeed;
    
    contrast_attention_z = Contrast(iTrig).att_pvalue;
    contrast_stimulus_z = Contrast(iTrig).stim_pvalue;
    
    idx_attention = find(contrast_attention_z<pcriterion);
    idx_stimulus = find(contrast_stimulus_z<pcriterion);
    
    disp(strcat('Cluster:',int2str(cluster),':','Attention:',num2str(length(idx_attention)),':','Stimulus:',num2str(length(idx_stimulus))));
    
end

save(strcat('Ignition-v7-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast.mat'),'Contrast');

end

function plotBinomialStatOnFreqMembershipContrastPerSeed

level = 'cluster';
load(strcat('Ignition-v6-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast.mat'));

func_vol = nifti('LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
func_vol.dat.fname = 'Z:\_DATA\Parcellation\758-Cluster\LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii';
func_img = func_vol.dat(:,:,:);

zthreshold = 2.3;
MNI_size = [91 109 91];
nTotalClusters = 758;

nTrig = length(Contrast);

for iTrig=1:nTrig
    
    thisSeed = Contrast(iTrig).thisSeed;
    
    stim_rs_pvalue = Contrast(iTrig).stim_rs_pvalue;
    stim_pv_pvalue = Contrast(iTrig).stim_pv_pvalue;
    
    stim_pv_sig = stim_pv_pvalue < stim_rs_pvalue;
    stim_rs_sig = stim_pv_pvalue > stim_rs_pvalue;
    
    stim_pv_z = zeros(size(stim_pv_sig));
    stim_rs_z = zeros(size(stim_rs_sig));
    
    stim_pv_z(stim_pv_sig) = abs(norminv(stim_pv_pvalue(stim_pv_sig))) * 1;
    stim_rs_z(stim_rs_sig) = abs(norminv(stim_rs_pvalue(stim_rs_sig))) * 1;
    
    att_at_pvalue = Contrast(iTrig).att_at_pvalue;
    att_pv_pvalue = Contrast(iTrig).att_pv_pvalue;
    
    att_at_sig = att_at_pvalue < att_pv_pvalue;
    att_pv_sig = att_at_pvalue > att_pv_pvalue;
    
    att_at_z = zeros(size(att_at_sig));
    att_pv_z = zeros(size(att_pv_sig));
    
    att_at_z(att_at_sig) = abs(norminv(att_at_pvalue(att_at_sig))) * 1;
    att_pv_z(att_pv_sig) = abs(norminv(att_pv_pvalue(att_pv_sig))) * 1;
    
    stim_pv_z(abs(stim_pv_z)<zthreshold) = 0;
    stim_rs_z(abs(stim_rs_z)<zthreshold) = 0;

    att_at_z(abs(att_at_z)<zthreshold) = 0;
    att_pv_z(abs(att_pv_z)<zthreshold) = 0;
    
    nifti_file = func_vol;
    offset = func_vol.dat.offset;
    scl_slope = func_vol.dat.scl_slope;
    scl_inter = func_vol.dat.scl_inter;
    dtype = 'FLOAT32';
    offset = 0;
    dim = func_vol.dat.dim;
    
    att_at_z_vol = zeros(MNI_size);
    att_pv_z_vol = zeros(MNI_size);
    stim_pv_z_vol = zeros(MNI_size);
    stim_rs_z_vol = zeros(MNI_size);
    
    for iCluster=1:nTotalClusters
        
        idx_voxels = find(func_img==iCluster);
        
        nVoxels = length(idx_voxels);
        
        for iVoxel=1:nVoxels
            
            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
            
            att_at_z_vol(idxx,idxy,idxz) = att_at_z(iCluster);
            att_pv_z_vol(idxx,idxy,idxz) = att_pv_z(iCluster);
            stim_pv_z_vol(idxx,idxy,idxz) = stim_pv_z(iCluster);
            stim_rs_z_vol(idxx,idxy,idxz) = stim_rs_z(iCluster);
            
        end
        
    end

    descrip = 'Attention-AT-Increase';
    fname = strcat('Ignition-v5-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast-',int2str(thisSeed),'-ATT-AT','.nii');
    input_data = att_at_z_vol; 
    real_save_image;
    
    descrip = 'Attention-PV-Decrease';
    fname = strcat('Ignition-v5-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast-',int2str(thisSeed),'-ATT-PV','.nii');
    input_data = att_pv_z_vol; 
    real_save_image;
    
    descrip = 'Stimulus-PV-Increase';
    fname = strcat('Ignition-v5-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast-',int2str(thisSeed),'-STIM-PV','.nii');
    input_data = stim_pv_z_vol; 
    real_save_image;
    
    descrip = 'Stimulus-RS-Decrease';
    fname = strcat('Ignition-v5-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast-',int2str(thisSeed),'-STIM-RS','.nii');
    input_data = stim_rs_z_vol; 
    real_save_image;

end

end

function plotBinomialStatOnFreqMembershipContrastPerSeedGreenRedBlueGray

level = 'cluster';
load(strcat('Ignition-v7-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast.mat'));

func_vol = nifti('LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
func_vol.dat.fname = 'Z:\_DATA\Parcellation\758-Cluster\LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii';
func_img = func_vol.dat(:,:,:);

zthreshold = 2.3;
MNI_size = [91 109 91];
nTotalClusters = 758;

nTrig = length(Contrast);

idx_rs = 1;
idx_pv = 2;
idx_at = 3;

idx_green = 1;
idx_red = 2;
idx_blue = 3;
idx_grey = 0;

f = 0.05;

for iTrig=1:nTrig
    
    thisSeed = Contrast(iTrig).thisSeed;
    
    total_membership = Contrast(iTrig).total_membership_norm;
    
    m_rs_mem = mean(squeeze(total_membership(idx_rs,:)));
    m_pv_mem = mean(squeeze(total_membership(idx_pv,:)));
    m_at_mem = mean(squeeze(total_membership(idx_at,:)));
    
    stim_rs_pvalue = Contrast(iTrig).stim_rs_pvalue;
    stim_pv_pvalue = Contrast(iTrig).stim_pv_pvalue;
    
    stim_pv_sig = stim_pv_pvalue < stim_rs_pvalue;
    stim_rs_sig = stim_pv_pvalue > stim_rs_pvalue;
    
    stim_pv_z = zeros(size(stim_pv_sig));
    stim_rs_z = zeros(size(stim_rs_sig));
    
    stim_pv_z(stim_pv_sig) = abs(norminv(stim_pv_pvalue(stim_pv_sig))) * 1;
    stim_rs_z(stim_rs_sig) = abs(norminv(stim_rs_pvalue(stim_rs_sig))) * 1;
    
    att_at_pvalue = Contrast(iTrig).att_at_pvalue;
    att_pv_pvalue = Contrast(iTrig).att_pv_pvalue;
    
    att_at_sig = att_at_pvalue < att_pv_pvalue;
    att_pv_sig = att_at_pvalue > att_pv_pvalue;
    
    att_at_z = zeros(size(att_at_sig));
    att_pv_z = zeros(size(att_pv_sig));
    
    att_at_z(att_at_sig) = abs(norminv(att_at_pvalue(att_at_sig))) * 1;
    att_pv_z(att_pv_sig) = abs(norminv(att_pv_pvalue(att_pv_sig))) * 1;
    
% % %     stim_pv_z(abs(stim_pv_z)<zthreshold) = 0;
% % %     stim_rs_z(abs(stim_rs_z)<zthreshold) = 0;
% % % 
% % %     att_at_z(abs(att_at_z)<zthreshold) = 0;
% % %     att_pv_z(abs(att_pv_z)<zthreshold) = 0;
    
    nifti_file = func_vol;
    offset = func_vol.dat.offset;
    scl_slope = func_vol.dat.scl_slope;
    scl_inter = func_vol.dat.scl_inter;
    dtype = 'FLOAT32';
    offset = 0;
    dim = func_vol.dat.dim;
    
    att_at_z_vol = zeros(MNI_size);
    att_pv_z_vol = zeros(MNI_size);
    stim_pv_z_vol = zeros(MNI_size);
    stim_rs_z_vol = zeros(MNI_size);
    
    att_at = zeros(1,nTotalClusters);
    att_pv = zeros(1,nTotalClusters);
    stim_pv = zeros(1,nTotalClusters);
    stim_rs = zeros(1,nTotalClusters);
    
    for iCluster=1:nTotalClusters
        
        idx_voxels = find(func_img==iCluster);
        
        nVoxels = length(idx_voxels);
        
        for iVoxel=1:nVoxels
            
            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
            
            if att_at_z(iCluster) < zthreshold && total_membership(idx_at,iCluster) >= m_at_mem
                
                att_at_z_vol(idxx,idxy,idxz) = idx_green;
                
            elseif att_at_z(iCluster) < zthreshold && total_membership(idx_at,iCluster) < m_at_mem
                
                att_at_z_vol(idxx,idxy,idxz) = idx_grey;
                
            elseif att_at_z(iCluster) >= zthreshold && total_membership(idx_at,iCluster) < m_at_mem
                
                att_at_z_vol(idxx,idxy,idxz) = idx_blue;
                
                att_at(iCluster) = idx_blue;
                
            elseif att_at_z(iCluster) >= zthreshold && total_membership(idx_at,iCluster) >= m_at_mem
                
                att_at_z_vol(idxx,idxy,idxz) = idx_red;
                
                att_at(iCluster) = idx_red;
                
            end
            
            
            if att_pv_z(iCluster) < zthreshold && total_membership(idx_pv,iCluster) > m_pv_mem
                
                att_pv_z_vol(idxx,idxy,idxz) = idx_green;
                
            elseif att_pv_z(iCluster) < zthreshold && total_membership(idx_pv,iCluster) < m_pv_mem
                
                att_pv_z_vol(idxx,idxy,idxz) = idx_grey;
                
            elseif att_pv_z(iCluster) > zthreshold && total_membership(idx_pv,iCluster) < m_pv_mem
                
                att_pv_z_vol(idxx,idxy,idxz) = idx_blue;
                
                att_pv(iCluster) = idx_blue;
                
            elseif att_pv_z(iCluster) > zthreshold && total_membership(idx_pv,iCluster) > m_pv_mem
                
                att_pv_z_vol(idxx,idxy,idxz) = idx_red;
                
                att_pv(iCluster) = idx_red;
                
            end
            
            
            if stim_pv_z(iCluster) < zthreshold && total_membership(idx_pv,iCluster) > m_pv_mem
                
                stim_pv_z_vol(idxx,idxy,idxz) = idx_green;
                
            elseif stim_pv_z(iCluster) < zthreshold && total_membership(idx_pv,iCluster) < m_pv_mem
                
                stim_pv_z_vol(idxx,idxy,idxz) = idx_grey;
                
            elseif stim_pv_z(iCluster) > zthreshold && total_membership(idx_pv,iCluster) < m_pv_mem
                
                stim_pv_z_vol(idxx,idxy,idxz) = idx_blue;
                
                stim_pv(iCluster) = idx_blue;
                
            elseif stim_pv_z(iCluster) > zthreshold && total_membership(idx_pv,iCluster) > m_pv_mem
                
                stim_pv_z_vol(idxx,idxy,idxz) = idx_red;
                
                stim_pv(iCluster) = idx_red;
                
            end
            
            
            if stim_rs_z(iCluster) < zthreshold && total_membership(idx_rs,iCluster) > m_rs_mem
                
                stim_rs_z_vol(idxx,idxy,idxz) = idx_green;
                
            elseif stim_rs_z(iCluster) < zthreshold && total_membership(idx_rs,iCluster) < m_rs_mem
                
                stim_rs_z_vol(idxx,idxy,idxz) = idx_grey;
                
            elseif stim_rs_z(iCluster) > zthreshold && total_membership(idx_rs,iCluster) < m_rs_mem
                
                stim_rs_z_vol(idxx,idxy,idxz) = idx_blue;
                
                stim_rs(iCluster) = idx_blue;
                
            elseif stim_rs_z(iCluster) > zthreshold && total_membership(idx_rs,iCluster) > m_rs_mem
                
                stim_rs_z_vol(idxx,idxy,idxz) = idx_red;
                
                stim_rs(iCluster) = idx_red;
                
            end
          
        end
        
    end
    
    fid = fopen('Membersbhip-v7.txt','a+');
    
    disp(strcat('CLUSTER:',int2str(thisSeed)));
    fprintf(fid,'%s\n',strcat('CLUSTER:',int2str(thisSeed)));
    
    u_labels = getROILabels(find(att_at==idx_blue));
    disp(strcat('att_at-LESS:',strjoin(u_labels,':')));
    fprintf(fid,'%s\n',strcat('att_at-LESS:',strjoin(u_labels,':')));
    
    u_labels = getROILabels(find(att_at==idx_red));
    disp(strcat('att_at-MORE:',strjoin(u_labels,':')));
    fprintf(fid,'%s\n',strcat('att_at-MORE:',strjoin(u_labels,':')));
    
    u_labels = getROILabels(find(att_pv==idx_blue));
    disp(strcat('att_pv-LESS:',strjoin(u_labels,':')));
    fprintf(fid,'%s\n',strcat('att_pv-LESS:',strjoin(u_labels,':')));
    
    u_labels = getROILabels(find(att_pv==idx_red));
    disp(strcat('att_pv-MORE:',strjoin(u_labels,':')));
    fprintf(fid,'%s\n',strcat('att_pv-MORE:',strjoin(u_labels,':')));
    
    u_labels = getROILabels(find(stim_rs==idx_blue));
    disp(strcat('stim_rs-LESS:',strjoin(u_labels,':')));
    fprintf(fid,'%s\n',strcat('stim_rs-LESS:',strjoin(u_labels,':')));
    
    u_labels = getROILabels(find(stim_rs==idx_red));
    disp(strcat('stim_rs-MORE:',strjoin(u_labels,':')));
    fprintf(fid,'%s\n',strcat('stim_rs-MORE:',strjoin(u_labels,':')));
    
    u_labels = getROILabels(find(stim_pv==idx_blue));
    disp(strcat('stim_pv-LESS:',strjoin(u_labels,':')));
    fprintf(fid,'%s\n',strcat('stim_pv-LESS:',strjoin(u_labels,':')));
    
    u_labels = getROILabels(find(stim_pv==idx_red));
    disp(strcat('stim_pv-MORE:',strjoin(u_labels,':')));
    fprintf(fid,'%s\n',strcat('stim_pv-MORE:',strjoin(u_labels,':')));
    
    fclose(fid);

    descrip = 'Attention-AT-Increase';
    fname = strcat('Ignition-v7-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast-',int2str(thisSeed),'-ATT-AT','.nii');
    input_data = att_at_z_vol; 
    real_save_image;
    
    descrip = 'Attention-PV-Decrease';
    fname = strcat('Ignition-v7-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast-',int2str(thisSeed),'-ATT-PV','.nii');
    input_data = att_pv_z_vol; 
    real_save_image;
    
    descrip = 'Stimulus-PV-Increase';
    fname = strcat('Ignition-v7-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast-',int2str(thisSeed),'-STIM-PV','.nii');
    input_data = stim_pv_z_vol; 
    real_save_image;
    
    descrip = 'Stimulus-RS-Decrease';
    fname = strcat('Ignition-v7-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast-',int2str(thisSeed),'-STIM-RS','.nii');
    input_data = stim_rs_z_vol; 
    real_save_image;

end

end

function plotBinomialStatOnFreqMembershipContrast_v2

level = 'cluster';
load(strcat('Ignition-v6-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast.mat'));

func_vol = nifti('LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
func_vol.dat.fname = 'Z:\_DATA\Parcellation\758-Cluster\LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii';
func_img = func_vol.dat(:,:,:);

zthreshold = 2.3;
MNI_size = [91 109 91];
nTotalClusters = 758;

nTrig = length(Contrast);

idx_rs = 1;
idx_pv = 2;
idx_at = 3;

idx_green = 1;
idx_red = 2;
idx_blue = 3;
idx_grey = 0;

f = 0.05;
normal = 1/8;

for iTrig=1:nTrig
    
    thisSeed = Contrast(iTrig).thisSeed;
    
    total_membership = Contrast(iTrig).total_membership_norm;
    
    m_rs_mem = mean(squeeze(total_membership(idx_rs,:)));
    m_pv_mem = mean(squeeze(total_membership(idx_pv,:)));
    m_at_mem = mean(squeeze(total_membership(idx_at,:)));
    
    stim_rs_pvalue = Contrast(iTrig).stim_rs_pvalue;
    stim_pv_pvalue = Contrast(iTrig).stim_pv_pvalue;
    
% % %     stim_pv_sig = stim_pv_pvalue < stim_rs_pvalue;
% % %     stim_rs_sig = stim_pv_pvalue > stim_rs_pvalue;
% % %     
% % %     stim_pv_z = zeros(size(stim_pv_sig));
% % %     stim_rs_z = zeros(size(stim_rs_sig));
    
    stim_pv_z = abs(norminv(stim_pv_pvalue));
    stim_rs_z = abs(norminv(stim_rs_pvalue));
    
    att_at_pvalue = Contrast(iTrig).att_at_pvalue;
    att_pv_pvalue = Contrast(iTrig).att_pv_pvalue;
    
% % %     att_at_sig = att_at_pvalue < att_pv_pvalue;
% % %     att_pv_sig = att_at_pvalue > att_pv_pvalue;
% % %     
% % %     att_at_z = zeros(size(att_at_sig));
% % %     att_pv_z = zeros(size(att_pv_sig));
    
    att_at_z = abs(norminv(att_at_pvalue));
    att_pv_z = abs(norminv(att_pv_pvalue));
    
% % %     stim_pv_z(abs(stim_pv_z)<zthreshold) = 0;
% % %     stim_rs_z(abs(stim_rs_z)<zthreshold) = 0;
% % % 
% % %     att_at_z(abs(att_at_z)<zthreshold) = 0;
% % %     att_pv_z(abs(att_pv_z)<zthreshold) = 0;
    
    nifti_file = func_vol;
    offset = func_vol.dat.offset;
    scl_slope = func_vol.dat.scl_slope;
    scl_inter = func_vol.dat.scl_inter;
    dtype = 'FLOAT32';
    offset = 0;
    dim = func_vol.dat.dim;
    
    att_z_vol = zeros(MNI_size);
    stim_z_vol = zeros(MNI_size);
    
    att = zeros(1,nTotalClusters);
    stim = zeros(1,nTotalClusters);
    
    for iCluster=1:nTotalClusters
        
        idx_voxels = find(func_img==iCluster);
        
        nVoxels = length(idx_voxels);
        
        for iVoxel=1:nVoxels
            
            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
            
            %%% ATTENTION CONTRAST
            
            if att_at_z(iCluster) >= zthreshold && att_at_z(iCluster) > att_pv_z(iCluster) %%% && total_membership(idx_at,iCluster) >= normal
                
                att_z_vol(idxx,idxy,idxz) = idx_red;
                
                att(iCluster) = idx_red;
                
            elseif att_pv_z(iCluster) >= zthreshold && att_at_z(iCluster) < att_pv_z(iCluster) %%% && total_membership(idx_pv,iCluster) >= normal
                
                att_z_vol(idxx,idxy,idxz) = idx_blue;
                
                att(iCluster) = idx_blue;
                
% % %             elseif att_at_z(iCluster) >= zthreshold && att_at_z(iCluster) > att_pv_z(iCluster) && total_membership(idx_at,iCluster) < normal
% % %                 
% % %                 att_z_vol(idxx,idxy,idxz) = idx_grey;
% % %                 
% % %             elseif att_pv_z(iCluster) >= zthreshold && att_at_z(iCluster) < att_pv_z(iCluster) && total_membership(idx_pv,iCluster) < normal
% % %                 
% % %                 att_z_vol(idxx,idxy,idxz) = idx_grey;
                
            elseif att_at_z(iCluster) < zthreshold && att_pv_z(iCluster) < zthreshold && total_membership(idx_at,iCluster) >= normal && total_membership(idx_pv,iCluster) >= normal
                
                att_z_vol(idxx,idxy,idxz) = idx_green;
                
            elseif att_at_z(iCluster) < zthreshold && att_pv_z(iCluster) < zthreshold && total_membership(idx_at,iCluster) < normal && total_membership(idx_pv,iCluster) < normal
                
                att_z_vol(idxx,idxy,idxz) = idx_grey;
                
            end
            
            if stim_pv_z(iCluster) >= zthreshold && stim_pv_z(iCluster) > stim_rs_z(iCluster) %% && total_membership(idx_pv,iCluster) >= normal
                
                stim_z_vol(idxx,idxy,idxz) = idx_red;
                
                stim(iCluster) = idx_red;
                
            elseif stim_rs_z(iCluster) >= zthreshold && stim_pv_z(iCluster) < stim_rs_z(iCluster) %% && total_membership(idx_rs,iCluster) >= normal
                
                stim_z_vol(idxx,idxy,idxz) = idx_blue;
                
                stim(iCluster) = idx_blue;
                
% %             elseif stim_pv_z(iCluster) >= zthreshold && stim_pv_z(iCluster) > stim_rs_z(iCluster) && total_membership(idx_pv,iCluster) < normal
% %                 
% %                 stim_z_vol(idxx,idxy,idxz) = idx_grey;
% %                 
% %             elseif stim_rs_z(iCluster) >= zthreshold && stim_pv_z(iCluster) < stim_rs_z(iCluster) && total_membership(idx_rs,iCluster) < normal
% %                 
% %                 stim_z_vol(idxx,idxy,idxz) = idx_grey;
                
            elseif stim_pv_z(iCluster) < zthreshold && stim_rs_z(iCluster) < zthreshold && total_membership(idx_pv,iCluster) >= normal && total_membership(idx_rs,iCluster) >= normal
                
                stim_z_vol(idxx,idxy,idxz) = idx_green;
                
            elseif stim_pv_z(iCluster) < zthreshold && stim_rs_z(iCluster) < zthreshold && total_membership(idx_pv,iCluster) < normal && total_membership(idx_rs,iCluster) < normal
                
                stim_z_vol(idxx,idxy,idxz) = idx_grey;
                
            end
            
        end
        
    end
    
    disp(strcat('CLUSTER:',int2str(thisSeed)));
    
    u_labels = getROILabels(find(att==idx_blue));
    disp(strcat('att-LESS:',strjoin(u_labels,':')));
    
    u_labels = getROILabels(find(att==idx_red));
    disp(strcat('att-MORE:',strjoin(u_labels,':')));
    
    u_labels = getROILabels(find(stim==idx_blue));
    disp(strcat('stim-LESS:',strjoin(u_labels,':')));
    
    u_labels = getROILabels(find(stim==idx_red));
    disp(strcat('stim-MORE:',strjoin(u_labels,':')));
    
    descrip = 'Attention-AT-Increase';
    fname = strcat('Ignition-v6-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast-',int2str(thisSeed),'-ATT','.nii');
    input_data = att_z_vol; 
    real_save_image;
    
% % %     descrip = 'Attention-PV-Decrease';
% % %     fname = strcat('Ignition-v6-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast-',int2str(thisSeed),'-ATT-PV','.nii');
% % %     input_data = att_pv_z_vol; 
% % %     real_save_image;
    
    descrip = 'Stimulus-PV-Increase';
    fname = strcat('Ignition-v6-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast-',int2str(thisSeed),'-STIM','.nii');
    input_data = stim_z_vol; 
    real_save_image;
    
% % %     descrip = 'Stimulus-RS-Decrease';
% % %     fname = strcat('Ignition-v6-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast-',int2str(thisSeed),'-STIM-RS','.nii');
% % %     input_data = stim_rs_z_vol; 
% % %     real_save_image;

end

end

function doFreqMembershipMedianContrastPerSeed

level = 'cluster';

load('Ignition-v10-cluster-All-0.04-0.07-checkSUBGRAPH_PEREVENT_CORBETTA_cluster_simpe_model_oneNORM.mat');

nTotalClusters = 758;
pcriterion = 0.01;

N = nTotalClusters;
max_flag = 1;

idx_RestingState = 1;
idx_PassiveViewing = 2;
idx_Attention = 3;

p_attention = zeros(1,N);
h_attention = zeros(1,N);
p_stimulus = zeros(1,N);
h_stimulus = zeros(1,N);

for iTrig=1:length(trig)
    
    tmp = trig(iTrig).stat{idx_RestingState+1,44};
    my_samples(:,idx_RestingState,:) = squeeze(tmp(:,max_flag,:));
    
    tmp = trig(iTrig).stat{idx_PassiveViewing+1,44};
    my_samples(:,idx_PassiveViewing,:) = squeeze(tmp(:,max_flag,:));
    
    tmp = trig(iTrig).stat{idx_Attention+1,44};
    my_samples(:,idx_Attention,:) = squeeze(tmp(:,max_flag,:));
    
    trig(iTrig).h_attention = zeros(1,N);
    trig(iTrig).h_stimulus = zeros(1,N);
    
    trig(iTrig).p_attention = zeros(1,N);
    trig(iTrig).p_stimulus = zeros(1,N);
    
    for iIg=1:N

       disp(strcat('Ig:',int2str(iIg)));

       trig(iTrig).p_attention(iIg) = ranksum(squeeze(my_samples(:,idx_Attention,iIg)),squeeze(my_samples(:,idx_PassiveViewing,iIg)));

       if median(squeeze(my_samples(:,idx_Attention,iIg))) > median(squeeze(my_samples(:,idx_PassiveViewing,iIg))); trig(iTrig).h_attention(iIg) = -1; end
       if median(squeeze(my_samples(:,idx_Attention,iIg))) < median(squeeze(my_samples(:,idx_PassiveViewing,iIg))); trig(iTrig).h_attention(iIg) = 1; end

       trig(iTrig).p_stimulus(iIg) = ranksum(squeeze(my_samples(:,idx_PassiveViewing,iIg)),squeeze(my_samples(:,idx_RestingState,iIg)));

       if median(squeeze(my_samples(:,idx_PassiveViewing,iIg))) > median(squeeze(my_samples(:,idx_RestingState,iIg))); trig(iTrig).h_stimulus(iIg) = -1; end
       if median(squeeze(my_samples(:,idx_PassiveViewing,iIg))) < median(squeeze(my_samples(:,idx_RestingState,iIg))); trig(iTrig).h_stimulus(iIg) = 1; end

    end

end

for iTrig=1:length(trig)
    
    trig(iTrig).contrast_attention = zeros(1,N);
    trig(iTrig).contrast_stimulus = zeros(1,N);

    trig(iTrig).contrast_attention_z = zeros(1,N);
    trig(iTrig).contrast_stimulus_z = zeros(1,N);

    for iIg=1:N

        trig(iTrig).contrast_attention(iIg) = trig(iTrig).p_attention(iIg);
        trig(iTrig).contrast_stimulus(iIg) = trig(iTrig).p_stimulus(iIg);

    end

    trig(iTrig).contrast_attention(trig(iTrig).contrast_attention>pcriterion) = 0;
    trig(iTrig).contrast_stimulus(trig(iTrig).contrast_stimulus>pcriterion) = 0;

    for iIg=1:N

        if trig(iTrig).contrast_attention(iIg) ~= 0; trig(iTrig).contrast_attention_z(iIg) = norminv(trig(iTrig).contrast_attention(iIg)); end
        if trig(iTrig).contrast_stimulus(iIg) ~= 0; trig(iTrig).contrast_stimulus_z(iIg) = norminv(trig(iTrig).contrast_stimulus(iIg)); end

        trig(iTrig).contrast_attention_z(iIg) = trig(iTrig).h_attention(iIg) * trig(iTrig).contrast_attention_z(iIg);
        trig(iTrig).contrast_stimulus_z(iIg) = trig(iTrig).h_stimulus(iIg) * trig(iTrig).contrast_stimulus_z(iIg);

    end

end

save(strcat('Ignition-v10','-',level,'-','Contrast-','checkSUBGRAPH_PEREVENT_CORBETTA_cluster_simpe_model_oneNORM','.mat'),'trig');

%%% DISPLAY CLUSTERS ON CONTRAST

for iTrig=1:length(trig)
   
    cluster = trig(iTrig).stat{2,1};
    
    contrast_attention_z = trig(iTrig).contrast_attention_z;
    contrast_stimulus_z = trig(iTrig).contrast_stimulus_z;
    
    idx_attention = find(contrast_attention_z);
    idx_stimulus = find(contrast_stimulus_z);
    
    disp(strcat('Cluster:',int2str(cluster),':','Attention:',num2str(idx_attention),':','Stimulus:',num2str(idx_stimulus)));
    
end

end

function checkSUBGRAPH_PEREVENT(level)

load(strcat('Ignition-v5-',level,'-All-0.04-0.07-biggest_subgraphs.mat'));

load('Z:\_PAPERS\ignition\Events\analysis\Ignition-v1-cluster-Contrast.mat');

nTotalRuns = 32;
nConditions = 3;
nTotalClusters = 758;
nTotalFlags = 4;

idx_resting = 1;
idx_passive = 2;
idx_attention = 3;

idx_att = find(contrast_attention_z);
idx_stim = find(contrast_stimulus_z);

all_idx_contrast = [idx_att,idx_stim];

for idx_clu=1:length(all_idx_contrast)
    
    thisSeed = all_idx_contrast(idx_clu);
    
    for iCondition=1:nConditions
        
        trig(idx_clu).stat{iCondition,1} = thisSeed;
        
        freq_membership = zeros(nTotalRuns,nTotalFlags,nTotalClusters);
    
        for iRun=1:nTotalRuns

            nEvents(iRun) = size(biggest_subgraphs(iRun,iCondition,thisSeed).perEvent,1);
            
            for iFlag=1:nTotalFlags
                
                membership = zeros(1,nTotalClusters);
                real_membership = zeros(nEvents(iRun),nTotalClusters);
                
                for iEvent=1:nEvents(iRun)
                    
                    membership = membership + biggest_subgraphs(iRun,iCondition,thisSeed).perEvent{iEvent,iFlag}.membership;
                    
                    real_membership(iEvent,:) = biggest_subgraphs(iRun,iCondition,thisSeed).perEvent{iEvent,iFlag}.membership;
                    
                end
                
                membership = membership ./ nTotalFlags;
                
                freq_membership(iRun,iFlag,:) = membership;
                
                reliability(iRun,iFlag) = getReliability(real_membership);
                
            end

        end

        trig(idx_clu).stat{iCondition,2} = mean(nEvents);
        trig(idx_clu).stat{iCondition,3} = std(nEvents);
        
        trig(idx_clu).stat{iCondition,4} = squeeze(mean(freq_membership,1));
        
        trig(idx_clu).stat{iCondition,5} = squeeze(mean(reliability,1));
        
        clear nEvents
    
    end
    
end

save(strcat('Ignition-v5-',level,'-All-0.04-0.07-','checkReliability','.mat'),'trig');


end

function checkSUBGRAPH_PEREVENT_CORBETTA_voxel(level)

disp('...loading biggest subgraph');
load(strcat('Ignition-v5-',level,'-All-0.04-0.07-biggest_subgraphs.mat'));

disp('...loading contrast');
load('Z:\_PAPERS\ignition\Events\analysis\Ignition-v1-cluster-Contrast.mat');

load('Z:\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat');

nROIs = 90;

disp('...functional networks');

networks = {'DAN' 'VAN' 'VIS' 'AUD' 'LAN' 'FPC' 'SMN' 'DMN'};
folder_nets = 'Z:\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually\';

DAN = nifti(strcat(folder_nets,'DAN-bin.nii'));
DAN.dat.fname = strcat(folder_nets,'DAN-bin.nii');
DAN_img = DAN.dat(:,:,:);

VAN = nifti(strcat(folder_nets,'VAN-bin.nii'));
VAN.dat.fname = strcat(folder_nets,'VAN-bin.nii');
VAN_img = VAN.dat(:,:,:);

VIS = nifti(strcat(folder_nets,'VIS-bin.nii'));
VIS.dat.fname = strcat(folder_nets,'VIS-bin.nii');
VIS_img = VIS.dat(:,:,:);

AUD = nifti(strcat(folder_nets,'AUD-bin.nii'));
AUD.dat.fname = strcat(folder_nets,'AUD-bin.nii');
AUD_img = AUD.dat(:,:,:);

LAN = nifti(strcat(folder_nets,'LAN-bin.nii'));
LAN.dat.fname = strcat(folder_nets,'LAN-bin.nii');
LAN_img = LAN.dat(:,:,:);

FPC = nifti(strcat(folder_nets,'FPC-bin.nii'));
FPC.dat.fname = strcat(folder_nets,'FPC-bin.nii');
FPC_img = FPC.dat(:,:,:);

SMN = nifti(strcat(folder_nets,'SMN-bin.nii'));
SMN.dat.fname = strcat(folder_nets,'SMN-bin.nii');
SMN_img = SMN.dat(:,:,:);

DMN = nifti(strcat(folder_nets,'DMN-bin.nii'));
DMN.dat.fname = strcat(folder_nets,'DMN-bin.nii');
DMN_img = DMN.dat(:,:,:);

folder_fun = 'Z:\_DATA\Parcellation\758-Cluster\';

MD = nifti(strcat(folder_fun,'LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii'));
MD.dat.fname = strcat(folder_fun,'LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
MD_img = MD.dat(:,:,:);

nTotalRuns = 32;
nTotalVoxels = 160990;
nConditions = 3;
nTotalClusters = 758;
nTotalFlags = 4;

MNI_size = [91 109 91];

idx_resting = 1;
idx_passive = 2;
idx_attention = 3;

idx_att = find(contrast_attention_z);
idx_stim = find(contrast_stimulus_z);

all_idx_contrast = [idx_att,idx_stim];

for idx_clu=1:length(all_idx_contrast)
    
    thisSeed = all_idx_contrast(idx_clu);
    
    for iCondition=1:nConditions
        
        trig(idx_clu).stat{iCondition+1,1} = thisSeed;
        
        freq_membership = zeros(nTotalRuns,nTotalFlags,nTotalClusters);
    
        for iRun=1:nTotalRuns
            
            disp(strcat(int2str(idx_clu),':',int2str(iCondition),':',int2str(iRun)));

            nEvents(iRun) = size(biggest_subgraphs(iRun,iCondition,thisSeed).perEvent,1);
            
            for iFlag=1:nTotalFlags
                
                membership = zeros(1,nTotalClusters);
                real_membership = zeros(nEvents(iRun),nTotalClusters);
                
                for iEvent=1:nEvents(iRun)
                    
                    membership = membership + biggest_subgraphs(iRun,iCondition,thisSeed).perEvent{iEvent,iFlag}.membership;
                    
                    real_membership(iEvent,:) = biggest_subgraphs(iRun,iCondition,thisSeed).perEvent{iEvent,iFlag}.membership;
                    
                end
                
                membership = membership ./ nTotalFlags;
                
                % ROI = getMD758Clusters;
                
                mem_DAN = zeros(1,nTotalVoxels);
                mem_VAN = zeros(1,nTotalVoxels);
                mem_VIS = zeros(1,nTotalVoxels);
                mem_LAN = zeros(1,nTotalVoxels);
                mem_SMN = zeros(1,nTotalVoxels);
                mem_DMN = zeros(1,nTotalVoxels);
                mem_FPC = zeros(1,nTotalVoxels);
                mem_AUD = zeros(1,nTotalVoxels);
                
                real_DAN = zeros(nEvents(iRun),nTotalVoxels);
                real_VAN = zeros(nEvents(iRun),nTotalVoxels);
                real_VIS = zeros(nEvents(iRun),nTotalVoxels);
                real_LAN = zeros(nEvents(iRun),nTotalVoxels);
                real_SMN = zeros(nEvents(iRun),nTotalVoxels);
                real_DMN = zeros(nEvents(iRun),nTotalVoxels);
                real_FPC = zeros(nEvents(iRun),nTotalVoxels);
                real_AUD = zeros(nEvents(iRun),nTotalVoxels);
                
                iiVoxel = 0;
                iiCluster = 0;
                for iROI=1:nROIs
                    
                    nClusters = ROI(iROI).nClusters;
                    
                    for iCluster=1:nClusters
                        
                        iiCluster = iiCluster + 1;
                        
                        idx_voxels = ROI(iROI).clusters(iCluster).idx_voxels;
                        
                        nVoxels = length(idx_voxels);
                        
                        for iVoxel=1:nVoxels
                            
                            iiVoxel = iiVoxel + 1;
                            
                            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
                            
                            if DAN_img(idxx,idxy,idxz)
                                
                                mem_DAN(iiVoxel) = membership(iiCluster);
                                real_DAN(:,iiVoxel) = real_membership(:,iiCluster);
                                
                            end
                            
                            if VAN_img(idxx,idxy,idxz)
                                
                                mem_VAN(iiVoxel) = membership(iiCluster);
                                real_VAN(:,iiVoxel) = real_membership(:,iiCluster);
                                
                            end
                            
                            if VIS_img(idxx,idxy,idxz)
                                
                                mem_VIS(iiVoxel) = membership(iiCluster);
                                real_VIS(:,iiVoxel) = real_membership(:,iiCluster);
                                
                            end
                            
                            if LAN_img(idxx,idxy,idxz)
                                
                                mem_LAN(iiVoxel) = membership(iiCluster);
                                real_LAN(:,iiVoxel) = real_membership(:,iiCluster);
                                
                            end
                            
                            if FPC_img(idxx,idxy,idxz)
                                
                                mem_FPC(iiVoxel) = membership(iiCluster);
                                real_FPC(:,iiVoxel) = real_membership(:,iiCluster);
                                
                            end
                            
                            if DMN_img(idxx,idxy,idxz)
                                
                                mem_DMN(iiVoxel) = membership(iiCluster);
                                real_DMN(:,iiVoxel) = real_membership(:,iiCluster);
                                
                            end
                            
                            if SMN_img(idxx,idxy,idxz)
                                
                                mem_SMN(iiVoxel) = membership(iiCluster);
                                real_SMN(:,iiVoxel) = real_membership(:,iiCluster);
                                
                            end
                            
                            if AUD_img(idxx,idxy,idxz)
                                
                                mem_AUD(iiVoxel) = membership(iiCluster);
                                real_AUD(:,iiVoxel) = real_membership(:,iiCluster);
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                freq_membership_DAN(iRun,iFlag,:) = mem_DAN;
                freq_membership_VAN(iRun,iFlag,:) = mem_VAN;
                freq_membership_VIS(iRun,iFlag,:) = mem_VIS;
                freq_membership_LAN(iRun,iFlag,:) = mem_LAN;
                freq_membership_DMN(iRun,iFlag,:) = mem_DMN;
                freq_membership_FPC(iRun,iFlag,:) = mem_FPC;
                freq_membership_SMN(iRun,iFlag,:) = mem_SMN;
                freq_membership_AUD(iRun,iFlag,:) = mem_AUD;
                
                reliability_DAN(iRun,iFlag) = getReliability(real_DAN);
                reliability_VAN(iRun,iFlag) = getReliability(real_VAN);
                reliability_VIS(iRun,iFlag) = getReliability(real_VIS);
                reliability_LAN(iRun,iFlag) = getReliability(real_LAN);
                reliability_DMN(iRun,iFlag) = getReliability(real_DMN);
                reliability_FPC(iRun,iFlag) = getReliability(real_FPC);
                reliability_SMN(iRun,iFlag) = getReliability(real_SMN);
                reliability_AUD(iRun,iFlag) = getReliability(real_AUD);
                
            end

        end

        trig(idx_clu).stat{iCondition+1,2} = mean(nEvents);
        trig(idx_clu).stat{iCondition+1,3} = std(nEvents);
        
        trig(idx_clu).stat{1,4} = 'mean-freq-DAN';
        trig(idx_clu).stat{1,5} = 'mean-freq-VAN';
        trig(idx_clu).stat{1,6} = 'mean-freq-VIS';
        trig(idx_clu).stat{1,7} = 'mean-freq-LAN';
        trig(idx_clu).stat{1,8} = 'mean-freq-FPC';
        trig(idx_clu).stat{1,9} = 'mean-freq-SMN';
        trig(idx_clu).stat{1,10} = 'mean-freq-AUD';
        trig(idx_clu).stat{1,11} = 'mean-freq-DMN';
        
        trig(idx_clu).stat{iCondition+1,4} = squeeze(mean(freq_membership_DAN,1));
        trig(idx_clu).stat{iCondition+1,5} = squeeze(mean(freq_membership_VAN,1));
        trig(idx_clu).stat{iCondition+1,6} = squeeze(mean(freq_membership_VIS,1));
        trig(idx_clu).stat{iCondition+1,7} = squeeze(mean(freq_membership_LAN,1));
        trig(idx_clu).stat{iCondition+1,8} = squeeze(mean(freq_membership_FPC,1));
        trig(idx_clu).stat{iCondition+1,9} = squeeze(mean(freq_membership_SMN,1));
        trig(idx_clu).stat{iCondition+1,10} = squeeze(mean(freq_membership_AUD,1));
        trig(idx_clu).stat{iCondition+1,11} = squeeze(mean(freq_membership_DMN,1));
        
        trig(idx_clu).stat{1,12} = 'mean-real-DAN';
        trig(idx_clu).stat{1,13} = 'mean-real-VAN';
        trig(idx_clu).stat{1,14} = 'mean-real-VIS';
        trig(idx_clu).stat{1,15} = 'mean-real-LAN';
        trig(idx_clu).stat{1,16} = 'mean-real-FPC';
        trig(idx_clu).stat{1,17} = 'mean-real-SMN';
        trig(idx_clu).stat{1,18} = 'mean-real-AUD';
        trig(idx_clu).stat{1,19} = 'mean-real-DMN';
        
        trig(idx_clu).stat{iCondition+1,12} = squeeze(mean(reliability_DAN,1));
        trig(idx_clu).stat{iCondition+1,13} = squeeze(mean(reliability_VAN,1));
        trig(idx_clu).stat{iCondition+1,14} = squeeze(mean(reliability_VIS,1));
        trig(idx_clu).stat{iCondition+1,15} = squeeze(mean(reliability_LAN,1));
        trig(idx_clu).stat{iCondition+1,16} = squeeze(mean(reliability_FPC,1));
        trig(idx_clu).stat{iCondition+1,17} = squeeze(mean(reliability_SMN,1));
        trig(idx_clu).stat{iCondition+1,18} = squeeze(mean(reliability_AUD,1));
        trig(idx_clu).stat{iCondition+1,19} = squeeze(mean(reliability_DMN,1));
        
        trig(idx_clu).stat{1,20} = 'std-freq-DAN';
        trig(idx_clu).stat{1,21} = 'std-freq-VAN';
        trig(idx_clu).stat{1,22} = 'std-freq-VIS';
        trig(idx_clu).stat{1,23} = 'std-freq-LAN';
        trig(idx_clu).stat{1,24} = 'std-freq-FPC';
        trig(idx_clu).stat{1,25} = 'std-freq-SMN';
        trig(idx_clu).stat{1,26} = 'std-freq-AUD';
        trig(idx_clu).stat{1,27} = 'std-freq-DMN';
        
        trig(idx_clu).stat{iCondition+1,20} = squeeze(std(freq_membership_DAN,0,1));
        trig(idx_clu).stat{iCondition+1,21} = squeeze(std(freq_membership_VAN,0,1));
        trig(idx_clu).stat{iCondition+1,22} = squeeze(std(freq_membership_VIS,0,1));
        trig(idx_clu).stat{iCondition+1,23} = squeeze(std(freq_membership_LAN,0,1));
        trig(idx_clu).stat{iCondition+1,24} = squeeze(std(freq_membership_FPC,0,1));
        trig(idx_clu).stat{iCondition+1,25} = squeeze(std(freq_membership_SMN,0,1));
        trig(idx_clu).stat{iCondition+1,26} = squeeze(std(freq_membership_AUD,0,1));
        trig(idx_clu).stat{iCondition+1,27} = squeeze(std(freq_membership_DMN,0,1));
        
        trig(idx_clu).stat{1,28} = 'std-real-DAN';
        trig(idx_clu).stat{1,29} = 'std-real-VAN';
        trig(idx_clu).stat{1,30} = 'std-real-VIS';
        trig(idx_clu).stat{1,31} = 'std-real-LAN';
        trig(idx_clu).stat{1,32} = 'std-real-FPC';
        trig(idx_clu).stat{1,33} = 'std-real-SMN';
        trig(idx_clu).stat{1,34} = 'std-real-AUD';
        trig(idx_clu).stat{1,35} = 'std-real-DMN';
        
        trig(idx_clu).stat{iCondition+1,28} = squeeze(std(reliability_DAN,0,1));
        trig(idx_clu).stat{iCondition+1,29} = squeeze(std(reliability_VAN,0,1));
        trig(idx_clu).stat{iCondition+1,30} = squeeze(std(reliability_VIS,0,1));
        trig(idx_clu).stat{iCondition+1,31} = squeeze(std(reliability_LAN,0,1));
        trig(idx_clu).stat{iCondition+1,32} = squeeze(std(reliability_FPC,0,1));
        trig(idx_clu).stat{iCondition+1,33} = squeeze(std(reliability_SMN,0,1));
        trig(idx_clu).stat{iCondition+1,34} = squeeze(std(reliability_AUD,0,1));
        trig(idx_clu).stat{iCondition+1,35} = squeeze(std(reliability_DMN,0,1));
        
        clear nEvents
    
    end
    
end

save(strcat('Ignition-v5-',level,'-All-0.04-0.07-','checkSUBGRAPH_PEREVENT_CORBETTA','.mat'));

end

function get3DVolumeReliability(level)

load(strcat('Ignition-v5-',level,'-All-0.04-0.07-','checkReliability','.mat'));

func_vol = nifti('LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
func_vol.dat.fname = 'Z:\_DATA\Parcellation\758-Cluster\LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii';
func_clu = func_vol.dat(:,:,:);

idx_resting = 1;
idx_passive = 2;
idx_attention = 3;

Condition_label = {'RS' 'PV' 'AT'};

MNI_size = [91 109 91];

nConditions = 3;
nTotalFlags = 4;
nTotalClusters = 758;

nSeed = length(trig);

for iSeed=1:nSeed
    
    thisSeed = trig(iSeed).stat{1,1};
    
    for iCondition=1:nConditions
        
        members = trig(iSeed).stat{iCondition,4};
        
        for iFlag=1:nTotalFlags
            
            thisFlagsMembers = squeeze(members(iFlag,:));
            
            volume = zeros(MNI_size);
            
            for iClu=1:nTotalClusters
        
                idx_voxels = find(func_clu == iClu);
        
                nVoxels = length(idx_voxels);
    
                for iVoxel=1:nVoxels

                    [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
                    volume(idxx,idxy,idxz) = thisFlagsMembers(iClu);

                end

            end
    
            nifti_file = func_vol;
            offset = func_vol.dat.offset;
            scl_slope = func_vol.dat.scl_slope;
            scl_inter = func_vol.dat.scl_inter;
            dtype = 'FLOAT32';
            offset = 0;
            dim = func_vol.dat.dim;

            descrip = 'Percentage-Of-Membership';
            fname = strcat('Ignition-v5-',level,'-TriggeredBy-',int2str(thisSeed),'-',Condition_label{iCondition},'-flag-',int2str(iFlag),'.nii');
            input_data = volume; 
            real_save_image;
            
        end
      
    end
    
end


end

function plotReliability(level)

load(strcat('Ignition-v5-',level,'-All-0.04-0.07-','checkReliability','.mat'));

nConditions = 3;
condition_color_A = {'b', 'k', 'r'};
condition_color_B = {'b*', 'k*', 'r*'};
nSeed = length(trig);

min_real = 0.3;
max_real = 0.7;

flags = 1:4;

for iSeed=1:nSeed
    
    thisSeed = trig(iSeed).stat{1,1};
    
    f = figure;
    
    allReal = [];
    for iCondition=1:nConditions
        
        real = trig(iSeed).stat{iCondition,5};
        
        allReal = [allReal,real];
        
        plot(flags,real,condition_color_A{iCondition});
        hold on
        plot(flags,real,condition_color_B{iCondition});
        hold on
        ylim([min_real max_real]);
        xlim([flags(1)-1 flags(end)+1]);
        
    end
    
    print(f,strcat('Ignition-v5-',level,'-All-0.04-0.07-',int2str(thisSeed),'-checkReliability','.epsc'),'-depsc');
    
end

end

function checkEvents

nRuns = 32;

idx_clu = 108;
idx_condition = 3;

eventsMatrix = [];

for iRun=1:nRuns
    
    ee = squeeze(Events(iRun,idx_condition,idx_clu,:));
    
    idx_ee = find(ee);
    
    for iEE=1:length(idx_ee)
    
        eventsMatrix = [eventsMatrix, Run(iRun).Cond(idx_condition).T(idx_ee(iEE)).phasematrix(:,idx_clu)];
            
    end
    
end

m_events = mean(eventsMatrix,2);
m_clusters = mean(eventsMatrix,1);

end

function checkFLAGS

nTotalRuns = 32;
seed = 108;
nTotalFlags = 4;

Flag_Color{2,1} = 'b-';
Flag_Color{2,2} = 'b-';
Flag_Color{2,3} = 'b-';
Flag_Color{2,4} = 'b-';

Flag_Color{3,1} = 'r-';
Flag_Color{3,2} = 'r-';
Flag_Color{3,3} = 'r-';
Flag_Color{3,4} = 'r-';


for iFlag=1:nTotalFlags
    
    f = figure;
    
    for iCondition=2:3

        for iRun=1:nTotalRuns

            nevents2 = Run(iRun).Cond(iCondition).nevents;

            Flag(iFlag,iRun,iCondition) = mean(squeeze(IntegStim2(iRun,iCondition,seed,iFlag,1:nevents2(seed))));

        end

        plot(squeeze(Flag(iFlag,:,iCondition)),Flag_Color{iCondition,iFlag});
        hold on

    end
    
    print(f,strcat('Flag-',int2str(iFlag),'-iCond-',int2str(iCondition),'.eps'),'-depsc');

end

end

function plotSUBGRAPH_PEREVENT_CORBETTA_cluster(level)

load('Z:\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat');

load(strcat('Ignition-v9-',level,'-All-0.04-0.07-','checkSUBGRAPH_PEREVENT_CORBETTA_cluster_simple_model_NORM','.mat'));

load('Z:\MD758-Corbetta-density-simple-model.mat');

MNI_size = [91 109 91];

mem_OUT = double( mem_DAN | mem_VAN | mem_VIS | mem_LAN | mem_FPC | mem_DMN | mem_SMN | mem_AUD );
mem_OUT = double(~mem_OUT);

idx_rs = 1;
idx_pv = 2;
idx_at = 3;

nFlags = 4;
nConditions = 3;
nTotalVoxels = 160990;
nNets = 8;
nTrig = 14;

net_label = {'DAN' 'VAN' 'VIS' 'LAN' 'FPC' 'SMN' 'AUD' 'DMN'};

color(1,:) = [0,0,128]./255;
color(2,:) = [0.5,0,0.9];
color(6,:) = [0,191,255]./255;
color(3,:) = [0,100,0]./255;
color(5,:) = [255,255,0]./255;
color(4,:) = [0.91,0.41,0.17];
color(8,:) = [255,0,0]./255;
color(7,:) = [255,222,173]./255;
color(9,:) = [0,200,0]./255;

nVoxelsNet(1) = sum(mem_DAN);
nVoxelsNet(2) = sum(mem_VAN);
nVoxelsNet(3) = sum(mem_VIS);
nVoxelsNet(4) = sum(mem_LAN);
nVoxelsNet(5) = sum(mem_FPC);
nVoxelsNet(6) = sum(mem_SMN);
nVoxelsNet(7) = sum(mem_AUD);
nVoxelsNet(8) = sum(mem_DMN);
nVoxelsNet(9) = sum(mem_OUT);

for idx_clu=1:nTrig
    
    seed = trig(idx_clu).stat{1+1,1};
       
     for iCondition=1:nConditions
            
         for idx_net=1:nNets+1
             
           d = squeeze(mean(trig(idx_clu).m_e_sub(iCondition,:,:),2));  
           [s,i] = sort(d);
           max_Flag = i(end);
           
           vec = trig(idx_clu).stat{iCondition+1,idx_net+3};
           max_n_vec = vec(max_Flag,:);
         
           summary(idx_clu,iCondition,idx_net) = sum(max_n_vec);
            
         end
         
         
         
    end

end

nColors = 9;

for idx_clu=1:nTrig

    f = figure;
    iiPlot = 0;
    for iCondition=1:nConditions
 
        thisSeed = trig(idx_clu).stat{2,1};

        iiPlot = iiPlot + 1;

        title(int2str(thisSeed));
        
        iiBar = 0;
        subplot(1,3,iiPlot);
        for idx_net=1:nNets+1
            
            iiBar = iiBar + 1;

            b = bar(iiBar,summary(idx_clu,iCondition,idx_net));
            set(b,'FaceColor',color(idx_net,:));
            ylim([0 max(summary(:))]);
            hold on
 
        end

    end
    
    print(f,strcat('Ignition-v9-cluster_simple_model_oneNORM-',int2str(thisSeed),'.eps'),'-depsc');

end

func_vol = nifti(strcat('LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii'));
func_vol.dat.fname = 'Z:\_DATA\Parcellation\758-Cluster\LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii';
func_vol_img = func_vol.dat(:,:,:);

nROI = 90;
for idx_clu=1:nTrig

    for iCondition=1:nConditions
 
        thisSeed = trig(idx_clu).stat{2,1};
        
        d = squeeze(mean(trig(idx_clu).m_e_sub(iCondition,:,:),2));  
        [s,i] = sort(d);
        max_Flag = i(end);
        
        vec = trig(idx_clu).stat{iCondition+1,13};
        max_n_vec = vec(max_Flag,:);
        
        volume = zeros(MNI_size);
        
        iiCluster = 0;
        for iROI=1:nROI
            
            nClusters = ROI(iROI).nClusters;
            
            for iClu=1:nClusters
                
                iiCluster = iiCluster + 1;
                
                idx_voxels = ROI(iROI).clusters(iClu).idx_voxels;
                
                for iVoxel=1:length(idx_voxels)
                    
                    [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
                    
                    volume(idxx,idxy,idxz) = max_n_vec(iiCluster);
                    
                end
                
            end
            
        end
        
        nifti_file = func_vol;
        offset = func_vol.dat.offset;
        scl_slope = func_vol.dat.scl_slope;
        scl_inter = func_vol.dat.scl_inter;

        dtype = 'FLOAT32';
        offset = 0;

        dim = func_vol.dat.dim;

        descrip = 'Membership';

        fname = strcat('Ignition-v9-Membership-',int2str(iCondition),'-',int2str(thisSeed),'.nii');
        input_data = volume; 
        real_save_image;
   
    end
    
end

end

function checkSUBGRAPH_PEREVENT_CORBETTA_cluster(level)

disp('...loading biggest subgraph');
load(strcat('Ignition-v5-',level,'-All-0.04-0.07-biggest_subgraphs.mat'));

disp('...loading contrast');
load('Z:\_PAPERS\ignition\Events\analysis\Ignition-v1-cluster-Contrast.mat');

load('Z:\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\FC-Voxels-AAL-ROI-corr-KMeans\FC-Voxels-AAL-ROI-corr-KMeans-Info-Mean-TS-corrected.mat');

nROIs = 90;

disp('...functional networks');

networks = {'DAN' 'VAN' 'VIS' 'AUD' 'LAN' 'FPC' 'SMN' 'DMN'};

load('MD758-Corbetta-density-simple-model.mat');

mem_OUT = double( mem_DAN | mem_VAN | mem_VIS | mem_LAN | mem_FPC | mem_DMN | mem_SMN | mem_AUD );
mem_OUT = double(~mem_OUT);

nTotalRuns = 32;
nTotalVoxels = 160990;
nConditions = 3;
nTotalClusters = 758;
nTotalFlags = 4;

MNI_size = [91 109 91];

idx_resting = 1;
idx_passive = 2;
idx_attention = 3;

idx_att = find(contrast_attention_z);
idx_stim = find(contrast_stimulus_z);

all_idx_contrast = [idx_att,idx_stim];

for idx_clu=1:length(all_idx_contrast)
    
    thisSeed = all_idx_contrast(idx_clu);
    
    for iCondition=1:nConditions
        
        trig(idx_clu).stat{iCondition+1,1} = thisSeed;
        
        freq_membership = zeros(nTotalRuns,nTotalFlags,nTotalClusters);
    
        for iRun=1:nTotalRuns
            
            disp(strcat(int2str(idx_clu),':',int2str(iCondition),':',int2str(iRun)));

            nEvents(iRun) = size(biggest_subgraphs(iRun,iCondition,thisSeed).perEvent,1);
            
            for iFlag=1:nTotalFlags
                
                membership = zeros(1,nTotalClusters);
                real_membership = zeros(nEvents(iRun),nTotalClusters);
                
                for iEvent=1:nEvents(iRun)
                    
                    membership = membership + biggest_subgraphs(iRun,iCondition,thisSeed).perEvent{iEvent,iFlag}.membership;
                    
                    real_membership(iEvent,:) = biggest_subgraphs(iRun,iCondition,thisSeed).perEvent{iEvent,iFlag}.membership;
                    
                    subgraph_size(iCondition,iRun,iFlag,iEvent) = length(biggest_subgraphs(iRun,iCondition,thisSeed).perEvent{iEvent,iFlag}.subgraph);
                    
                end
                
                membership = membership ./ nEvents(iRun);

% % %                 freq_membership_DAN(iRun,iFlag,:) = ( ( mem_DAN ) ./ sum(mem_DAN) ) .* ( ( membership ) ./ sum(membership) );
% % %                 freq_membership_VAN(iRun,iFlag,:) = ( ( mem_VAN ) ./ sum(mem_VAN) ) .* ( ( membership ) ./ sum(membership) );
% % %                 freq_membership_VIS(iRun,iFlag,:) = ( ( mem_VIS ) ./ sum(mem_VIS) ) .* ( ( membership ) ./ sum(membership) );
% % %                 freq_membership_LAN(iRun,iFlag,:) = ( ( mem_LAN ) ./ sum(mem_LAN) ) .* ( ( membership ) ./ sum(membership) );
% % %                 freq_membership_DMN(iRun,iFlag,:) = ( ( mem_DMN ) ./ sum(mem_DMN) ) .* ( ( membership ) ./ sum(membership) );
% % %                 freq_membership_FPC(iRun,iFlag,:) = ( ( mem_FPC ) ./ sum(mem_FPC) ) .* ( ( membership ) ./ sum(membership) );
% % %                 freq_membership_SMN(iRun,iFlag,:) = ( ( mem_SMN ) ./ sum(mem_SMN) ) .* ( ( membership ) ./ sum(membership) );
% % %                 freq_membership_AUD(iRun,iFlag,:) = ( ( mem_AUD ) ./ sum(mem_AUD) ) .* ( ( membership ) ./ sum(membership) );
% % %                 freq_membership_OUT(iRun,iFlag,:) = ( ( mem_OUT ) ./ sum(mem_OUT) ) .* ( ( membership ) ./ sum(membership) );

                freq_membership_DAN(iRun,iFlag,:) = ( ( mem_DAN ) ./ sum(mem_DAN) ) .* ( ( membership ) );
                freq_membership_VAN(iRun,iFlag,:) = ( ( mem_VAN ) ./ sum(mem_VAN) ) .* ( ( membership )  );
                freq_membership_VIS(iRun,iFlag,:) = ( ( mem_VIS ) ./ sum(mem_VIS) ) .* ( ( membership ) );
                freq_membership_LAN(iRun,iFlag,:) = ( ( mem_LAN ) ./ sum(mem_LAN) ) .* ( ( membership )  );
                freq_membership_DMN(iRun,iFlag,:) = ( ( mem_DMN ) ./ sum(mem_DMN) ) .* ( ( membership )  );
                freq_membership_FPC(iRun,iFlag,:) = ( ( mem_FPC ) ./ sum(mem_FPC) ) .* ( ( membership )  );
                freq_membership_SMN(iRun,iFlag,:) = ( ( mem_SMN ) ./ sum(mem_SMN) ) .* ( ( membership ) );
                freq_membership_AUD(iRun,iFlag,:) = ( ( mem_AUD ) ./ sum(mem_AUD) ) .* ( ( membership )  );
                freq_membership_OUT(iRun,iFlag,:) = ( ( mem_OUT ) ./ sum(mem_OUT) ) .* ( ( membership )  );
                
                freq_membership_ALL(iRun,iFlag,:) = membership;
                
                for iEvent=1:nEvents(iRun)
                    
                    real_DAN(iEvent,:) = real_membership(iEvent,:) .* mem_DAN;
                    real_VAN(iEvent,:) = real_membership(iEvent,:) .* mem_VAN;
                    real_SMN(iEvent,:) = real_membership(iEvent,:) .* mem_SMN;
                    real_LAN(iEvent,:) = real_membership(iEvent,:) .* mem_LAN;
                    real_VIS(iEvent,:) = real_membership(iEvent,:) .* mem_VIS;
                    real_FPC(iEvent,:) = real_membership(iEvent,:) .* mem_FPC;
                    real_DMN(iEvent,:) = real_membership(iEvent,:) .* mem_DMN;
                    real_AUD(iEvent,:) = real_membership(iEvent,:) .* mem_AUD;
                    real_OUT(iEvent,:) = real_membership(iEvent,:) .* mem_OUT;
                    
                    real_ALL(iEvent,:) = real_membership(iEvent,:);
                    
                end
                
                reliability_DAN(iRun,iFlag) = getReliability(real_DAN > 0);
                reliability_VAN(iRun,iFlag) = getReliability(real_VAN > 0);
                reliability_VIS(iRun,iFlag) = getReliability(real_VIS > 0);
                reliability_LAN(iRun,iFlag) = getReliability(real_LAN > 0);
                reliability_DMN(iRun,iFlag) = getReliability(real_DMN > 0);
                reliability_FPC(iRun,iFlag) = getReliability(real_FPC > 0);
                reliability_SMN(iRun,iFlag) = getReliability(real_SMN > 0);
                reliability_AUD(iRun,iFlag) = getReliability(real_AUD > 0);
                reliability_OUT(iRun,iFlag) = getReliability(real_OUT > 0);
                
                reliability_ALL(iRun,iFlag) = getReliability(real_ALL > 0);
                
                clear real_DAN real_VAN real_SMN real_LAN real_VIS real_FPC real_DMN real_AUD real_OUT real_ALL
                
            end
            
        end
        
        trig(idx_clu).stat{iCondition+1,2} = mean(nEvents);
        trig(idx_clu).stat{iCondition+1,3} = std(nEvents);
        
        trig(idx_clu).stat{1,4} = 'mean-freq-DAN';
        trig(idx_clu).stat{1,5} = 'mean-freq-VAN';
        trig(idx_clu).stat{1,6} = 'mean-freq-VIS';
        trig(idx_clu).stat{1,7} = 'mean-freq-LAN';
        trig(idx_clu).stat{1,8} = 'mean-freq-FPC';
        trig(idx_clu).stat{1,9} = 'mean-freq-SMN';
        trig(idx_clu).stat{1,10} = 'mean-freq-AUD';
        trig(idx_clu).stat{1,11} = 'mean-freq-DMN';
        trig(idx_clu).stat{1,12} = 'mean-freq-OUT';
        trig(idx_clu).stat{1,13} = 'mean-freq-ALL';
        
        trig(idx_clu).stat{iCondition+1,4} = squeeze(mean(freq_membership_DAN,1));
        trig(idx_clu).stat{iCondition+1,5} = squeeze(mean(freq_membership_VAN,1));
        trig(idx_clu).stat{iCondition+1,6} = squeeze(mean(freq_membership_VIS,1));
        trig(idx_clu).stat{iCondition+1,7} = squeeze(mean(freq_membership_LAN,1));
        trig(idx_clu).stat{iCondition+1,8} = squeeze(mean(freq_membership_FPC,1));
        trig(idx_clu).stat{iCondition+1,9} = squeeze(mean(freq_membership_SMN,1));
        trig(idx_clu).stat{iCondition+1,10} = squeeze(mean(freq_membership_AUD,1));
        trig(idx_clu).stat{iCondition+1,11} = squeeze(mean(freq_membership_DMN,1));
        trig(idx_clu).stat{iCondition+1,12} = squeeze(mean(freq_membership_OUT,1));
        trig(idx_clu).stat{iCondition+1,13} = squeeze(mean(freq_membership_ALL,1));
        
        trig(idx_clu).stat{1,14} = 'mean-real-DAN';
        trig(idx_clu).stat{1,15} = 'mean-real-VAN';
        trig(idx_clu).stat{1,16} = 'mean-real-VIS';
        trig(idx_clu).stat{1,17} = 'mean-real-LAN';
        trig(idx_clu).stat{1,18} = 'mean-real-FPC';
        trig(idx_clu).stat{1,19} = 'mean-real-SMN';
        trig(idx_clu).stat{1,20} = 'mean-real-AUD';
        trig(idx_clu).stat{1,21} = 'mean-real-DMN';
        trig(idx_clu).stat{1,22} = 'mean-real-OUT';
        trig(idx_clu).stat{1,23} = 'mean-real-ALL';
        
        trig(idx_clu).stat{iCondition+1,14} = squeeze(mean(reliability_DAN,1));
        trig(idx_clu).stat{iCondition+1,15} = squeeze(mean(reliability_VAN,1));
        trig(idx_clu).stat{iCondition+1,16} = squeeze(mean(reliability_VIS,1));
        trig(idx_clu).stat{iCondition+1,17} = squeeze(mean(reliability_LAN,1));
        trig(idx_clu).stat{iCondition+1,18} = squeeze(mean(reliability_FPC,1));
        trig(idx_clu).stat{iCondition+1,19} = squeeze(mean(reliability_SMN,1));
        trig(idx_clu).stat{iCondition+1,20} = squeeze(mean(reliability_AUD,1));
        trig(idx_clu).stat{iCondition+1,21} = squeeze(mean(reliability_DMN,1));
        trig(idx_clu).stat{iCondition+1,22} = squeeze(mean(reliability_OUT,1));
        trig(idx_clu).stat{iCondition+1,23} = squeeze(mean(reliability_ALL,1));
        
        trig(idx_clu).stat{1,24} = 'std-freq-DAN';
        trig(idx_clu).stat{1,25} = 'std-freq-VAN';
        trig(idx_clu).stat{1,26} = 'std-freq-VIS';
        trig(idx_clu).stat{1,27} = 'std-freq-LAN';
        trig(idx_clu).stat{1,28} = 'std-freq-FPC';
        trig(idx_clu).stat{1,29} = 'std-freq-SMN';
        trig(idx_clu).stat{1,30} = 'std-freq-AUD';
        trig(idx_clu).stat{1,31} = 'std-freq-DMN';
        trig(idx_clu).stat{1,32} = 'std-freq-OUT';
        trig(idx_clu).stat{1,33} = 'std-freq-ALL';
        
        trig(idx_clu).stat{iCondition+1,24} = squeeze(std(freq_membership_DAN,0,1));
        trig(idx_clu).stat{iCondition+1,25} = squeeze(std(freq_membership_VAN,0,1));
        trig(idx_clu).stat{iCondition+1,26} = squeeze(std(freq_membership_VIS,0,1));
        trig(idx_clu).stat{iCondition+1,27} = squeeze(std(freq_membership_LAN,0,1));
        trig(idx_clu).stat{iCondition+1,28} = squeeze(std(freq_membership_FPC,0,1));
        trig(idx_clu).stat{iCondition+1,29} = squeeze(std(freq_membership_SMN,0,1));
        trig(idx_clu).stat{iCondition+1,30} = squeeze(std(freq_membership_AUD,0,1));
        trig(idx_clu).stat{iCondition+1,31} = squeeze(std(freq_membership_DMN,0,1));
        trig(idx_clu).stat{iCondition+1,32} = squeeze(std(freq_membership_OUT,0,1));
        trig(idx_clu).stat{iCondition+1,33} = squeeze(std(freq_membership_ALL,0,1));
        
        trig(idx_clu).stat{1,34} = 'std-real-DAN';
        trig(idx_clu).stat{1,35} = 'std-real-VAN';
        trig(idx_clu).stat{1,36} = 'std-real-VIS';
        trig(idx_clu).stat{1,37} = 'std-real-LAN';
        trig(idx_clu).stat{1,38} = 'std-real-FPC';
        trig(idx_clu).stat{1,39} = 'std-real-SMN';
        trig(idx_clu).stat{1,40} = 'std-real-AUD';
        trig(idx_clu).stat{1,41} = 'std-real-DMN';
        trig(idx_clu).stat{1,42} = 'std-real-OUT';
        trig(idx_clu).stat{1,43} = 'std-real-ALL';
        
        trig(idx_clu).stat{iCondition+1,34} = squeeze(std(reliability_DAN,0,1));
        trig(idx_clu).stat{iCondition+1,35} = squeeze(std(reliability_VAN,0,1));
        trig(idx_clu).stat{iCondition+1,36} = squeeze(std(reliability_VIS,0,1));
        trig(idx_clu).stat{iCondition+1,37} = squeeze(std(reliability_LAN,0,1));
        trig(idx_clu).stat{iCondition+1,38} = squeeze(std(reliability_FPC,0,1));
        trig(idx_clu).stat{iCondition+1,39} = squeeze(std(reliability_SMN,0,1));
        trig(idx_clu).stat{iCondition+1,40} = squeeze(std(reliability_AUD,0,1));
        trig(idx_clu).stat{iCondition+1,41} = squeeze(std(reliability_DMN,0,1));
        trig(idx_clu).stat{iCondition+1,42} = squeeze(std(reliability_OUT,0,1));
        trig(idx_clu).stat{iCondition+1,43} = squeeze(std(reliability_ALL,0,1));
        
        trig(idx_clu).stat{1,44} = 'freq-membership-ALL';
        trig(idx_clu).stat{iCondition+1,44} = freq_membership_ALL;
        
        m_e_sub = squeeze(mean(subgraph_size,4));
        
        trig(idx_clu).m_e_sub = m_e_sub;

        clear nEvents reliability_DAN reliability_VAN reliability_VIS reliability_LAN reliability_FPC reliability_SMN reliability_DMN reliability_AUD freq_membership_ALL
        
    end
    
    clear subgraph_size
    
end

save(strcat('Ignition-v10-',level,'-All-0.04-0.07-','checkSUBGRAPH_PEREVENT_CORBETTA_cluster_simpe_model_oneNORM','.mat'),'trig');

end

function getSUBGRAPH_PHASE

load('Ignition-v6-cluster-All-0.04-0.07-biggest_subgraphs.mat');
load('Ignition-v6-cluster-All-0.04-0.07-allOthers.mat');

nFlags = 4;
nRuns = 32;
idx_clusters = [90 108 112 242 339 549 604 605 606 607 608 623 627 699 101 229 317 527 666 727];
idx_rs = 1;
idx_pv = 2;
idx_at = 3;
max_flag = 1;
nTotalClusters = 758;
nConditions = 3;

iiFull = 0;
for iCluster=1:length(idx_clusters)
    
    this_cluster = idx_clusters(iCluster);
    
    SeedGraph(iCluster).this_cluster = this_cluster;
    
    for iCondition=1:nConditions
  
        clear graph 
        
        iiEvent = 0;
        for iRun=1:nRuns

            nEvents = size(biggest_subgraphs(iRun,iCondition,this_cluster).perEvent,1);
            
            [s,i] = sort(mean(squeeze(IntegStim2(iRun,iCondition,this_cluster,:,1:nEvents)),2),'Descend');
            
            max_flag = i(1);

            for iEvent=1:nEvents

                iiEvent = iiEvent + 1;

                graph(iiEvent,1:nTotalClusters) = 0;
                graph(iiEvent,biggest_subgraphs(iRun,iCondition,this_cluster).perEvent{iEvent,max_flag}.subgraph) = 1;
                
                tmp_graph_event(iEvent) = length(biggest_subgraphs(iRun,iCondition,this_cluster).perEvent{iEvent,max_flag}.subgraph);
                
            end
            
            [s,i] = sort(tmp_graph_event,'Descend');
            
            max_event = i(1);
            
            max_graph(iRun,1:nTotalClusters) = 0;
            max_graph(iRun,biggest_subgraphs(iRun,iCondition,this_cluster).perEvent{max_event,max_flag}.subgraph) = 1;
            
%             if sum(max_graph(iRun,:)) == nTotalClusters; 
%                 
%                 disp(strcat(int2str(this_cluster),':',int2str(iCondition),':',int2str(iRun),':',int2str(max_flag),':',int2str(max_event))); 
%                 
%                 iiFull = iiFull + 1;
%                 
%                 full_graph(iiFull,1) = this_cluster;
%                 full_graph(iiFull,2) = iCondition;
%                 full_graph(iiFull,3) = iRun;
%                 full_graph(iiFull,4) = max_flag;
%                 full_graph(iiFull,5) = max_event;
%             
%             end
            
            clear tmp_graph_event

        end
       
        SeedGraph(iCluster).Cond(iCondition).graph = graph;
        SeedGraph(iCluster).Cond(iCondition).max_graph = max_graph;
        
        clear graph max_graph
    
    end
    
end

for iCluster=1:length(idx_clusters)
    
    this_cluster = idx_clusters(iCluster);
    
    SeedGraph(iCluster).this_cluster = this_cluster;
    
    for iCondition=1:nConditions
  
        clear graph_per_flag 
        
        for iFlag=1:nFlags
            
            iiEvent = 0;
            for iRun=1:nRuns

                nEvents = size(biggest_subgraphs(iRun,iCondition,this_cluster).perEvent,1);
               
                for iEvent=1:nEvents
                    
                    iiEvent = iiEvent + 1;
                
                    graph_per_flag(iFlag).graph(iiEvent,1:nTotalClusters) = 0;
                    graph_per_flag(iFlag).graph(iiEvent,biggest_subgraphs(iRun,iCondition,this_cluster).perEvent{iEvent,iFlag}.subgraph) = 1;
                
                end
                
            end
                
        end
        
        SeedGraph(iCluster).Cond(iCondition).graph_per_flag = graph_per_flag;
        
    end
                
end

save('Ignition-v7-cluster-All-0.04-0.07-SeedGraph.mat','SeedGraph');

end

function plotGraphConcat

load('Ignition-v7-cluster-All-0.04-0.07-SeedGraph.mat');
load('Ignition-v6-cluster-All-0.04-0.07-allOthers.mat');

nConditions = 3;
nTotalClusters = 758;
MarkerSize = 12;

LineWidth = 4;

idx_rs = 1;
idx_pv = 2;
idx_at = 3;

cond_color{1} = 'b-';
cond_color{2} = 'r-';
cond_color{3} = 'k-';

for iSeed=1:length(SeedGraph)
    
    this_cluster = SeedGraph(iSeed).this_cluster;
    
    for iCondition=1:nConditions
    
        % graph = SeedGraph(iSeed).Cond(iCondition).max_graph;
        graph = SeedGraph(iSeed).Cond(iCondition).graph;
        
        iiEvent = size(graph,1);
        
        disp(strcat(int2str(this_cluster),':',int2str(iiEvent)));
        
        graph = double(~graph);

        f = figure;
        graph = graph'; imagesc(graph);
        hold on
% % %         for xG=1:iiEvent
% % %             for yG=1:nTotalClusters
% % %                 if graph(xG,yG); plot(xG+1,yG+1,'k*','MarkerSize',MarkerSize); hold on; end
% % %             end
% % %         end
        plot(10,this_cluster,'b*','MarkerSize',MarkerSize);
        colormap(gray);
        xlim([1 iiEvent]);
        ylim([1 nTotalClusters]);
        set(gca,'YTick',[],'XTick',[]);
        % box off
        % axis tight
        axis off
        print(f,strcat('Graph-',int2str(this_cluster),'-',int2str(iCondition),'.eps'),'-depsc');
      
        for iEvent=1:iiEvent
            this_seed_graph_size(iEvent) = sum(SeedGraph(iSeed).Cond(iCondition).graph(iEvent,:));
        end
        
        for iEvent=1:iiEvent
            this_seed_graph_size_per_run(iCondition,iEvent) = sum(SeedGraph(iSeed).Cond(iCondition).graph(iEvent,:));
        end
        
        f = figure;
        
        plot(1:iiEvent,this_seed_graph_size,cond_color{iCondition},'LineWidth',LineWidth);
        hold on
        plot(1:0.01:iiEvent,mean(this_seed_graph_size),cond_color{iCondition});
        xlim([1 iiEvent]);
        ylim([1 nTotalClusters]);
        set(gca,'YTick',[],'XTick',[]);
        % box off
        % axis tight
        axis off
        print(f,strcat('Graph-Size-',int2str(this_cluster),'-',int2str(iCondition),'.eps'),'-depsc');
        
        clear this_seed_graph_size
        
    end

   plotDistributionDensity(squeeze(this_seed_graph_size_per_run(idx_rs,:)),squeeze(this_seed_graph_size_per_run(idx_pv,:)),'RS','PV',strcat('Distribution of Graph Size','-',int2str(this_cluster),'-','RS'),'Graph Size','Frequency',cond_color{idx_rs},'k-',0,0);
  
   plotDistributionDensity(squeeze(this_seed_graph_size_per_run(idx_pv,:)),squeeze(this_seed_graph_size_per_run(idx_rs,:)),'PV','RS',strcat('Distribution of Graph Size','-',int2str(this_cluster),'-','PV'),'Graph Size','Frequency',cond_color{idx_pv},'k-',0,1);
    
   plotDistributionDensity(squeeze(this_seed_graph_size_per_run(idx_rs,:)),squeeze(this_seed_graph_size_per_run(idx_pv,:)),'RS','PV',strcat('Distribution of Graph Size','-',int2str(this_cluster),'-','RS-PV'),'Graph Size','Frequency',cond_color{idx_rs},cond_color{idx_pv},0,0);
    
   
   plotDistributionDensity(squeeze(this_seed_graph_size_per_run(idx_at,:)),squeeze(this_seed_graph_size_per_run(idx_pv,:)),'AT','PV',strcat('Distribution of Graph Size','-',int2str(this_cluster),'-','AT-PV'),'Graph Size','Frequency',cond_color{idx_at},cond_color{idx_pv},0,0);

    stim_label = strcat('Distribution of Graph Size','-',int2str(this_cluster),'-','RS-PV');
    att_label = strcat('Distribution of Graph Size','-',int2str(this_cluster),'-','AT-PV');
    
    stimulus = [squeeze(mevokedinteg2(:,idx_rs,this_cluster)),squeeze(mevokedinteg2(:,idx_pv,this_cluster))];
    attention = [squeeze(mevokedinteg2(:,idx_at,this_cluster)),squeeze(mevokedinteg2(:,idx_pv,this_cluster))];
    
    plotDistributionBoxPlot(stimulus,{'RS' 'PV'},stim_label);
    plotDistributionBoxPlot(attention,{'AT' 'PV'},att_label);

end

end

function plotFullGraphTS

load('Ignition-v5-cluster-All-0.04-0.07-SeedGraph.mat');
load('Z:\_PAPERS\ignition\Events\analysis\Ignition-v5-cluster-All-0.04-0.07-allOthers.mat');

nTR = 150;
nTotalClusters = 758;
nSubjects = 8;
nROI = 90;
nRuns = 4;
nTotalRuns = 32;
nConditions = 3; %%% 1 = RestingState, 2 = PassiveViewing, 3 = Attentive Tracking
TR = 2; %sampling interval(TR)  
nTotalVoxels = 160990;
MNI_size = [91 109 91];

integration_cutoff = 10;
integration_cutoff2 = 9;

disp('getting Data...');

% all_data = getFormattedData('cluster','All');

nFullGraphs = size(full_graph,1);

for iFullGraph=1:nFullGraphs
    
    clear this_data
   
    iCond = full_graph(iFullGraph,2);
    iRun = full_graph(iFullGraph,3);
    
    events_ts = squeeze(Events(iRun,iCond,:,:));
    
    s_events_ts = sum(events_ts,1);
    
    max(s_events_ts)
    
end

s_tise = std(tise,0,1);
m_tise = mean(tise,1);
ms_tise = s_tise + m_tise;


end

function getGraphReliability

load('Ignition-v7-cluster-All-0.04-0.07-SeedGraph.mat');

nConditions = 3;
nTotalClusters = 758;
MarkerSize = 12;

nFlags = 4;

idx_rs = 1;
idx_pv = 2;
idx_at = 3;

cond_color{1} = 'b-';
cond_color{2} = 'r-';
cond_color{3} = 'k-';

cond_color_star{1} = 'b*';
cond_color_star{2} = 'r*';
cond_color_star{3} = 'k*';

for iSeed=1:length(SeedGraph)
    
    this_cluster = SeedGraph(iSeed).this_cluster;
    
    SeedGraphReal(iSeed).this_cluster = SeedGraph(iSeed).this_cluster;
    
    for iCondition=1:nConditions
   
        for iFlag=1:nFlags
        
            graph = SeedGraph(iSeed).Cond(iCondition).graph_per_flag(iFlag).graph;
           
            SeedGraphReal(iSeed).Cond(iCondition).real(iFlag) = getReliability(graph);
            
        end
      
    end
    
end

save('Ignition-v7-cluster-All-0.04-0.07-SeedGraphReliability.mat','SeedGraphReal');

end

function plotGraphReliability

load('Ignition-v7-cluster-All-0.04-0.07-SeedGraphReliability.mat');

nConditions = 3;
nFlags = 4;

cond_color{1} = 'b-';
cond_color{2} = 'r-';
cond_color{3} = 'k-';

cond_color_star{1} = 'b*';
cond_color_star{2} = 'r*';
cond_color_star{3} = 'k*';

MarkerSize = 12;
LineWidth = 4;

idx_rs = 1;
idx_pv = 2;
idx_at = 3;

for iSeed=1:length(SeedGraphReal)

    this_cluster = SeedGraphReal(iSeed).this_cluster;

    %for iCondition=1:nConditions
        
        f = figure;

        plot(1:nFlags,SeedGraphReal(iSeed).Cond(idx_rs).real,cond_color{idx_rs},'MarkerSize',MarkerSize,'LineWidth',LineWidth);
        hold on
        plot(1:nFlags,SeedGraphReal(iSeed).Cond(idx_rs).real,cond_color_star{idx_rs},'MarkerSize',MarkerSize,'LineWidth',LineWidth);
        
        plot(1:nFlags,SeedGraphReal(iSeed).Cond(idx_pv).real,cond_color{idx_pv},'MarkerSize',MarkerSize,'LineWidth',LineWidth);
        hold on
        plot(1:nFlags,SeedGraphReal(iSeed).Cond(idx_pv).real,cond_color_star{idx_pv},'MarkerSize',MarkerSize,'LineWidth',LineWidth);

        set(gca,'XTick',1:nFlags);
        xlim([0 nFlags+1]);
        ylim([min( [min(SeedGraphReal(iSeed).Cond(idx_rs).real(:)),min(SeedGraphReal(iSeed).Cond(idx_pv).real(:))] ) max( [max(SeedGraphReal(iSeed).Cond(idx_rs).real(:)),max(SeedGraphReal(iSeed).Cond(idx_pv).real(:))] ) ]);
        set(gca,'YTick',[min( [min(SeedGraphReal(iSeed).Cond(idx_rs).real(:)),min(SeedGraphReal(iSeed).Cond(idx_pv).real(:))] ) max( [max(SeedGraphReal(iSeed).Cond(idx_rs).real(:)),max(SeedGraphReal(iSeed).Cond(idx_pv).real(:))] ) ]);
        % axis off
        axis tight

        print(f,strcat('Reliability-',int2str(this_cluster),'-','RS-PV','-perFlag.eps'),'-depsc');
    
    %end
    
end

end

function u_labels = getROILabels(idx_clusters)

ROI_info = getROIInfoOnClusters;

labels = cell.empty;
nROIs = 90;

iLabel = 0;
for iClu=1:length(idx_clusters)
    
    idx = idx_clusters(iClu);
  
    for iROI=1:nROIs

        if idx >= ROI_info{iROI,4} && idx <= ROI_info{iROI,5}
            
            iLabel = iLabel + 1;
            
            labels{iLabel} = ROI_info{iROI,2};
            
        end

    end
    
end

if ~isempty(labels)
    u_labels = unique(labels);
else
    u_labels = cell.empty;
end

end

function plot3DclustersFromContrast

level = 'cluster';
load(strcat('Ignition-v6-',level,'-All-0.04-0.07-biggest_subgraphs-Contrast.mat'));

func_vol = nifti('LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
func_vol.dat.fname = 'Z:\_DATA\Parcellation\758-Cluster\LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii';
func_img = func_vol.dat(:,:,:);

zthreshold = 2.3;
MNI_size = [91 109 91];
nTotalClusters = 758;

nTrig = length(Contrast);

for iTrig=1:nTrig
    
    thisSeed = Contrast(iTrig).thisSeed;
    
    nifti_file = func_vol;
    offset = func_vol.dat.offset;
    scl_slope = func_vol.dat.scl_slope;
    scl_inter = func_vol.dat.scl_inter;
    dtype = 'FLOAT32';
    offset = 0;
    dim = func_vol.dat.dim;
    
    volume = zeros(MNI_size);
    
     for iCluster=1:nTotalClusters
        
        idx_voxels = find(func_img==iCluster);
        
        nVoxels = length(idx_voxels);
        
        for iVoxel=1:nVoxels
            
            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
            
            if iCluster==thisSeed
                
                volume(idxx,idxy,idxz) = iCluster;
                
            end
            
        end
        
     end
     
    descrip = 'Cluster';
    fname = strcat(int2str(thisSeed),'.nii');
    input_data = volume; 
    real_save_image;
     
end
    
end
