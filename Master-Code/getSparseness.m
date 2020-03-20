function getSparseness(spassDataSetName,samplingFrequency,startTime,endTime,binSize)

%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function loads a data set of electrophysiology coming from Spass
% software format and computes the Sparseness parameter along trials
% following the definition as described on the following reference paper [1,2].
%
% INPUT:
%      1. spassDataSetName: the name of Spass data set file.
%      2. samplingFrenquency: the sampling frenquency used during
%      electrophysiological recordings.
%      3. startTime: the start time of the trial in milliseconds.
%      4. endTime: the end time of the trial in milliseconds.
%      5. binSize: the size of the bin to compute spike rate, in
%      milliseconds.
%
% OUTPUT:
%      1.sparseness: structure with three variables:
%               a. gallant2002Sparseness: Sparseness paremeter as described
%               on references papers, for each condition.
%               b. trialsSpikes: the spikes for all trials in each
%               condition.
%               c. nTrials: the amount of trials, per condition.
%
% [1] Vinje, W., & Gallant, J. (2000). Sparse coding and decorrelation in 
% primary visual cortex during natural vision. Science (New York, NY), 287(5456), 1273.
%
% [2] Vinje, W., & Gallant, J. (2002). Natural stimulation of the nonclassical receptive 
% field increases information transmission efficiency in V1. Journal of Neuroscience, 22(7), 2904.
%
%
% Author: João V. Dornas, joaodornas@gmail.com, 03/2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

disp('BEGIN');

latency_startTime = 500;

latency_endTime = latency_startTime + 300;

p = 10000;

idx = strfind(spassDataSetName,'\');

registro = spassDataSetName(idx(length(idx))+1:length(spassDataSetName));

%%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Spass...');

    Spass = load(strcat(spassDataSetName));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% LOAD CONDITIONS TRIALS LABELS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Conditions Trials Labels...');

    nConditions = max(Spass.stimIds);
    
    for iCond=1:nConditions
   
        trialsLabels(iCond).label = find(Spass.stimIds == iCond); 
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%% SET RESOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Set Resolution...');

    spikeTimes = Spass.spike_times;
    
    spikeTimes = spikeTimes./ samplingFrequency;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% BEGIN CONDITIONS LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Begin Conditions Loop...');

    for iCond=1:nConditions
       
       disp(strcat('Begin Condition:',int2str(iCond))); 
        
       trialsSpikes = nan(size(spikeTimes(trialsLabels(iCond).label,:)));
       
       trialsSpikes = spikeTimes(trialsLabels(iCond).label,:);
        
       nTrials = size(trialsSpikes,1);
       
       alltrials = trialsSpikes;

       alltrials = reshape(alltrials.',[],1);

       alltrials = alltrials.'; 
       
       latency_time_data = load(strcat(registro,'-latency.mat'));
       latency_time = latency_time_data.latency_time;
       
       spikesVector = reshape(trialsSpikes.',[],1);

       spikesVector = spikesVector.';
       
       spikesVector = sort(spikesVector);
       
       spikesVector = spikesVector(spikesVector>0);

       spikesVector = spikesVector(spikesVector>(startTime/1000) & spikesVector<(endTime/1000)) - latency_time;
       
       nBins = (endTime - startTime)/binSize;
       
       for iBin=1:nBins
          
           nSpikes = length(spikesVector(spikesVector>=((iBin-1)*binSize/1000 + startTime/1000) & spikesVector<(iBin*binSize/1000 + startTime/1000)));
           
           sparseData(iCond).binRate(iBin) = nSpikes/(nTrials*binSize/1000);           
           sparseData(iCond).binRateSq(iBin) = ( sparseData(iCond).binRate(iBin) )^2;  
                     
       end
  
       disp('...Calculate Sparseness');
       
       sparseData(iCond).gallant2002Sparseness = ( 1 - ( mean(sparseData(iCond).binRate)^2 / ( mean(sparseData(iCond).binRate)^2 + std(sparseData(iCond).binRate)^2 ) ) / ( 1 - ( 1 / nBins ) ) ) ;
       
       sparseData(iCond).trialsSpikes = trialsSpikes;
       sparseData(iCond).nTrials = nTrials;
           
       disp(strcat('End Condition . ',int2str(iCond)));

    end  
    
disp('End Conditions Loop');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Save data...');
       
    save(strcat(registro,'-sparseness.mat'),'sparseData','samplingFrequency','startTime','endTime','binSize');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('END');
    
toc
    
end