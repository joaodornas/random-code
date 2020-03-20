% function sparsenessAlongAllTrials(date,site_index,channel,registro,video_index,start_time,end_time,bin_size,nConditions)

function getSparseness(spassDataSetName,samplingFrequency,startTime,endTime,totalTime,binSize)

%%% DESCRIPTION
%
% This function loads a data set of electrophysiology coming from Spass
% software format and computes the Sparseness parameter along trials
% following the definition as described on the following reference paper [1].
%
% [1] Vinje, W., & Gallant, J. (2000). Sparse coding and decorrelation in 
% primary visual cortex during natural vision. Science (New York, NY), 287(5456), 1273.
%
%
% Author: João V. Dornas, joaodornas@gmail.com, 03/2017.

tic;

disp('BEGIN');

%%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Spass...');

    Spass = load(strcat(spassDataSetName,'.mat'));

    allSpikeTrainsMatrix = Spass.spike_times ./ samplingFrequency;
    
    nTrials = size(allSpikeTrainsMatrix,1);
    
    for iTrial=1:nTrials
        
        trains(iTrial).vectors = allSpikeTrainsMatrix(iTrial,:);
        
        trains(iTrial).vectors = trains(iTrial).vectors(trains(iTrial).vectors>0);
        
        trains(iTrial).vectors = trains(iTrial).vectors(trains(iTrial).vectors>=(startTime/1000) & trains(iTrial).vectors<=(endTime/1000));
        
    end
    
        
    allSpikeTrainsVector = 0;
    
    for iTrial=1:nTrials
        
        allSpikeTrainsVector = [allSpikeTrainsVector (trains(iTrial).vectors + (totalTime*(iTrial-1)))];
        
    end
    
    allSpikeTrainsVector(1) = [];
    
    nBins = nTrials * (endTime - startTime) / (binSize/1000);
       
    for iBin=1:nBins
          
        spikes = length(allSpikeTrainsVector(allSpikeTrainsVector>=((iBin-1)*binSize/1000 + startTime/1000) & allSpikeTrainsVector<(iBin*binSize/1000 + startTime/1000)));

        sparseData.binRate(iBin) = spikes/(binSize/1000);

    end
    
    sparseData.SparsenessInTime = 0;
       
    for iBin=1:length(sparseData.binRate)
       
        sparseData.SparsenessInTime = [sparseData.SparsenessInTime ( 1 - ( mean(sparseData.binRate(1:iBin))^2 / ( mean(sparseData.binRate(1:iBin))^2 + std(sparseData.binRate(1:iBin))^2 ) ) / 1 - ( 1 / nBins ) )];
       
    end
    
    
    plot(1:length(sparseData.SparsenessInTime),sparseData.SparsenessInTime)  
    
end