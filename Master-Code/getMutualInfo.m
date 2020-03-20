
function getMutualInfo(spassDataSetName,samplingFrequency,startTime,endTime,binSize)

%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function loads a data set of electrophysiology coming from Spass
% software format and computes the Mutual Information for each condition 
% folloewing the definition as described on the reference paper [1]. For
% statistical inference about the Entropies calculation, we used the MA 
% bound [2].
%
% INPUT:
%      1. spassDataSetName: the name of Spass data set file.
%      2. samplingFrenquency: the sampling frenquency used during
%      electrophysiological recordings.
%      3. startTime: the start time of the trial in milliseconds.
%      4. endTime: the end time of the trial in milliseconds.
%      5. latencyStartTime: the start point of latency.
%      6. latencyTime: the latency duration.
% 
% OUTPUT:
%      1.info: structure with one substructure per condition and one
%      variable:
%               1.a condition: substructure for each condition with the
%               following variables:
%               
%               In a spike train a "Word" is defined as the Spike Rate
%               inside a bin of a specific size. So, for instance, for a
%               bin size of 1 ms we can have two words: 0 and 1.
%
%               For each size of bin, we have the following variables:
%
%                	WordMatrix: a matrix with all words in a condition
%                               across time (columns) and trials (rows).
%                	WordsUnique: unique words over all trials.
%                	nWordsUnique: amount of unique words.
%                   Frequencies: the frequency of each word over all trials.
%                   CondWordsUnique: unique words conditioned to a point in
%                                   time.
%                   nCondWordsUnique: amount to unique words conditioned to
%                                   a point in time.
%                   CondFrequencies: the frequencies of each conditioned
%                   word for each point in time.
%
%                	HResponse: the Entropy over all trials
%                	HResponseMAbound: the MA Entropy
%                	HCondResponse: the Entropy conditioned to a point in
%                                   time.
%                	HCondResponseMAbound: the MA conditoned Entropy
%
%                	meanHCondResponse: mean conditional Entropy
%                	stdHCondResponse: standard deviation of conditional
%                                     Entropy
%                	meanHCondResponseMAbound: the mean of MA conditional 
%                                              entropies
%                                              
%                	stdHCondResponseMAbound: the standard deviation of 
%                                            MA conditional Entropies
%
%                	MutualInfo: the Mutual Information 
%                	CodingEfficiency: the Coding Efficiency
%
%               Considering all values across bins sizes, we have:
%
%                   MaxMutualInfoBits: the maximum value of Mutual
%                                      Information
%                   MaxMutualInfoResolution: the bin size of the maximum
%                                      value of Mutual Information
%                   MaxMutualInfoCodingEfficiency: the Coding Efficiency
%                                      for the bin size with the maximum 
%                                      value of Mutual Information
%                	MaxCodingEfficiency: the maximum value for Coding
%                	                     Efficiency
%                   MaxCodingEfficiencyResolution: the bin size of the
%                                        maximum value of Coding Efficiency
%            
%               1.b sizes_of_a_bin: the sizes of each bin in miliseconds.
%
%
%
% [1] Rob R. de Ruyter van Steveninck, Georey D. Lewen, S. P. Strong,
%     R. Koberle and William Bialek, Science 275, 1805 (1997).
% [2] S. K. Ma, J. Stat. Phys. 26, 221 (1981).
%
% Author: Joao V. Dornas, joaodornas@gmail.com, 04/2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = strfind(spassDataSetName,'\');

registro = spassDataSetName(idx(length(idx))+1:length(spassDataSetName));

latency_start_time = 500;

latency_end_time = latency_start_time + 300;

%% BIN SIZE

sizesOfABin = 1:(endTime-startTime);

nBinsSizes = length(sizesOfABin);

%% CONDITIONS
    
Spass = load(strcat(spassDataSetName));
    
nConditions = Spass.parameters.nconditions;
    
for iCondition=1:nConditions
   
    trialsLabel(iCondition).labels = find(Spass.stimIds == iCondition); 
    
end
    
spikeTimes = Spass.spike_times ./ samplingFrequency;
    
%% EACH LOOP FOR A CONDITION

for iCondition=1:nConditions

    disp(strcat('Condition:',int2str(iCondition)));

    allTrials = spikeTimes(trialsLabel(iCondition).labels(:),:);

    nTrials = min(size(allTrials,1),size(allTrials,2));
    
    latency_time_data = load(strcat(registro,'-latency.mat'));
    latency_time = latency_time_data.latency_time;
    
    for t=1:nTrials

        Trial(t).spike_train = allTrials(t,:);

    	a_e = Trial(t).spike_train(Trial(t).spike_train>=(0/1000) & Trial(t).spike_train<(latency_start_time/1000));

        after_ae = Trial(t).spike_train(Trial(t).spike_train>=(latency_start_time/1000) & Trial(t).spike_train<(endTime/1000)) - latency_time;

        Trial(t).spike_train = [];

        Trial(t).spike_train = [a_e, after_ae];

        Trial(t).spike_train = Trial(t).spike_train(Trial(t).spike_train>0);

        Trial(t).spike_train = Trial(t).spike_train(Trial(t).spike_train>=(startTime/1000) & Trial(t).spike_train<=(endTime/1000));
    
        nSpikes(t) = numel(Trial(t).spike_train);

    end

    maxNSpikes = max(nSpikes);

    spike_trains = zeros(t,maxNSpikes);

    for t=1:nTrials
    
        spike_trains(t,:) = [Trial(t).spike_train, zeros(1,maxNSpikes - numel(Trial(t).spike_train))];
    
    end

    spikeTrains = spike_trains;
    
    WordMatrix(1:nBinsSizes) = struct('wm',[]);
    WordsUnique(1:nBinsSizes) = struct('w',[]);
    nWordsUnique(1:nBinsSizes) = struct('nw',[]);
    Frequencies(1:nBinsSizes) = struct('f',[]);
    CondWordsUnique(1:nBinsSizes) = struct('cw',[]);
    nCondWordsUnique(1:nBinsSizes) = struct('cnw',[]);
    CondFrequencies(1:nBinsSizes) = struct('cf',[]);

    HResponse = zeros(1,nBinsSizes);
    HResponseMAbound = zeros(1,nBinsSizes);
    HCondResponse(1:nBinsSizes) = struct('hcr',[]);
    HCondResponseMAbound(1:nBinsSizes) = struct('hcr',[]);

    meanHCondResponse = zeros(1,nBinsSizes);
    stdHCondResponse = zeros(1,nBinsSizes);
    meanHCondResponseMAbound = zeros(1,nBinsSizes);
    stdHCondResponseMAbound = zeros(1,nBinsSizes);

    MutualInfo = zeros(1,nBinsSizes);
    CodingEfficiency = zeros(1,nBinsSizes);

    for bin=1:nBinsSizes

        disp(strcat('binSize:',int2str(sizesOfABin(bin))));

        nBins = ceil( ( endTime - startTime ) / sizesOfABin(bin) );

        disp(strcat('getWordMatrix'));

        WordMatrix(bin).wm = getWordMatrix(spikeTrains,nTrials,nBins,sizesOfABin(bin),startTime);

        disp(strcat('getWordFrequencies'));

        [WordsUnique(bin).w, nWordsUnique(bin).nw, Frequencies(bin).f, CondWordsUnique(bin).cw, nCondWordsUnique(bin).cnw, CondFrequencies(bin).cf] = getWordFrequencies(WordMatrix(bin).wm);

        disp('getEntropy');

        [HResponse(bin), HResponseMAbound(bin)] = getEntropy(Frequencies(bin).f);

        disp('getCondEntropy');

        [HCondResponse(bin).hcr, HCondResponseMAbound(bin).hcr] = getCondEntropy(CondFrequencies(bin).cf); 

        meanHCondResponse(bin) = mean(HCondResponse(bin).hcr);
        stdHCondResponse(bin) = std(HCondResponse(bin).hcr);

        meanHCondResponseMAbound(bin) = mean(HCondResponseMAbound(bin).hcr);
        stdHCondResponseMAbound(bin) = std(HCondResponseMAbound(bin).hcr);

        MutualInfo(bin) = HResponse(bin) - meanHCondResponse(bin);
        CodingEfficiency(bin) = ( ( HResponse(bin) - meanHCondResponse(bin) ) / HResponse(bin) ) * 100;

    end

    info.Condition(iCondition).WordMatrix = WordMatrix;

    info.Condition(iCondition).WordsUnique = WordsUnique;
    info.Condition(iCondition).nWordsUnique = nWordsUnique;
    info.Condition(iCondition).Frequencies = Frequencies;
    info.Condition(iCondition).CondWordsUnique = CondWordsUnique;
    info.Condition(iCondition).nCondWordsUnique = nCondWordsUnique;
    info.Condition(iCondition).CondFrequencies = CondFrequencies;

    info.Condition(iCondition).HResponse = HResponse;
    info.Condition(iCondition).HResponseMAbound = HResponseMAbound;

    info.Condition(iCondition).HCondResponse = HCondResponse;
    info.Condition(iCondition).HCondResponseMAbound = HCondResponseMAbound;

    info.Condition(iCondition).meanHCondResponse = meanHCondResponse;
    info.Condition(iCondition).stdHCondResponse = stdHCondResponse;

    info.Condition(iCondition).meanHCondResponseMAbound = meanHCondResponseMAbound;
    info.Condition(iCondition).stdHCondResponseMAbound = stdHCondResponseMAbound;

    info.Condition(iCondition).MutualInfo = MutualInfo;
    info.Condition(iCondition).CodingEfficiency = CodingEfficiency;

    for b=1:length(sizesOfABin)

        MutualInfo(b) = info.Condition(iCondition).MutualInfo(b);
        CodingEfficiency(b) = info.Condition(iCondition).CodingEfficiency(b);

    end

    info.Condition(iCondition).MaxMutualInfoBits = max(MutualInfo);

    info.Condition(iCondition).MaxMutualInfoResolution = find(MutualInfo==max(MutualInfo));

    info.Condition(iCondition).MaxMutualInfoCodingEfficiency = CodingEfficiency(find(MutualInfo==max(MutualInfo)));

    info.Condition(iCondition).MaxCodingEfficiency = max(CodingEfficiency);

    info.Condition(iCondition).MaxCodingEfficiencyResolution = find(CodingEfficiency==max(CodingEfficiency));


    clear WordMatrix;
    clear WordsUnique;
    clear nWordsUnique;
    clear Frequencies;
    clear CondWordsUnique;
    clear nCondWordsUnique;
    clear CondFrequencies;

    clear HResponse;
    clear HResponseMAbound;
    clear HCondResponse;
    clear HCondResponseMAbound;

    clear meanHCondResponse;
    clear stdHCondResponse;
    clear meanHCondResponseMAbound;
    clear stdHCondResponseMAbound;

    clear MutualInfo;
    clear CodingEfficiency;

end

info.sizes_of_a_bin = sizesOfABin;

save(strcat(registro,'-mutual-information.mat'),'info','samplingFrequency','startTime','endTime','binSize');
    
end

%%% Gets Conditional Entropy 
function [ HCond, HCondMAbound ] = getCondEntropy(CondFrequencies)

for i=1:length(CondFrequencies)

    frequencies = CondFrequencies(i).freq;

    [ H, HMAbound ] = getEntropy(frequencies);
    
    HCond(i) = H;
    
    HCondMAbound(i) = HMAbound;
    
    clear H;
    clear HMAbound;
    
end


end

%%% Gets Entropy
function [ H, HMAbound ] = getEntropy(frequencies)

H = 0;

for f=1:length(frequencies)

   H = H - frequencies(f)*log2(frequencies(f));

end

sumOfFrequencies = 0;
for f=1:length(frequencies)

   sumOfFrequencies = sumOfFrequencies + frequencies(f)^2;

end

HMAbound = - log2(sumOfFrequencies);

end

%%% Gets Words and Their Frequencies
function [WordsUnique, nWordsUnique, frequencies, CondWordsUnique, nCondWordsUnique, CondFrequencies] = getWordFrequencies(WordMatrix)
   
nTotalWords = size(WordMatrix,1) * size(WordMatrix,2);

ConditionalWords = size(WordMatrix,2);

allWords = reshape(WordMatrix,1,[]);

[WordsUnique, nWordsUnique] = count_unique(allWords);

frequencies = nWordsUnique ./ nTotalWords;

for i=1:size(WordMatrix,2)
    
    binGroup = WordMatrix(:,i);
    
    [CondWordsUnique(i).words, nCondWordsUnique(i).nWords] = count_unique(binGroup);
    
    CondFrequencies(i).freq = nCondWordsUnique(i).nWords./ size(WordMatrix,1) ;
    
end

end

%%% Builds the Word Matrix
function WordMatrix = getWordMatrix(spikeTrains,nTrials,nBins,binSize,startTime)

scale = 1000;

nWords = nBins;

WordMatrix = [];

for t=1:nTrials
    
   train = spikeTrains(t,:);
   
   sentence = [];
   
   for w=1:nWords
       
       word = numel(train(train>=startTime/scale + (w-1)*binSize/scale & train<=startTime/scale + (w)*binSize/scale));
       
       sentence = [sentence, word];
       
       word = 0;
       
   end
   
   WordMatrix = [WordMatrix; sentence];
    
end

end

%%% Remove Latency from All Spike Trains
function spikeTrains = takeLatencyOff(allTrials,nTrials,startTime,endTime,latencyStartTime,latencyTime)

scale = 1000;

for t=1:nTrials

    Trial(t).spikeTrain = allTrials(t,:);

    a_e = Trial(t).spikeTrain(Trial(t).spikeTrain>=(0/scale) & Trial(t).spikeTrain<(latencyStartTime/scale));

    after_ae = Trial(t).spikeTrain(Trial(t).spikeTrain>=(latencyStartTime/scale) & Trial(t).spikeTrain<(endTime/scale)) - latencyTime/scale;

    Trial(t).spikeTrain = [];

    Trial(t).spikeTrain = [a_e, after_ae];

    Trial(t).spikeTrain = Trial(t).spikeTrain(Trial(t).spikeTrain>0);

    Trial(t).spikeTrain = Trial(t).spikeTrain(Trial(t).spikeTrain>=(startTime/scale) & Trial(t).spikeTrain<=(endTime/scale));
    
    nSpikes(t) = numel(Trial(t).spikeTrain);

end

maxNSpikes = max(nSpikes);

spikeTrains = zeros(t,maxNSpikes);

for t=1:nTrials
    
    spikeTrains(t,:) = [Trial(t).spikeTrain, zeros(1,maxNSpikes - numel(Trial(t).spikeTrain))];
    
end

end


