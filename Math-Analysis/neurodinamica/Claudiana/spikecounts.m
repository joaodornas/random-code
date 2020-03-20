function spikeCountMatrix=spikecounts(spike_times,stimIds,parameters)

%function: spikeCountMatrix=spikecounts(spike_times,stimIds,parameters)
%creates a matrix for ROC analysis containing spike counts within the 
%stimulation period, where rows are conditions and columns are trials 
%(which will be assigned NaN if non existent). 
%Inputs are default Dynamo formats for spike time matrix (spike_times),
%stimIds vector and parameters structure, containing the following fields: 
%nconditions (default 11), samplingFreq (default 32000 Hz), baseline_duration
%(default 1000 ms), stim_duration (default 4000 ms)
%
%by Lucas Pinto in June 2009
%modified by Pedro Vieira in October 2009

%get function parameters or assign default values
if nargin<2
    disp('error: i need more input arguments!');
elseif nargin == 2
    nconditions=48;
    samplingFreq=32000;
    baseline_duration=500;
    stim_duration=2000;
elseif nargin == 3
    nconditions=parameters.nconditions;
    samplingFreq=parameters.samplingFreq;
    baseline_duration=parameters.baseline_duration;
    stim_duration=parameters.stim_duration;
end

%initialize spike count Matrix, where rows are conditions and columns are
%trials, which will be assigned NaN if non existent
nTrialsPerCond=numel(stimIds)/nconditions;
spikeCountMatrix=nan(nconditions,nTrialsPerCond);

%count spikes
for i=1:nconditions
    
    %create matrix for condition i
    respMatrix=spike_times(find(stimIds==i),:)./samplingFreq;

    ntrials_iCond=size(respMatrix,1);
    
    %for each trial count # spikes within situmulation period
    for j=1:ntrials_iCond
        
        spikes_ij=find(respMatrix(j,:)>baseline_duration/1000 & respMatrix(j,:)<=(baseline_duration+stim_duration)/1000);
        
        if ~isempty(spikes_ij)
            spikeCountMatrix(i,j)=numel(spikes_ij);
        else spikeCountMatrix(i,j)=0;
        end
    end
        %check for conditions with ntrials_iCond > nTrialsPerCond and 
        %replace with NaN
        %maxNtrials=size(spikeCountMatrix,2);
        %if ntrials_iCond < maxNtrials
        %    spikeCountMatrix(t,ntrials_iCond+1:maxNtrials)=nan;
        %end
        
   
    
end
%check for conditions with ntrials_iCond > nTrialsPerCond and 
%replace with NaN
for t=1:nconditions
    respMatrix=spike_times(find(stimIds==t),:)./samplingFreq;
    ntrials_iCond=size(respMatrix,1);
    maxNtrials=size(spikeCountMatrix,2);
    if ntrials_iCond < maxNtrials
            spikeCountMatrix(t,ntrials_iCond+1:maxNtrials)=nan;
    end
end