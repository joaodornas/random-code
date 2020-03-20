% function shapeForRasterandPSTH(spike_times,trialsize,trialsPerCondition,samplingFreq)
% 
% % shapes matrices for each condition for raster plots. Inputs are spike_time, 
% % which is a matrix with raw spike times (from labview vi's, trial
% % size in ms, trials per condition and sampling frequency. Defaults are
% % 7000, 10 and 32000, respectively.
% % shapeForRasterandPSTH(spike_times,trialsize,trialsPerCondition,samplingFreq)
% 
% if nargin==1
%     trialsize=7000;
%     trialsPerCondition=10;
%     samplingFreq=32000;
% elseif nargin==2
%     trialsPerCondition=10;
%     samplingFreq=32000;
% elseif nargin==3
%     samplingFreq=32000;
% end

ntrials=size(spike_times,1);
nconditions=ntrials/trialsPerCondition;
maxnspikes=size(spike_times,2);


%%%%%%%%%%%%% organize in condition matrices %%%%%%%%%%%%%%%%%%
for i=1:nconditions
    a=find(stimIds==i);
    respMatrix_condition=spike_times(a,:);
    raster_condition=nan(size(a,2),trialsize);
    spiketimes=respMatrix_condition;
    b=find(spiketimes==-1);
    spiketimes(b)=NaN;
    respMatrix_condition=respMatrix_condition./samplingFreq;
    spiketimes=spiketimes./samplingFreq;
    for j=1:size(a,2)
        cols=find(respMatrix_condition(j,:)>0);
        times=respMatrix_condition(j,cols);
        times=round(times*1000);
        raster_condition(j,times)=1*j;
    end
    spiketimes_condition=nanmean(spiketimes)';
    c=find(find(~isnan(spiketimes_condition)));
    spiketimes_condition=spiketimes_condition(c);
    spiketimes_condition=sort(spiketimes_condition,'ascend');
    assignin('base',strcat('spiketimes_condition_',int2str(i)),spiketimes_condition);
    assignin('base',strcat('raster_condition_',int2str(i)),raster_condition);
end