function [psth]=makepsth(spike_times,stimIds,parameters,plotPsth)

%function [psth]=makepsth(spike_times,stimIds,parameters,plotPsth)
%returns a structure containing the psth for each stimulation condition
%
%Inputs: spike_times is a standard matrix from Dinamo containing spike times
%in sampling units, where each row is a trial and each column is a spike; 
%stimIds is also a Dinamo matrix containing the stim label for each trial;
%parameters is a structure which should contain the fields:
%binsize: size of the psth bin in ms, default 10
%trialsize: duration of the trial in ms, default 7000
%nconditions: number of stimulation conditions, default 36
%samplingFreq: sampling frequency in Hz, default 32000
%
%plotPsth is a string. 'on' if user wants to plot psths, 'off' otherwise
%by Lucas Pinto in March 2009

%get input arguments and/or set defaults
if nargin<2
    'error: i need more input arguments!'
elseif nargin == 2
    binsize=10;
    trialsize=7000;
    nconditions=36;
    samplingFreq=32000;
    plotPsth='off';
elseif nargin == 3
    binsize=parameters.binsize;
    trialsize=parameters.trialsize;
    nconditions=parameters.nconditions;
    samplingFreq=parameters.samplingFreq;
    plotPsth='off';
elseif nargin == 4
    binsize=parameters.binsize;
    trialsize=parameters.trialsize;
    nconditions=parameters.nconditions;
    samplingFreq=parameters.samplingFreq;
end

%psth for each condition
for i=1:nconditions
    a=find(stimIds==i);
    respMatrix=spike_times(a,:)./samplingFreq;
    psth_i=zeros(1,trialsize/binsize);

    ntrials_i=size(respMatrix,1);

    for n=0:binsize:trialsize-1

        count=find(respMatrix>n/1000 & respMatrix<=(n+binsize)/1000);

        if isempty(count)==1
            psth_i((n/binsize)+1)=0;
        else psth_i((n/binsize)+1)=(size(count,1)/ntrials_i)*(1000/binsize);
        end

    end

    %include in psth data structure
    eval(['psth.condition' num2str(i) '=psth_i;']);

end

%plot if plotPsth='on'
switch plotPsth
    case 'on'
        figure;
        for p=1:nconditions
            subplot(6,6,p);
            plot(eval(['psth.condition' num2str(p)]));
        end
end