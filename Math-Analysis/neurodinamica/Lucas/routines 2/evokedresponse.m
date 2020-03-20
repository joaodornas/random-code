function [h,p]=evokedresponse(spike_times,stimIds,parameters)

%function [h,p]=evokedresponse(spike_times,stimIds,parameters)
%determines if evoked response is significantly greater than baseline for
%each condition, using wilcoxon's sign rank test for non-parametric
%distributions or paired t test for normal distributions, as determined by
%the lilliefors normality test
%returns two matrices, one containing boolean h, which is 0 if null
%hypothesis (h0) cannot be discarded and 1 if response is significant (h1)
%
%Inputs: spike_times is a standard matrix from Dinamo containing spike times
%in sampling units, where each row is a trial and each column is a spike time;
%stimIds is also a Dinamo matrix containing the stim label for each trial;
%parameters is a structure which should contain the fields:
%binsize: size of the psth bin in ms, default 10
%trialsize: duration of the trial in ms, default 7000
%nconditions: number of stimulation conditions, default 36
%samplingFreq: sampling frequency in Hz, default 32000
%baseline_duration: duration of pre-stimulation period in ms, default 1000
%alpha: desired level of significance, default 0.05
%by Lucas Pinto in March 2009

%get input arguments and/or set defaults
if nargin<2
    error('error: i need more input arguments!')
elseif nargin == 2
    binsize=10;
    trialsize=7000;
    nconditions=36;
    samplingFreq=32000;
    baseline_duration=1000;
    alpha=0.05;
elseif nargin == 3
    binsize=parameters.binsize;
    trialsize=parameters.trialsize;
    nconditions=parameters.nconditions;
    samplingFreq=parameters.samplingFreq;
    baseline_duration=parameters.baseline_duration;
    alpha=parameters.alpha;
end

%initialize matrices
h=nan(nconditions,1);
p=h;

%determine significance for each condition
for i=1:nconditions
    a=find(stimIds==i);
    respMatrix=spike_times(a,:)./samplingFreq;

    ntrials_i=size(respMatrix,1);
    meanActivity_matrix=zeros(ntrials_i,2);

    for l=1:ntrials_i

        %calculate psth for each trial
        psth_i=zeros(1,trialsize/binsize);

        for n=0:binsize:trialsize-1

            count=find(respMatrix(l,:)>n/1000 & respMatrix(l,:)<=(n+binsize)/1000);

            if isempty(count)==1
                psth_i((n/binsize)+1)=0;
            else psth_i((n/binsize)+1)=size(count,1)*(1000/binsize);
            end

        end

        %determine trial mean for baseline and evoked response of equal
        %time length
        meanActivity_matrix(l,:)=[mean(psth_i(1:baseline_duration/binsize)) mean(psth_i(baseline_duration/binsize+1:baseline_duration/binsize*2))];

    end
    
    if numel(meanActivity_matrix(:,1))<4 && numel(meanActivity_matrix(:,2))<4
        %sign rank if numel<4, in which case lillietest is not possible
        [p_i,h_i]=signrank(meanActivity_matrix(:,1),meanActivity_matrix(:,2),alpha);
    else
        %test for normality
        [h_lillie_base]=lillietest(meanActivity_matrix(:,1));
        [h_lillie_resp]=lillietest(meanActivity_matrix(:,2));

        %t test or sign rank depending on normality
        if h_lillie_base==0 && h_lillie_resp==0
            try
                [h_i,p_i]=ttest(meanActivity_matrix(:,1),meanActivity_matrix(:,2),alpha);
            catch
                [p_i,h_i]=signrank(meanActivity_matrix(:,1),meanActivity_matrix(:,2),alpha);
            end
        else
            [p_i,h_i]=signrank(meanActivity_matrix(:,1),meanActivity_matrix(:,2),alpha);
        end
    end
    
    h(i)=h_i;
    p(i)=p_i;

end