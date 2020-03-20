
function getReliability(spassDataSetName,samplingFrequency,startTime,endTime,binSize)

%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function loads a data set of electrophysiology coming from Spass
% software format and computes the Reliability parameter along trials
% for each condition following the definition as described on the reference 
% paper [1].
%
% INPUT:
%      1. spassDataSetName: the name of Spass data set file.
%      2. samplingFrenquency: the sampling frenquency used during
%      electrophysiological recordings.
%      3. startTime: the start time of the trial in milliseconds.
%      4. endTime: the end time of the trial in milliseconds.
%
% OUTPUT:
%      1.reliabilityPerCondition: structure with four variables, per condition and
% resolution:
%               a. Resolution: the resolution in seconds used to convolve
%               the spike train
%               b. reliability: Reliability parameter as described
%               c. all_trains_conv: the spike trains convolved
%               d. all_trains_delta: the spike trains converted to a binary
%               sequence
%
% [1] Schreiber, S., Fellous, J. M., Whitmer, D., Tiesinga, P., & Sejnowski, T. J. (2003). 
% A new correlation-based measure of spike timing reliability. Neurocomputing, 52-54, 925?931. 
%
% Author: Joao V. Dornas, joaodornas@gmail.com, 03/2017.
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

%%% SET RESOLUTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale = 10^3;
shift_cost = [0 2.^(-4:9)];
resolutions = ( 1./shift_cost );
resolutions(1) = [];
nResolutions = length(resolutions);  
sigma = resolutions.*scale;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DEFINE TOTAL TIME  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

totalTime = endTime - startTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% BEGIN CONDITIONS LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Begin Conditions Loop...');

for iResolution=1:nResolutions
    
    disp(strcat('Begin Resolution:',int2str(iResolution))); 
    
    for iCond=1:nConditions
       
       disp(strcat('Begin Condition:',int2str(iCond))); 
 
       trialsSpikes = spikeTimes(trialsLabels(iCond).label,:);
        
       nTrials = size(trialsSpikes,1);
       
       alltrials = trialsSpikes;

       alltrials = reshape(alltrials.',[],1);

       alltrials = alltrials.'; 
       
       latency_time_data = load(strcat(registro,'-latency.mat'));
       latencyTime = latency_time_data.latency_time;
       
       disp(strcat('nTrials:',int2str(nTrials)));
       
       all_trains_conv = zeros(nTrials,totalTime);
       all_trains_delta = zeros(nTrials,totalTime);
       
       for iTrial=1:nTrials
          
           [spike_train, delta_functions] = getSpikesConvolved(trialsSpikes(iTrial,:),sigma(iResolution),totalTime,startTime,scale,latencyTime);
           
           all_trains_conv(iTrial,:) = spike_train;
           all_trains_delta(iTrial,:) = delta_functions;

       end
       
       reliabilityPerCondition(iCond).resolution(iResolution).resolution = resolutions(iResolution);
       
       reliabilityPerCondition(iCond).resolution(iResolution).reliability = getReliabilityLocal(all_trains_conv);

       reliabilityPerCondition(iCond).resolution(iResolution).all_trains_conv = all_trains_conv;
                
       reliabilityPerCondition(iCond).resolution(iResolution).all_trains_delta = all_trains_delta;
       
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Save data...');
       
    save(strcat(registro,'-reliability.mat'),'reliabilityPerCondition','samplingFrequency','startTime','endTime','binSize');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc

end

function [spikes_convolved, delta_functions] = getSpikesConvolved(trial,sigma,totalTime,startTime,scale,latencyTime)

%%%   CONVOLVE SPIKE TRAIN WITH A GAUSSIAN FUNCTION   %%%%%%%%%%%%%%%%%%%%%

    disp('Convolve Spike Train...');
    
    trial = int32(trial*scale);
    trial = trial - startTime - latencyTime;

    trial(trial<=0) = [];
    trial(trial>totalTime) = [];

    delta_functions = zeros(1,totalTime);
    
    delta_functions(trial(:)) = 1;
    
    gauss_duration = 3*sigma;
    t = (-gauss_duration/2):1:(gauss_duration/2);
        
    gaussian = (1/(sqrt(2*pi*(sigma^2)))) * exp(-(t.^2)/(2*(sigma)^2));
    
    spikes_convolved = conv(delta_functions,gaussian);
    
    spikes_convolved = spikes_convolved./max(spikes_convolved);
    
    spikes_convolved = spikes_convolved(ceil(length(gaussian)/2):end-floor(length(gaussian)/2));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function reliability = getReliabilityLocal(all_trials)

%%%   CALCULATE RELIABILITY   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('Calculate Reliability...');
    
    reliability = 0;
    
    totalTrials = size(all_trials,1);
    
    for i=1:totalTrials
       
        for j=(i+1):totalTrials
        
            normai = norm(all_trials(i,:));
            normaj = norm(all_trials(j,:));
            
            reliability = reliability + dot(all_trials(i,:),all_trials(j,:))/(normai * normaj);
            
        end
        
    end

    reliability = reliability .* 2 ./ (totalTrials.*(totalTrials - 1));
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



