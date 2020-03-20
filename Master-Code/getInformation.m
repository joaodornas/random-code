function getInformation(spassDataSetName,samplingFrequency,startTime,endTime,binSize)

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

allCondTrials = [];
for iCond=1:nConditions

   disp(strcat('Begin Condition:',int2str(iCond))); 

   trialsSpikes = spikeTimes(trialsLabels(iCond).label,:);

   nTrials = size(trialsSpikes,1);

   latency_time_data = load(strcat(registro,'-latency.mat'));
   latency_time = latency_time_data.latency_time;
   
   Cond(iCond).allTrials = [];

   for j=1:nTrials

       spike_train = trialsSpikes(j,:);

       a_e = spike_train(spike_train>=(0/1000) & spike_train<(latency_startTime/1000));

       after_ae = spike_train(spike_train>=(latency_startTime/1000) & spike_train<(endTime/1000)) - latency_time;

       spike_train = [];

       spike_train = [a_e, after_ae];

       spike_train = spike_train(spike_train>0);

       spike_train = spike_train(spike_train>=(startTime/1000) & spike_train<=(endTime/1000));

       spike_train = sort(spike_train);
       
       sizeTrials = size(Cond(iCond).allTrials,2);
       
       if (sizeTrials < length(spike_train))
           
           Cond(iCond).allTrials = [Cond(iCond).allTrials,zeros(size(Cond(iCond).allTrials,1),length(spike_train) - sizeTrials)];
           
       elseif (sizeTrials > length(spike_train))
           
           spike_train = [spike_train,zeros(1,sizeTrials - length(spike_train))];
           
       end
       
       Cond(iCond).allTrials = [Cond(iCond).allTrials;spike_train];

   end
   
   x(iCond) = size(Cond(iCond).allTrials,1);
   y(iCond) = size(Cond(iCond).allTrials,2);
   
end

[a, Ix] = sort(x,'descend');
[b, Iy] = sort(y,'descend');

allCondTrials = [];

for iCond=1:nConditions
    
    if b > size(Cond(iCond).allTrials,2)
    
        Cond(iCond).allTrials = [Cond(iCond).allTrials,zeros(size(Cond(iCond).allTrials,1),b-size(Cond(iCond).allTrials,2))];
    
    end
    
    if a > size(Cond(iCond).allTrials,1)
    
        Cond(iCond).allTrials = [Cond(iCond).allTrials;zeros(a-size(Cond(iCond).allTrials,1),size(Cond(iCond).allTrials,2))];
    
    end
        
    allCondTrials(:,:,iCond) = Cond(iCond).allTrials;
    
end
   
   R = allCondTrials;

   opts.nt = nTrials;
   opts.method = 'dr';
   opts.bias = 'pt';
   opts.btsp = 500;
   opts.verbose = false;

   disp('Begin Information Theory...');
   
   [X, Y] = information(R, opts, 'I', 'Ish');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Save data...');
       
    save(strcat(registro,'-information.mat'),'X','Y','opts','samplingFrequency','startTime','endTime','binSize');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

