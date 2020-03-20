function getKernel(spassDataSetName,samplingFrequency,start_time,end_time,binSize)

latency_start_time = 500;

latency_end_time = latency_start_time + 300;

p = 1000;

idx = strfind(spassDataSetName,'\');

registro = spassDataSetName(idx(length(idx))+1:length(spassDataSetName));

invert = 0;

%%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Spass...');

    Spass = load(char(spassDataSetName));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% LOAD CONDITIONS TRIALS LABELS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Conditions Trials Labels...');

nConditions = max(Spass.stimIds);

    for iCondition=1:nConditions
   
        trials_label(iCondition).labels = find(Spass.stimIds == iCondition); 
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%% SET RESOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Set Resolution...');

    spike_times = Spass.spike_times ./ samplingFrequency;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iCondition=1:(nConditions + 1)
    
    if iCondition == (nConditions + 1)
        
        j = 2;
        
    else
        
        j = iCondition;
        
    end
          
        trials_spikes = spike_times(trials_label(j).labels(:),:); 
        
        alltrials = trials_spikes;

        alltrials = reshape(alltrials.',[],1);

        alltrials = alltrials.'; 

        %latency_time_data = load(strcat(registro,'-latency.mat'));
        %latency_time = latency_time_data.latency_time;

        nTrials(j) = min(size(trials_spikes,1),size(trials_spikes,2));

        allSpikes = [];

        %if (j == 2) && (invert), disp('Begin invert spike train...'); end
        
        for k=1:nTrials(j)

            spikes = trials_spikes(k,:);

            %a_e = spikes(spikes>=(start_time/1000) & spikes<(latency_start_time/1000));

            % = spikes(spikes>=(latency_start_time/1000) & spikes<(end_time/1000)) - latency_time;

            %spikes = [];

            %spikes = [a_e, after_ae];

            spikes = spikes(spikes>0);

            spikes = sort(spikes);

            spikes = spikes(spikes>=(start_time/1000) & spikes<=(end_time/1000));

            if (j == 2) && (invert)

                bin_size = 1;

                nBins = (end_time - start_time)/bin_size;

                spikes_inverted = invertSpikeTrain(spikes,start_time,nBins,bin_size);
                spikes = spikes_inverted.spikes;

            end

            allSpikes = [allSpikes, spikes];

        end

        allSpikes = reshape(allSpikes.',[],1);

        allSpikes = sort(allSpikes);

        allSpikes = allSpikes(allSpikes>0);
        
        time_points = linspace(start_time/1000, end_time/1000, ( end_time/1000 - start_time/1000 ) * p);
        
        [kernel(iCondition).density, kernel(iCondition).timePoints, kernel(iCondition).optimalBinWidth, kernel(iCondition).WBinsTested, kernel(iCondition).Cost, kernel(iCondition).confb95] = ssvkernel(allSpikes,time_points);
    
end

datakernel = struct('kernel',kernel,'start_time',start_time,'end_time',end_time,'p',p,'latency_time',latency_time);

save(strcat(registro(1:length(registro)-3),'-kernel-density-function-',int2str(start_time),'-',int2str(end_time),'.mat'),'datakernel');

end

