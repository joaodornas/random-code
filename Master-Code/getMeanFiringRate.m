function getMeanFiringRate(spassDataSetName,samplingFrequency,start_time,end_time,binSize)

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

    for i=1:nConditions
   
        trials_label(i).label = find(Spass.stimIds == i); 
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%% SET RESOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Set Resolution...');

    spike_times = Spass.spike_times ./ samplingFrequency;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% BEGIN CONDITIONS LOOP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Begin Conditions Loop...');

nConditions = max(Spass.stimIds);

    for i=1:nConditions
       
       disp(strcat('Begin Condition . ',int2str(i))); 
       
       trials_spikes = spike_times(trials_label(i).label,:);
        
       nTrials = size(trials_spikes,1);
       
       spikes_vector = reshape(trials_spikes.',[],1);
       
       spikes_vector = spikes_vector.';
       
       spikes_vector = sort(spikes_vector);
       
       spikes_vector = spikes_vector(spikes_vector>0);

       spikes_vector = spikes_vector(spikes_vector>(start_time/1000) & spikes_vector<(end_time/1000));
       
       meanFR(i) = ( length(spikes_vector) / ( nTrials ) ) / ( (end_time/1000) - (start_time/1000) ) ;
        
       nBins = ( (end_time/1000) - (start_time/1000) ) / (binSize/1000);
       
        for k=1:nBins
          
                spikes = length(spikes_vector(spikes_vector>=((k-1)*binSize/1000 + start_time/1000) & spikes_vector<(k*binSize/1000 + start_time/1000)));
           
                binRate(k) = spikes/((nTrials)*(binSize/1000));
                
                condition(i).binRate(k) = binRate(k);

        end
       
       minFR(i) = min(binRate);
       maxFR(i) = max(binRate);
       desvio(i) = std(binRate);
      
       
%        capable = @(x) foomean(x);
%        bootStrap(i).ci = bootci(10000,{capable,binRate},'type','bca');
       
       if (i == 2) && (invert)

            disp('...invert spike train .');
            
            binSize_resolution = 1;
            
            nBins_new = (end_time - start_time)/binSize_resolution;
            spikes_inverted = invertSpikeTrain(spikes_vector,start_time,nBins_new,binSize_resolution);
            spikes_vector_inverted = [];
            spikes_vector_inverted = spikes_inverted.spikes;
            
            spikes_vector_inverted = sort(spikes_vector_inverted);
       
            spikes_vector_inverted = spikes_vector_inverted(spikes_vector_inverted>0);

            spikes_vector_inverted = spikes_vector_inverted(spikes_vector_inverted>(start_time/1000) & spikes_vector_inverted<(end_time/1000));
       
            meanFR(4) = ( length(spikes_vector_inverted) / ( nTrials ) ) / ( (end_time/1000) - (start_time/1000) ) ;
            
            nBins = ( (end_time/1000) - (start_time/1000) ) / (binSize/1000);
       
            binRate = [];
            
            for k=1:nBins
          
                spikes = length(spikes_vector_inverted(spikes_vector_inverted>=((k-1)*binSize/1000 + start_time/1000) & spikes_vector_inverted<(k*binSize/1000 + start_time/1000)));
           
                binRate(k) = spikes/((nTrials)*(binSize/1000));
                
                condition(4).binRate(k) = binRate(k);

            end
            
            minFR(4) = min(binRate);
            maxFR(4) = max(binRate);
            desvio(4) = std(binRate);
       
%             capable = @(x) foomean(x);
%             bootStrap(i).ci = bootci(10000,{capable,binRate},'type','bca');
       
       end
       
    end

    
save(strcat(registro(1:length(registro)-3),'-mean-firing-rate-',int2str(start_time),'-',int2str(end_time),'.mat'));


end
