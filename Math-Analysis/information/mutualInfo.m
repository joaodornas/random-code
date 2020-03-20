
function mutualInfo(protocols,start_stimulus_time,end_stimulus_time,folder)

%% BIN SIZE

sizes_of_a_bin = 1:(end_stimulus_time-start_stimulus_time);

nBinsSizes = length(sizes_of_a_bin);

%% LATENCY TIME

latency_start_time = start_stimulus_time;

latency_end_time = latency_start_time + 300;

p = 10000;

%% PROTOCOLS

if strcmp(protocols(end-2:end),'mat')
    
    registros = protocols;
    nReg = 1;
    
elseif strcmp(protocols(end-2:end),'txt')
    
    registros = importdata(protocols);
    nReg = length(registros);

end

%% EACH LOOP FOR A PROTOCOL

matlabpool('AttachedFiles',{'gpuGetWordMatrix.m','gpuGetWordFrequencies.m','gpuGetCondEntropy.m','count_unique.m'});

for r=1:nReg
    
    if nReg == 1
        
        protocol = registros;
        name = registros(2:end-4);
        
    else
        
        protocol = registros{r};
        name = registros{r}(2:end-4);
        
    end
    
    disp(strcat('Protocol:',name));
    
    Spass = load(char(protocol));
    
    nConditions = Spass.parameters.nconditions;
    
    for i=1:nConditions
   
        trials_label(i).labels = find(Spass.stimIds == i); 
    
    end
    
    spike_times = Spass.spike_times ./ 32000;
    
    %% EACH LOOP FOR A CONDITION
    
    for i=1:nConditions
        
        disp(strcat('Condition:',int2str(i)));
        
        all_trials = spike_times(trials_label(i).labels(:),:);
        
        nTrials = min(size(all_trials,1),size(all_trials,2));
       
        latency_time = getLatency(name,i);
        
        for WO=0:1
        
            disp(strcat('WOLatency:',int2str(WO)));
            
            spike_trains = takeLatencyOff(all_trials,nTrials,start_stimulus_time,end_stimulus_time,latency_start_time,latency_time*WO);
            
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
            
            parfor bin=1:nBinsSizes

                disp(strcat('bin_size:',int2str(sizes_of_a_bin(bin))));

                nBins = ceil( ( end_stimulus_time - start_stimulus_time ) / sizes_of_a_bin(bin) );
                
                disp(strcat('getWordMatrix'));
                
                WordMatrix(bin).wm = gpuGetWordMatrix(spike_trains,nTrials,nBins,sizes_of_a_bin(bin),start_stimulus_time);
                
                disp(strcat('getWordFrequencies'));
                
                [WordsUnique(bin).w, nWordsUnique(bin).nw, Frequencies(bin).f, CondWordsUnique(bin).cw, nCondWordsUnique(bin).cnw, CondFrequencies(bin).cf] = gpuGetWordFrequencies(WordMatrix(bin).wm);

                disp('getEntropy');

                [HResponse(bin), HResponseMAbound(bin)] = gpuGetEntropy(Frequencies(bin).f);
                
                disp('getCondEntropy');
                
                [HCondResponse(bin).hcr, HCondResponseMAbound(bin).hcr] = gpuGetCondEntropy(CondFrequencies(bin).cf); 
                
                meanHCondResponse(bin) = mean(HCondResponse(bin).hcr);
                stdHCondResponse(bin) = std(HCondResponse(bin).hcr);

                meanHCondResponseMAbound(bin) = mean(HCondResponseMAbound(bin).hcr);
                stdHCondResponseMAbound(bin) = std(HCondResponseMAbound(bin).hcr);

                MutualInfo(bin) = HResponse(bin) - meanHCondResponse(bin);
                CodingEfficiency(bin) = ( ( HResponse(bin) - meanHCondResponse(bin) ) / HResponse(bin) ) * 100;

            end
            
                info.Condition(i).Latency(WO+1).WordMatrix = WordMatrix;
                
                info.Condition(i).Latency(WO+1).WordsUnique = WordsUnique;
                info.Condition(i).Latency(WO+1).nWordsUnique = nWordsUnique;
                info.Condition(i).Latency(WO+1).Frequencies = Frequencies;
                info.Condition(i).Latency(WO+1).CondWordsUnique = CondWordsUnique;
                info.Condition(i).Latency(WO+1).nCondWordsUnique = nCondWordsUnique;
                info.Condition(i).Latency(WO+1).CondFrequencies = CondFrequencies;
                
                info.Condition(i).Latency(WO+1).HResponse = HResponse;
                info.Condition(i).Latency(WO+1).HResponseMAbound = HResponseMAbound;
                
                info.Condition(i).Latency(WO+1).HCondResponse = HCondResponse;
                info.Condition(i).Latency(WO+1).HCondResponseMAbound = HCondResponseMAbound;
                
                info.Condition(i).Latency(WO+1).meanHCondResponse = meanHCondResponse;
                info.Condition(i).Latency(WO+1).stdHCondResponse = stdHCondResponse;

                info.Condition(i).Latency(WO+1).meanHCondResponseMAbound = meanHCondResponseMAbound;
                info.Condition(i).Latency(WO+1).stdHCondResponseMAbound = stdHCondResponseMAbound;

                info.Condition(i).Latency(WO+1).MutualInfo = MutualInfo;
                info.Condition(i).Latency(WO+1).CodingEfficiency = CodingEfficiency;
            
            for b=1:length(sizes_of_a_bin)

                MutualInfo(b) = info.Condition(i).Latency(WO+1).MutualInfo(b);
                CodingEfficiency(b) = info.Condition(i).Latency(WO+1).CodingEfficiency(b);

            end
            
            info.Condition(i).Latency(WO+1).MaxMutualInfoBits = max(MutualInfo);

            info.Condition(i).Latency(WO+1).MaxMutualInfoResolution = find(MutualInfo==max(MutualInfo));

            info.Condition(i).Latency(WO+1).MaxMutualInfoCodingEfficiency = CodingEfficiency(find(MutualInfo==max(MutualInfo)));

            info.Condition(i).Latency(WO+1).MaxCodingEfficiency = max(CodingEfficiency);

            info.Condition(i).Latency(WO+1).MaxCodingEfficiencyResolution = find(CodingEfficiency==max(CodingEfficiency));
            
            
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

    end
    
    info.sizes_of_a_bin = sizes_of_a_bin;
        
    save(strcat(folder,name,'-mutual-information'),'info');
    
end

matlabpool close;
        
end


