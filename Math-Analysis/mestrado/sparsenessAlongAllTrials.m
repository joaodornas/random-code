function sparsenessAlongAllTrials(date,site_index,channel,registro,video_index,start_time,end_time,bin_size,nConditions)

tic

disp('BEGIN');

%%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Spass...');

    Spass = load(strcat('_',registro,'-','v',int2str(video_index),'.mat'));

    
    allSpikeTrainsMatrix = Spass.spike_times ./ 32000;
    
    nTrials = size(allSpikeTrainsMatrix,1);
    
    for i=1:nTrials
        
        trains(i).vectors = allSpikeTrainsMatrix(i,:);
        
        trains(i).vectors = trains(i).vectors(trains(i).vectors>0);
        
        trains(i).vectors = trains(i).vectors(trains(i).vectors>=(start_time/1000) & trains(i).vectors<=(end_time/1000));
        
    end
    
        
    allSpikeTrainsVector = 0;
    
    for i=1:nTrials
        
        allSpikeTrainsVector = [allSpikeTrainsVector (trains(i).vectors + (10*(i-1)))];
        
    end
    
    allSpikeTrainsVector(1) = [];
    
    last_time = allSpikeTrainsVector(size(allSpikeTrainsVector,2));
    
     nBins = nTrials * (end_time - start_time) / (bin_size/1000);
       
    for k=1:nBins
          
        spikes = length(allSpikeTrainsVector(allSpikeTrainsVector>=((k-1)*bin_size/1000 + start_time/1000) & allSpikeTrainsVector<(k*bin_size/1000 + start_time/1000)));

        sparseData.binRate(k) = spikes/(bin_size/1000);

    end
    
     sparseData.SparsenessInTime = 0;
       
    for s=1:size(sparseData.binRate,2)
       
        sparseData.SparsenessInTime = [sparseData.SparsenessInTime ( 1 - ( mean(sparseData.binRate(1:s))^2 / ( mean(sparseData.binRate(1:s))^2 + std(sparseData.binRate(1:s))^2 ) ) / 1 - ( 1 / nBins ) )];
       
    end
    
    f = figure;

        filepath = strcat('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Centro-Contorno/Experimento - v1/Analise/',date,'/','sitio',int2str(site_index),'/',channel,'/','v',int2str(video_index));
        filepath = strcat(filepath,'/','sparseness','/','v',int2str(video_index),'-',date,'-','sitio',int2str(site_index),'-',channel);
        filepath = strcat(filepath,'-bin_size-',int2str(bin_size),'-sparseness-along-all-trials');
       
        save(filepath,'sparseData');
    
        plot(1:size(sparseData.SparsenessInTime,2),sparseData.SparsenessInTime)
        ylabel(strcat('Sparseness - Condition ',int2str(i)));
        xlabel('Time');

        print(f,'-djpeg',strcat(filepath,'-SparsenessAlongAllTials'));
        
    
end