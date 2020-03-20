function loadDataAndCompute(protocol,metric)

[pathToData,dataFolders] = loadPathToData(protocol);

samplingFrequency = 32000;
startTime = 500;
endTime = 9500;
binSize = 30;

for i=1:length(dataFolders)
    
    spassDataSetName = strcat(pathToData,dataFolders{i});

    switch metric
        
        case 'InformationToolbox'
            
            getInformation(spassDataSetName,samplingFrequency,startTime,endTime,binSize);
        
        case 'Kernel'
            
            getKernel(spassDataSetName,samplingFrequency,startTime,endTime,binSize);
        
        case 'Latency'
            
            getLatency(spassDataSetName,samplingFrequency,startTime,endTime,binSize);
    
        case 'Raster_Plot'
        
            getRasterPlot(spassDataSetName,samplingFrequency,startTime,endTime,binSize);
    
        case 'Fano_Factor'
        
            fanoFactorData = getfanoFactor(spassDataSetName,samplingFrequency,startTime,endTime,binSize);
        
        case 'Reliability'
        
            getReliability(spassDataSetName,samplingFrequency,startTime,endTime,binSize);
        
        case 'spass2Metric'
            
            switch protocol
    
                case 'VideosFB'
        
                    Conditions = 2;
                    
                    inverter_movie = 0;
            
                    spass2MetricSpace(spassDataSetName,Conditions,inverter_movie);
                    
                    inverter_movie = 1;
                    
                    spass2MetricSpace(spassDataSetName,Conditions,inverter_movie);
                    
                case 'VideosCC'
                    
                    Conditions = 2;
                    
                    inverter_movie = 0;
                    
                    spass2MetricSpace(spassDataSetName,Conditions,inverter_movie);
                    
            end
            
        case 'Metric_Space'
            
            nRepetitions = 60;
        
            getMetricSpace(spassDataSetName,startTime,endTime,nRepetitions);
        
        case 'Sparseness'
        
            getSparseness(spassDataSetName,samplingFrequency,startTime,endTime,binSize);
        
        case 'Direct_Method'
        
            
            getMutualInfo(spassDataSetName,samplingFrequency,startTime,endTime,binSize);
        
        case 'Power_Spectrum'
        
            powerSpectrumData = getPowerSpectrum(spassDataSetName,samplingFrequency,startTime,endTime,binSize);
        
        case 'Circular_Variance'
        
        circularVarianceData = getCircularVariance(spassDataSetName,samplingFrequency,startTime,endTime,binSize);
        
        case 'Spike_Count_Distributions'
        
            spikeCountDistributionsData = getSpikeCountDistributions(spassDataSetName,samplingFrequency,startTime,endTime,binSize);
        
        case 'InterSpike_Intervals'
        
            interSpikeIntervalsData = getInterSpikeIntervals(spassDataSetName,samplingFrequency,startTime,endTime,binSize);
        
    end
    
end
        



end

