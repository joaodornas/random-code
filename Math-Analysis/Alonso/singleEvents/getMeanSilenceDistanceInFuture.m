function [meanDistance, stdDistance] = getMeanSilenceDistanceInFuture(firstSpike,allTrials)


nTrials = max(size(allTrials,1),size(allTrials,2));

for n=1:nTrials
    
    spikes = allTrials(1,n).trial;
    
    spikes = spikes(spikes>firstSpike);
    
    spikes = sort(spikes);
    
    if numel(spikes) > 0
        
        spiketime(n) = spikes(1);
        
    else
        
        spiketime(n) = 0;
        
    end
    
end

spiketime(spiketime==0) = [];


for n=1:length(spiketime)
    
    allDistances(n) = abs(spiketime(n) - firstSpike);
    
end

meanDistance = mean(allDistances);
stdDistance = std(allDistances);
    
end


