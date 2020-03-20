function allSurrogatesPoisson

nTrials = 60;

time = 10000;

for r=2:6
    
    ratemax = r
    
    allTrials = spikeTrialsFromPoisson(ratemax,nTrials,time);
    
    save(strcat('/Volumes/Data/DATA/Surrogate/dados/','allTrials-surrogate-rate-',int2str(ratemax)),'allTrials');
    
end

end

