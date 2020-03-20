function allTrials = spikeTrialsFromPoisson(ratemax,nTrials,time)

    one = 1;
    
    refractory_time = 2;

    allTrials = zeros(nTrials,time);
    
    rate = poisspdf(one,ratemax);

    for n=1:nTrials
        
        allTrials(n,:) = getSpikeTrain(rate,time,refractory_time);
        
    end
    
    
    
function spike_train = getSpikeTrain(rate,time,refractory_time)    
    
    spike_train = zeros(1,time);

    for j=1:time

        if (j == 1) || (j == 2)

           M = 0;

        else

           M = refractory_time;

        end

        probability = rand(1); 

        if (probability < rate) && (spike_train(j-M) == 0)

           spike_train(j) = j/1000; 

        end

    end
    
end
         
end