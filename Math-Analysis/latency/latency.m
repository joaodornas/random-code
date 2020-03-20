function time = latency(trials,start_time,latency_end_time,p)

    trials = trials(trials>=0 & trials<=(latency_end_time/1000));
    
    time_points = linspace(1/1000, latency_end_time/1000, (latency_end_time/1000)*p);
    
    [kernel.density, kernel.timePoints, kernel.optimalBinWidth, kernel.WBinsTested, kernel.Cost, kernel.confb95] = ssvkernel(trials,time_points);
          
    timeHistoryBegin = (start_time/1000)*p + (20/1000)*p ;
    timeHistoryEnd = (start_time/1000)*p + (150/1000)*p + (20/1000)*p;
    
    Ypico = max(kernel.density(timeHistoryBegin:timeHistoryEnd));

    for k=timeHistoryBegin:timeHistoryEnd

        if kernel.density(k) == Ypico

            Xpico = k;

        end

    end

    ae = median(kernel.density(1:timeHistoryBegin));

    A = ( (Ypico - ae) / 2 ) + ae ;
    kernel.density = kernel.density';
    V = kernel.density(timeHistoryBegin:Xpico);
    XM = dsearchn(V,A);
    XM = XM + timeHistoryBegin;
    latencyPoints = kernel.timePoints(XM);
    
    time = (latencyPoints*1000) / p;
    
end
