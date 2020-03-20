function [spikes_convolved, delta_functions] = getSpikesConvolved(trial,sigma,total_time,start_time,end_time)
    
    trial = int32(trial*1000);
    trial = trial - start_time;

    trial(trial<=0) = [];
    trial(trial>total_time) = [];

    delta_functions = zeros(1,total_time);
    
    delta_functions(trial(:)) = 1;
    
    gauss_duration = 3*sigma;
    t = (-gauss_duration/2):1:(gauss_duration/2);
        
    gaussian = (1/(sqrt(2*pi*(sigma^2)))) * exp(-(t.^2)/(2*(sigma)^2));
    
    spikes_convolved = conv(delta_functions,gaussian);
    
    spikes_convolved = spikes_convolved./max(spikes_convolved);
    
    spikes_convolved = spikes_convolved(ceil(length(gaussian)/2):end-floor(length(gaussian)/2));

end

