function fractalFano(data,direction)

name = eval(strcat('data.FF','_',direction,'_data.','name'));

nBins = eval(strcat('length(data.FF','_',direction,'_data.','size_of_bins)'));

info = eval(strcat('data.FF','_',direction,'_data.','results','_',direction));

if strcmp(direction,'across_bins')
    
    fanos = info.FF_per_trial_across_bins;
    
    nTrials = size(fanos,2);
    
    for n=1:nBins
        
        logFanoBins(n) = sum(fanos(n,1:nTrials))/nTrials;
        
    end
    
    logFanoBins = log10(logFanoBins);
    
    f = figure;
    
    plot(log10(1:nBins),logFanoBins,'r');
    print(f,'-depsc',strcat('/Volumes/Data/DATA/Surrogate/FF-Fractal/',name,'-FF-Fractal_across_bins-',direction,'.eps'));
    
    logFano = struct('name',name,'logFanoBins',logFanoBins);
    
    save(strcat('/Volumes/Data/DATA/Surrogate/FF-Fractal/',name,'-FF-Fractal_across_bins-',direction,'.mat'),'logFano');
    
elseif strcmp(direction,'across_trials')
    
    fanos = info.FF_across_trials_per_time_bin;
    
    for n=1:nBins
        
        logFanoTrials(n) = sum(fanos(n,1:n))/n;
        
    end
    
    logFanoTrials = log10(logFanoTrials);
    
    f = figure;
    
    plot(log10(1:nBins),logFanoTrials,'r');
    print(f,'-depsc',strcat('/Volumes/Data/DATA/Surrogate/FF-Fractal/',name,'-FF-Fractal_across_trials-',direction,'.eps'));
    
    logFano = struct('name',name,'logFanoTrials',logFanoTrials);
    
    save(strcat('/Volumes/Data/DATA/Surrogate/FF-Fractal/',name,'-FF-Fractal_across_trials-',direction,'.mat'),'logFano');
    
end

end

