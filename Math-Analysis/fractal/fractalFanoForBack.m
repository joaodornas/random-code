function fractalFanoForBack(data,trials,direction)

name = eval(strcat('data.FF',trials,'_',direction,'.','name'));

nBins = eval(strcat('length(data.FF',trials,'_',direction,'.','size_of_bins)'));

info = eval(strcat('data.FF',trials,'_',direction,'.','results',trials,'_',direction));

if strcmp(direction,'across_bins')
    
    fanos = info.FF_per_trial_across_bins;
    
    nTrials = size(fanos,2);
    
    for n=1:nBins
        
        logFanoBins(n) = sum(fanos(n,1:nTrials))/nTrials;
        
    end
    
    logFanoBins = log10(logFanoBins);
    
    g = figure;
    
    X = log10(1:nBins);
    Y = logFanoBins;
    
    [X, Y] = removeInfZeros(X,Y);    
    [X, Y] = removeNegativos(X,Y);  

    if (length(X) < 2) && (length(Y) < 2)
        
        coeff(1) = 0;
        coeff(2) = 0;
        
    else
        
        f = fit(X.',Y.','power1');       
    
        coeff = coeffvalues(f);
    
    end
    
    plot(X,Y,'r');
    text(0.8*max(X),0.8*max(Y),num2str(coeff(2)),'FontSize',20);
    print(g,'-depsc',strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF-Fractal/',name,'-FF-Fractal_across_bins-',trials,'.eps'));
    
    logFano = struct('name',name,'logFanoBins',logFanoBins,'coeff',coeff);
    
    save(strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF-Fractal/',name,'-FF-Fractal_across_bins-',trials,'.mat'),'logFano');
    
elseif strcmp(direction,'across_trials')
    
    fanos = info.FF_across_trials_per_time_bin;
    
    for n=1:nBins
        
        logFanoTrials(n) = sum(fanos(n,1:n))/n;
        
    end
    
    logFanoTrials = log10(logFanoTrials);
    
    g = figure;
    
    X = log10(1:nBins);
    Y = logFanoTrials;
    
    [X, Y] = removeInfZeros(X,Y);   
    [X, Y] = removeNegativos(X,Y);  
    
     if (length(X) < 2) && (length(Y) < 2)
        
        coeff(1) = 'NaN';
        coeff(2) = 'NaN';
        
     else
        
        f = fit(X.',Y.','power1');
    
        coeff = coeffvalues(f);
    
     end
     
    plot(X,Y,'r');
    text(0.8*max(X),0.8*max(Y),num2str(coeff(2)),'FontSize',20);
    print(g,'-depsc',strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF-Fractal/',name,'-FF-Fractal_across_trials-',trials,'.eps'));
    
    logFano = struct('name',name,'logFanoTrials',logFanoTrials,'coeff',coeff);
    
    save(strcat('/Volumes/Data/DATA/Forward-Backward/Power-Law/Fano-Factor/FF-Fractal/',name,'-FF-Fractal_across_trials-',trials,'.mat'),'logFano');
    
end

end

