function plotDongFractal(mainPath,bar,protocol,Afor,Aback,Pfor,Pback)

mkdir(mainPath,protocol);

plotFractalSpace(mainPath,protocol,Afor,Pfor,'Forward',bar);
plotFractalSpace(mainPath,protocol,Aback,Pback,'Backward',bar);

plotFractalTime(mainPath,protocol,Afor,Pfor,'Forward',bar);
plotFractalTime(mainPath,protocol,Aback,Pback,'Backward',bar);


function plotFractalSpace(mainPath,protocol,Amplitude,Power,kind,bar)


    f = 1:(size(Amplitude,1)/2);
    w = 1:(size(Amplitude,2)/2);

    LX = log10(f);
    LPX = log10(f);

    Y = Amplitude(f,1);
    PY = Power(f,1);
    for i=2:size(Amplitude,2)

        Y = Y + Amplitude(f,i);
        PY = PY + Power(f,i);

    end
    
    Y = Y ./ size(Amplitude,2);
    PY = PY ./ size(Amplitude,2);

    LY = log10(Y);
    LPY = log10(PY);

    [LX, LY] = removeInfZeros(LX,LY);    
    [LX, LY] = removeNegativos(LX,LY); 

    [LPX, LPY] = removeInfZeros(LPX,LPY);    
    [LPX, LPY] = removeNegativos(LPX,LPY); 

    if (length(LX) < 2) && (length(LY) < 2)

        coeff(1) = 0;
        coeff(2) = 0;

    else
        
        
        if ~iscolumn(LY), LY = LY.'; end
        
        if ~iscolumn(LX), LX = LX.'; end

        f = fit(LX,LY,'power1');       

        coeff = coeffvalues(f);

    end

    g = figure;

    set(gcf,'Visible','off');
    plot(LX,LY,'r');
    text(0.8*max(LX),0.8*max(LY),num2str(coeff(2)),'FontSize',20);
    print(g,'-depsc',strcat(mainPath,protocol,bar,protocol,'-Amplitude-Espectro-Video-across-space-',kind,'.eps'));

    if (length(LPX) < 2) && (length(LPY) < 2)

        coeff(1) = 0;
        coeff(2) = 0;

    else

        if ~iscolumn(LPY), LPY = LPY.'; end
        
        if ~iscolumn(LPX), LPX = LPX.'; end
        
        f = fit(2.*LPX,LPY,'power1');       

        coeff = coeffvalues(f);

    end

    h = figure;

    set(gcf,'Visible','off');
    plot(2.*LPX,LPY,'r');
    text(0.8*max(LPX),0.8*max(LPY),num2str(coeff(2)),'FontSize',20);
    print(h,'-depsc',strcat(mainPath,protocol,bar,protocol,'-Power-Espectro-Video-across-space-',kind,'.eps'));


end


function plotFractalTime(mainPath,protocol,Amplitude,Power,kind,bar)


    f = 1:(size(Amplitude,1)/2);
    w = 1:(size(Amplitude,2)/2);

    LX = log10(w);
    LPX = log10(w);


    Y = Amplitude(1,w);
    PY = Power(1,w);
    for i=2:size(Amplitude,1)

        Y = Y + Amplitude(i,w);
        PY = Y + Power(i,w);

    end
    
    Y = Y ./ size(Amplitude,1);
    PY = PY ./ size(Amplitude,1);

    LY = log10(Y);
    LPY = log10(PY);

    [LX, LY] = removeInfZeros(LX,LY);    
    [LX, LY] = removeNegativos(LX,LY); 

    [LPX, LPY] = removeInfZeros(LPX,LPY);    
    [LPX, LPY] = removeNegativos(LPX,LPY); 

    if (length(LX) < 2) && (length(LY) < 2)

        coeff(1) = 0;
        coeff(2) = 0;

    else
     
        if ~iscolumn(LY), LY = LY.'; end
        
        if ~iscolumn(LX), LX = LX.'; end
        
        f = fit(LX,LY,'power1');       

        coeff = coeffvalues(f);

    end

    g = figure;

    set(gcf,'Visible','off');
    plot(LX,LY,'r');
    text(0.8*max(LX),0.8*max(LY),num2str(coeff(2)),'FontSize',20);
    print(g,'-depsc',strcat(mainPath,protocol,bar,protocol,'-Amplitude-Espectro-Video-across-time-',kind,'.eps'));

    if (length(LPX) < 2) && (length(LPY) < 2)

        coeff(1) = 0;
        coeff(2) = 0;

    else

 
        if ~iscolumn(LPY), LPY = LPY.'; end
        
        if ~iscolumn(LPX), LPX = LPX.'; end
        
        f = fit(2.*LPX,LPY,'power1');       

        coeff = coeffvalues(f);

    end

    h = figure;

    set(gcf,'Visible','off');
    plot(2.*LPX,LPY,'r');
    text(0.8*max(LPX),0.8*max(LPY),num2str(coeff(2)),'FontSize',20);
    print(h,'-depsc',strcat(mainPath,protocol,bar,protocol,'-Power-Espectro-Video-across-time-',kind,'.eps'));


end


end

