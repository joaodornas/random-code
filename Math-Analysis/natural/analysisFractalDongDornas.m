function analysisFractalDongDornas(protocol,bar,kind,mainPath)


Forwarddata = load(strcat(protocol,'-',kind,'Fourier-analysis-Forward.mat'));
   
Backwarddata = load(strcat(protocol,'-',kind,'Fourier-analysis-Backward.mat'));


Afor = abs(Forwarddata.data.RF);
Aback = abs(Backwarddata.data.RF);

Pfor = Afor.^2;
Pback = Aback.^2;

dim = size(Forwarddata);

if size(dim,2) == 2
    
    plotDongFractal(mainPath,bar,protocol,Afor,Aback,Pfor,Pback);
    
elseif size(dim,2) == 3
    
    plotDornasFractal(mainPath,bar,protocol,Afor,Aback,Pfor,Pback);
    
end



end

