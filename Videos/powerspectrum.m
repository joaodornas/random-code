function [amp, f] = powerspectrum(picture, nbins)

    if nargin < 2
        nbins = 100;            % Default No of frequency 'bins'
    end
    if nargin < 3
        lowcut = 2;             % Ignore lower 2% of curve when finding
    end                         % line of best fit.

    picgray = imread(picture);
    
    [nrows ncolumns] = size(picgray);

    [x,y] = meshgrid(1:ncolumns,1:nrows);

    % The following fiddles the origin to the correct position
    % depending on whether we have and even or odd size.  
    % In addition the values of x and y are normalised to +- 0.5
    if mod(ncolumns,2) == 0
      x = (x-ncolumns/2-1)/ncolumns;
    else
      x = (x-(ncolumns+1)/2)/(ncolumns-1);
    end
    if mod(nrows,2) == 0
      y = (y-nrows/2-1)/nrows;
    else
      y = (y-(nrows+1)/2)/(nrows-1);
    end

    radius = sqrt(x.^2 + y.^2);
    
    % Quantise radius to the desired number of frequency bins
    radius = round(radius/max(max(radius))*(nbins-1)) + 1;
    
    amp = zeros(1, nbins);
    fcount = ones(1, nbins);
    
    fourier = fft2(double(picgray));
    amplitude = fftshift(abs(fourier));
    
    % Now go through the spectrum and build the histogram of amplitudes
    % vs frequences.
    for r = 1:nrows
        
        for c = 1:ncolumns
            
            ind = radius(r,c); 
            amp(ind) = amp(ind) + amplitude(r,c);
            fcount(ind) = fcount(ind) + 1;
            
        end
        
    end
    
    % Average the amplitude at each frequency bin. We also add 'eps'
    % to avoid potential problems later on in taking logs.
    amp = amp./fcount + eps;
    
    % Generate corrected frequency scale for each quantised frequency bin.
    % Note that the maximum frequency is sqrt(0.5) 
    f = [nbins:-1:1]/nbins*sqrt(.5);
    
%     figure(1);
%     plot(f,amp,'b');
%     title('Histogram of amplitude vs frequency');
%     
%     figure(2);
%     loglog(f,amp,'b-');
%     xlabel('log10. Spatial Frequency (cycles/image)');
%     ylabel('log10. Relative Amplitude');
%     title('Relative Amplitude x Spatial Frequency');
% 
%     figure(3);
%     plot(log10(f),log10(amp),'b--');
%     xlabel('log10. Spatial Frequency (cycles/image)');
%     ylabel('log10. Relative Amplitude');
%     title('Relative Amplitude x Spatial Frequency');
% 
%     figure(4);
%     plot(-log10(f),log(amp),'b-.');
%     xlabel('log10. Spatial Frequency (cycles/image)');
%     ylabel('log10. Power Spectrum');
%     title('Power Spectrum x Spatial Frequency');

    figure(5);
    plot(-2*log10(f),log10(amp.^2),'b-');
    xlabel('log10. Frequencia Espacial (ciclos/imagem)');
    ylabel('log10. Espectro de Potencia');
    title('LEI DE POTENCIA');

%     % Find line of best fit (ignoring specified fraction of low frequency values) 
%     fst = round(nbins*lowcut/100);
%     p = polyfit(-2*log10(f(fst:end)), log10(amp(fst:end).^2),1);
%     y = exp(p(1)*log(f) + p(2));        % log(y) = p(1)*log(f) + p(2)
%     hold on, loglog(f, y,'Color',[1 0 0]);
% 
%     n = round(nbins/10);
%     text(f(n), amp(n), sprintf('Slope = %f.2',p(1)));
%     title('Histogram of log amplitude vs log frequency');
% 
%     slope = p(1);

end