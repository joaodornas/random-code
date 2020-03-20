function Fits_Gauss_Conditions = getGaussianFitsAutoCorrEvents(Correlations_Conditions,kind)

    condicoes = max(size(Correlations_Conditions,1),size(Correlations_Conditions,2));
    
    for c=1:condicoes
         
        nEvents = max(size(Correlations_Conditions(c).eventos,1),size(Correlations_Conditions(c).eventos,2));
        
        if nEvents >= 1

            for n=1:nEvents


                eval(strcat('atual_xcorr = Correlations_Conditions(c).eventos(n).xcorr_',kind,';'));
                eval(strcat('atual_correlogram = Correlations_Conditions(c).eventos(n).correlogram_',kind,';'));


                if mod(length(atual_xcorr),2) ~= 0

                        k = 1;
                        j = 0;

                else

                        k = 0;
                        j = 1;

                end

                if mod(length(atual_correlogram),2) ~= 0

                        w = 1;
                        u = 0;

                else

                        w = 0;
                        u = 1;


                end

                tracex = (-( length(atual_xcorr) - k )/2 + j):( ( length(atual_xcorr) - k )/2);
                tracecorr = (- ( length(atual_correlogram) - w )/2 + u):( ( length(atual_correlogram) - w )/2);

                h = 0.2;

                trace = min(tracex,tracecorr);

                [sigmax, mux, Ax] =  mygaussfit(trace,atual_xcorr,h);
                [sigma, mu, A] =  mygaussfit(trace,atual_correlogram,h);


                eval(strcat('Fits_Gauss_Conditions(c).eventos(n).xcorr_',kind,'.sigmax = sigmax;'));

                eval(strcat('Fits_Gauss_Conditions(c).eventos(n).xcorr_',kind,'.mux = mux;'));

                eval(strcat('Fits_Gauss_Conditions(c).eventos(n).xcorr_',kind,'.Ax = Ax;'));

                eval(strcat('Fits_Gauss_Conditions(c).eventos(n).xcorr_',kind,'.xcorr_norm = atual_xcorr ./ Ax ;'));

                eval(strcat('Fits_Gauss_Conditions(c).eventos(n).correlogram_',kind,'.sigma = sigma;'));

                eval(strcat('Fits_Gauss_Conditions(c).eventos(n).correlogram_',kind,'.mu = mu;'));

                eval(strcat('Fits_Gauss_Conditions(c).eventos(n).correlogram_',kind,'.A = A;'));

                eval(strcat('Fits_Gauss_Conditions(c).eventos(n).correlogram_',kind,'.correlogram_norm = atual_correlogram ./ A;'));

                F_correlogram = @(y_correlogram) A * exp(( - (y_correlogram - mu).^2) / (2*sigma^2) ) ;

                F_xcorr = @(y_xcorr) Ax * exp(( - (y_xcorr - mux).^2) / (2*sigmax^2) ) ;

                area_correlogram = integral(F_correlogram,min(trace),max(trace));

                area_xcorr = integral(F_xcorr,min(trace),max(trace));

                eval(strcat('Fits_Gauss_Conditions(c).eventos(n).xcorr_',kind,'.area_xcorr = area_xcorr;'));

                eval(strcat('Fits_Gauss_Conditions(c).eventos(n).correlogram_',kind,'.area_correlogram = area_correlogram;'));

                eval(strcat('Fits_Gauss_Conditions(c).eventos(n).correlogram_',kind,'.latency_correlogram = find(Fits_Gauss_Conditions(c).eventos(n).correlogram_',kind,'.correlogram_norm == max(Fits_Gauss_Conditions(c).eventos(n).correlogram_',kind,'.correlogram_norm));'));

                eval(strcat('Fits_Gauss_Conditions(c).eventos(n).xcorr_',kind,'.latency_xcorr = find(Fits_Gauss_Conditions(c).eventos(n).xcorr_',kind,'.xcorr_norm == max(Fits_Gauss_Conditions(c).eventos(n).xcorr_',kind,'.xcorr_norm));'));

            end

        else
            
            Fits_Gauss_Conditions(c).eventos = [];
            
        end
        
    end

end
