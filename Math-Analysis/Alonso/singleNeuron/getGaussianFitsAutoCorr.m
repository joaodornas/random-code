function Fits_Gauss_Conditions = getGaussianFitsAutoCorr(Correlations_Conditions,kind)

    condicoes = max(size(Correlations_Conditions,1),size(Correlations_Conditions,2));
    
    for c=1:condicoes
         
        if strcmp(kind,'Spike')
            
            eval(strcat('atual_xcorr = Correlations_Conditions(c).xcorr_',kind,';'));
            eval(strcat('atual_correlogram = Correlations_Conditions(c).correlogram_',kind,';'));
            
        else
        
            eval(strcat('atual_xcorr = Correlations_Conditions(c).xcorr_',kind,';'));
            eval(strcat('atual_correlogram = Correlations_Conditions(c).correlogram_',kind,';'));
            
        end

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

        [sigmax, mux, Ax] =  mygaussfit(tracex,atual_xcorr,h);
        [sigma, mu, A] =  mygaussfit(tracecorr,atual_correlogram,h);

        
        eval(strcat('Fits_Gauss_Conditions(c).xcorr_',kind,'.sigmax = sigmax;'));

        eval(strcat('Fits_Gauss_Conditions(c).xcorr_',kind,'.mux = mux;'));

        eval(strcat('Fits_Gauss_Conditions(c).xcorr_',kind,'.Ax = Ax;'));
        
        eval(strcat('Fits_Gauss_Conditions(c).xcorr_',kind,'.xcorr_norm = atual_xcorr ./ Ax ;'));
        
        eval(strcat('Fits_Gauss_Conditions(c).correlogram_',kind,'.sigma = sigma;'));

        eval(strcat('Fits_Gauss_Conditions(c).correlogram_',kind,'.mu = mu;'));

        eval(strcat('Fits_Gauss_Conditions(c).correlogram_',kind,'.A = A;'));

        eval(strcat('Fits_Gauss_Conditions(c).correlogram_',kind,'.correlogram_norm = atual_correlogram ./ A;'));
        
        F_correlogram = @(y_correlogram) A * exp(( - (y_correlogram - mu).^2) / (2*sigma^2) ) ;
        
        F_xcorr = @(y_xcorr) Ax * exp(( - (y_xcorr - mux).^2) / (2*sigmax^2) ) ;
        
        area_correlogram = integral(F_correlogram,-length(atual_correlogram)/2 + 1,(length(atual_correlogram)/2));
        
        area_xcorr = integral(F_xcorr,-length(atual_xcorr)/2 + 1,(length(atual_xcorr)/2));
        
        eval(strcat('Fits_Gauss_Conditions(c).xcorr_',kind,'.area_xcorr = area_xcorr;'));
        
        eval(strcat('Fits_Gauss_Conditions(c).correlogram_',kind,'.area_correlogram = area_correlogram;'));
         
        eval(strcat('Fits_Gauss_Conditions(c).correlogram_',kind,'.latency_correlogram = find(Fits_Gauss_Conditions(c).correlogram_',kind,'.correlogram_norm == max(Fits_Gauss_Conditions(c).correlogram_',kind,'.correlogram_norm));'));
        
        eval(strcat('Fits_Gauss_Conditions(c).xcorr_',kind,'.latency_xcorr = find(Fits_Gauss_Conditions(c).xcorr_',kind,'.xcorr_norm == max(Fits_Gauss_Conditions(c).xcorr_',kind,'.xcorr_norm));'));
        

    end

end
