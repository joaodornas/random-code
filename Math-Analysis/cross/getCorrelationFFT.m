function corr = getCorrelationFFT(A,B)



    sd_A = std(A);
    
    sd_B = std(B);
    
%     M_A = mean(A);
%     
%     M_B = mean(B);
   
%     A = A - M_A;
%     
%     B = B - M_B;

    dim_A = max(size(A,1),size(A,2));
    
    dim_B = max(size(B,1),size(B,2));
    
    if dim_A > dim_B, B = [ B, zeros(1, (dim_A - dim_B) ) ]; end
    
    if dim_B > dim_A, A = [ A, zeros(1, (dim_B - dim_A) ) ]; end
    
    dim = max(dim_A,dim_B);
    
    nfftA = 2^nextpow2(2*dim-1);
    
    nfftB = 2^nextpow2(2*dim-1);
    
    corr = ifft( fft(A,nfftA) .* conj(fft(B,nfftB)) );
    
    corr = [corr(end-dim+2:end) , corr(1:dim)];
    
    corr = corr ./ (sd_A * sd_B * dim);


%AUTO-CORRELATION FFT

%x = rand(100,1);
%len = length(x);

%# autocorrelation
%nfft = 2^nextpow2(2*len-1);
%r = ifft( fft(x,nfft) .* conj(fft(x,nfft)) );

%# rearrange and keep values corresponding to lags: -(len-1):+(len-1)
%r = [r(end-len+2:end) ; r(1:len)];

%# compare with MATLAB's XCORR output
%all( (xcorr(x)-r) < 1e-10 )

end

