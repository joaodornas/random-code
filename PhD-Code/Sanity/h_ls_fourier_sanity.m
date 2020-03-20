function h_ls_fourier_sanity

san = load('sanity.mat');

dataLength = 1000;

Ci = san.Ci(1:dataLength);
Di = san.Di(1:dataLength);
Trans = san.Tevent;

%gc_windowing(Ci,Di,{'Ci','Di'},Trans);
%gc_windowing(Di,Ci,{'Di','Ci'},Trans);

Le = length(Ci);
fs = 1000;
f = (0:Le-1)*(fs/Le);
nW = length(f);
TO = 1;
Tmax = Le;

N = 2;
L = nW/tan((pi/2)*((N - 1)/N));

%causal = Ci./Di;
causal = Di./Ci;

causal(isnan(causal)) = 0;
causal(isinf(causal)) = 0;

TS = fft(causal);

F = zeros(length(TS),1);

F(1:nW) = exp(-1*1i.*(1:nW).*TO).*TS(1:nW);

for k=1:nW 
    
    [alphas, aInt, aR] = regressorRealF(F,nW,f,k,N);
    [betas, bInt, bR] = regressorImagF(F,nW,f,k,N);
   
    C(k) = getCausal(alphas,betas);
    
end

figure
plot(1:nW,C);


function causalStrength = getCausal(alphas,betas)

    n = length(alphas);
    
    numerator = 0;
    denominador = 0;
    
    for l=1:n
        
        numerator = numerator + (alphas(l) - betas(l))^2;

        denominador = denominador + (alphas(l)^2 + betas(l)^2);
        
    end
    
    causalStrength = 1 - sqrt(numerator)/sqrt(2*denominador);
    
end

function [alphas, aInt, aR] = regressorRealF(F,nW,f,T,N)
        
        F = F';
        x = real(F(1:nW).*exp(1i.*f(1:nW)*T));
        
        for i=1:N
            
           for j=1:nW
               
               y(i,j) = real(basis(f(j),i));
            
           end
           
        end
        
        x = x';
        y = y';
        
        [alphas, aInt, aR] = regress(x,y);
end

function [betas, bInt, bR] = regressorImagF(F,nW,f,T,N)
        
        F = F';
        x = -imag(F(1:nW).*exp(1i*f(1:nW)*T));
        
        for i=1:N
            
           for j=1:nW
               
               y(i,j) = imag(basis(f(j),i));
            
           end
           
        end
        
        x = x';
        y = y';
        
        [betas, bInt, bR] = regress(x,y);
end

function phi = basis(w,n)
        
    phi = ((1 + 1i*w)^n)/(1 - 1i*w)^(n+1);
        
end




end

