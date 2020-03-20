function getinout = getinout(ruido,T,t)


in = [];
out = [];

for i=1:(t/T)
    
    in(i) = normrnd(0,1);
    out(i) = ( 1/4 + (3/4)*exp(-4*i*T) ) * in(i) + normrnd(0,0.02);    
    
end

getinout = [in;out];


end