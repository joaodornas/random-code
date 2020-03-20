function [ r,p ] = pearson( X,Y )


N = length(X);

A = 0;
for i=1:N
   
    A = A + ( X(i) - mean(X) )*( Y(i) - mean(Y) ) ;
    
end


B = 0;
for i=1:N
    
   B = B + ( X(i) - mean(X) )^2 ;
    
end

C = 0;
for i=1:N
    
   C = C + ( Y(i) - mean(Y) )^2 ;
    
end

r = A / ( sqrt(B) * sqrt(C) ) ;


t = r.*sqrt((N-2)./(1-r.^2)) ;

p = 2 * ( 1 - tcdf(abs(t),N-2) ) ;

end

