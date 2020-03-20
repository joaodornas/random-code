function M = vonMise(m,A,teta,preferida,anti_preferida,DI)

k = 2;
k_2 = 2;

if (DI < 0.5) || (DI > 0.9)
    
    A_2 = 0;
    
else
    
    A_2 = 0.7 * A;
    
end
    
    M = m + ( A * exp( k * ( cos((teta - preferida)*((2*pi)/360)) - 1 ) ) ) + ( A_2 * exp( k_2 * ( cos((teta - anti_preferida)*((2*pi)/360)) - 1 ) ) ) ; 

end

