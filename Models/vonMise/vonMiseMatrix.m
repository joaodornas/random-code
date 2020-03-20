function M = vonMiseMatrix(cdata,A,preferida,anti_preferida,DI)

MM = size(cdata,1);
NN = size(cdata,2);

M = 0;

k = 2;
k_2 = 2;

if (DI < 0.5) || (DI > 0.9)
    
    A_2 = 0;
    
else
    
    A_2 = 0.7 * A;
    
end

for mm=1:MM

    for nn=1:NN

        teta = cdata(mm,nn);

        M = M + ( A * exp( k * ( cos(teta - ((preferida)*((2*pi)/360))) - 1 ) ) ) + ( A_2 * exp( k_2 * ( cos(teta - ((anti_preferida)*((2*pi)/360))) - 1 ) ) ) ; 

    end

end

M = M / ( MM * NN );

end

