function E = gpuArrayFunTest(A,B)

%A = rand(N,N);

%B = rand(N,N);

C = A./B;

getSize = @size;

N = getSize(A,1);

D = zeros(N,N);

for i=1:N
   
    D(N,i) = A(i,N) + B(N,i);
    
end

E = D + C;

end

