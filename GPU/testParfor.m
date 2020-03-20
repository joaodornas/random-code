
n = 40;

for i=1:n
    
    X(1:n,i) = i;

end

X;

matlabpool open

disp('GPU');

tic

A = gpuArray(X);

Y = zeros(n,n);

B = gpuArray(Y);

parfor i=1:n
    
    i
    
    for j=1:n
        
        B(i,j) = A(i,j);
        
    end
    
end

C = gather(B);

C;

toc

disp('no GPU');

tic

A = X;

Y = zeros(n,n);

B = gpuArray(Y);

for i=1:n
    
    for j=1:n
        
        B(i,j) = A(i,j);
        
    end
    
end

B;

toc

matlabpool close
    

