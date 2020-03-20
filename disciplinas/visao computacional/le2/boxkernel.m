function kernel = boxkernel(n)


kernel = ones(n);

kernel = kernel * 1/(n*n);


end