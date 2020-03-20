function kernel = randomkernel(n)


kernel = rand(n);

kernel = kernel * 1/(sum(sum(kernel)));


end