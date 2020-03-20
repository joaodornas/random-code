function probability = F(t0, n, t )

tt = t - t0;

probability = 1 - exp(-n*tt);


end

