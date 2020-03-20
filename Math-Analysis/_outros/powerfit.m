function [a, b]=  powerfit(x,y,N)

% fits  y=a*x^b ;

% log(y) = log(a) + b*log(x);

ylog=log(y);

p = polyfit(log(x),log(y),N);

b = p(1);
a = exp(p(2));

end
