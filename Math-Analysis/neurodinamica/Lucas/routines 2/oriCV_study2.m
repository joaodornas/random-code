%double von mises:

%baseline+(pref_A*exp(pref_k*(cos(x-prefDir)-1)))+(anti_A*exp(anti_k*(cos(x-antiPrefDir)-1)))

%model parameters
baseline=0;
pref_A=10;
pref_k=2;
prefDir=deg2rad(90);
anti_A=0;
anti_k=0;
antiPrefDir=deg2rad(270);

%orientation circular variance (with integration)
syms x
a=int((baseline+(pref_A*exp(pref_k*(cos(x-prefDir)-1)))+(anti_A*exp(anti_k*(cos(x-antiPrefDir)-1)))).*exp(x*2*i),x,deg2rad(0),deg2rad(360));
b=int(baseline+(pref_A*exp(pref_k*(cos(x-prefDir)-1)))+(anti_A*exp(anti_k*(cos(x-antiPrefDir)-1))),x,deg2rad(0),deg2rad(360));
oriCV=1-(double(abs(a))/double(b));

%plot
angles=deg2rad(0:22.5:347.5);
response=baseline+(pref_A.*exp(pref_k.*(cos(angles-prefDir)-1)))+(anti_A.*exp(anti_k.*(cos(angles-antiPrefDir)-1)));
plot(angles,response);
