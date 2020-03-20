%model parameters
baseline=0;
pref_A=10;
pref_k=2.7;
prefDir=deg2rad(90);
anti_A=0;
anti_k=5;
antiPrefDir=deg2rad(270);

%orientation circular variance (with discrete sum)
angles=deg2rad(0:22.5:360);
response=baseline+(pref_A.*exp(pref_k.*(cos(angles-prefDir)-1)))+(anti_A.*exp(anti_k.*(cos(angles-antiPrefDir)-1)));

oriCV=1-(abs(sum(response.*exp(angles*2*i)))/(sum(response)));

plot(angles,response);