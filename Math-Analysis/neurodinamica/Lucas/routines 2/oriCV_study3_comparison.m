%model parameters
baseline=0;
pref_A=10;
%pref_k=2;
prefDir=deg2rad(90);
anti_A=0;
anti_k=0;
antiPrefDir=deg2rad(270);

oriCV_sum_vector=zeros(18,1);
oriCV_int_vector=zeros(18,1);

for n=1:18
    
    pref_k=0.1*(n+9);

    %orientation circular variance (with discrete sum)
    angles=deg2rad(0:22.5:360);
    response=baseline+(pref_A.*exp(pref_k.*(cos(angles-prefDir)-1)))+(anti_A.*exp(anti_k.*(cos(angles-antiPrefDir)-1)));
    oriCV_sum=1-(abs(sum(response.*exp(angles*2*i)))/(sum(response)));

    oriCV_sum_vector(n)=oriCV_sum;
    clear oriCV_sum

    %orientation circular variance (with integration)
    syms x
    a=int((baseline+(pref_A*exp(pref_k*(cos(x-prefDir)-1)))+(anti_A*exp(anti_k*(cos(x-antiPrefDir)-1)))).*exp(x*2*i),x,deg2rad(0),deg2rad(360));
    b=int(baseline+(pref_A*exp(pref_k*(cos(x-prefDir)-1)))+(anti_A*exp(anti_k*(cos(x-antiPrefDir)-1))),x,deg2rad(0),deg2rad(360));
    oriCV_int=1-(double(abs(a))/double(b));

    oriCV_int_vector(n)=oriCV_int;
    clear oriCV_int

end

plot(oriCV_int_vector,oriCV_sum_vector);
[r,p]=corr(oriCV_int_vector,oriCV_sum_vector)