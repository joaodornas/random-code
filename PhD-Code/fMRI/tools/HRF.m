

nTRs=170;

t=1:1:nTRs;
box=[ones(1,17),zeros(1,17),ones(1,17),zeros(1,17),ones(1,17),zeros(1,17),ones(1,17),zeros(1,17),ones(1,17),zeros(1,17)];

T0=0; n=4; lamda=0.5;
hrf=((t-T0).^(n-1)).*exp(-(t-T0)/lamda)/((lamda^n)*factorial(n-1));

B=conv(hrf,box)/10;

plot(zscore(ts));
hold on
plot(zscore(B(1:nTRs)),'r');

