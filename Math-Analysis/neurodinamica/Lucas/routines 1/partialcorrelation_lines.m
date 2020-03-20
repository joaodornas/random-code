n=36;
p=0.1;
dof=n-3;
t=tinv(1-p,dof);
Rs_vector=zeros(1,201);
syms Rs
R=solve('(((1/2*log((1+Rs)/(1-Rs)))-Zi)/sqrt(2/(dof)))=t',Rs);

Zi=0;
Rs_line=eval(R);

for i=1:201
Ri=((i-1)/100)-1;
Zi=1/2*log((1+Ri)./(1-Ri));
Rs_vector(i)=eval(R);
end

a=find(Rs_vector<=Rs_line);
Rs_vector(a)=Rs_line;
Rs_vector(201)=1;

figure;
Ri_vector=-1:0.01:1;
plot(Ri_vector,Rs_vector,Rs_vector,Ri_vector,Rind,Rspeed)