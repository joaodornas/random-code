function [classif_vector,t_vector]=partcorr(Ri,Rs,n,p,lineplot)

%%partcorr calculates and plots partial correlation for column vector Ri
%and Rs, with "n" conditions and p value "p". Lineplot is the option for
%plotting, and should be set to 'on' or 'off'. Defaults are n=36, p=0.1,
%lineplot='on'
%[classif_vector,t_vector]=partcorr(Ri,Rs,n,p,lineplot)

% set variables
if nargin==2
    n=36;
    p=0.1;
    lineplot='on';
elseif nargin==3
    p=0.1;
    lineplot='on';
elseif nargin==4
    lineplot='on';
end
    
dof=n-3;
t_line=tinv(1-p,dof);
neg_t_line=-t_line;

ncells=size(Ri,1);
t_vector=zeros(size(Ri));

%calculate t values and classify
for i=1:ncells
    
    Zi=1/2*log((1+Ri(i))./(1-Ri(i)));
    Zs=1/2*log((1+Rs(i))./(1-Rs(i)));
    t=(Zs-Zi)/sqrt(2/dof);
    t_vector(i)=t;
    
    if t>=t_line
        if i==1
            classif_vector=cellstr('speed');
        else classif_vector=[classif_vector;'speed'];
        end
    elseif t<=neg_t_line
        if i==1
            classif_vector=cellstr('independent');
        else classif_vector=[classif_vector;'independent'];
        end
    else if i==1
            classif_vector=cellstr('intermediate');
        else classif_vector=[classif_vector;'intermediate'];
        end
    end

end

%plot
switch lineplot
    case 'on'
        Rs_vector=zeros(1,201);
        syms Rs_line
        R=solve('(((1/2*log((1+Rs_line)/(1-Rs_line)))-Zi_line)/sqrt(2/(dof)))=t_line',Rs_line);
        
        Zi_line=0;
        Rs_zero=eval(R);
        
        for i=1:201
            Ri_line=((i-1)/100)-1;
            Zi_line=1/2*log((1+Ri_line)./(1-Ri_line));
            Rs_vector(i)=eval(R);
        end

        a=find(Rs_vector<=Rs_zero);
        Rs_vector(a)=Rs_zero;
        Rs_vector(201)=1;
        Ri_vector=-1:0.01:1;
        
        figure;
        plot(Ri_vector,Rs_vector,Rs_vector,Ri_vector,Ri,Rs);
end
