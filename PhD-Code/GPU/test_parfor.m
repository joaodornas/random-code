
N = 300;
M = 5000;
T1 = M;
T2 = M;

A = rand(N,M);
B = rand(N,M);

matlabpool OPEN

%% CPU
tic
t = [];
for i=1:T1
    
    X = A(:,i);
    
    for j=1:T2
        
        Y = B(:,j);
        
        %tic;
        rho = corr2(X,Y);
        %a_std = std(A);
        %t(j) = toc;
        
    end
    
end
t = toc;

disp(strcat('Mean toc For Loop:',num2str(mean(t))));

%% CPU - PARALLEL
tic
t = [];
for i=1:T1
    
    X = A(:,i);
    %Y = X.*2;
    
    parfor j=1:T2
        
        Y = B(:,j);
        
        %tic;
        rho = corr2(X,Y);
        %a_std = std(A);
        %t(j) = toc;
        
    end
    
end
t = toc;

disp(strcat('Mean toc Par For Loop:',num2str(mean(t))));

matlabpool CLOSE