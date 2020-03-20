
N = 300;
T1 = 1;
T2 = 100000;

A = rand(N,1);
B = rand(N,1);

gpu_A = gpuArray(A);
gpu_B = gpuArray(B);

matlabpool OPEN

%% CPU
tic
for i=1:T1
    
    for j=1:T2
        
        rho = corr2(A,B);
        %a_std = std(A);
        
    end
    
end
toc

%% CPU - PARALLEL
tic
for i=1:T1
    
    parfor j=1:T2
        
        rho = corr2(A,B);
        %a_std = std(A);
        
    end
    
end
toc

%% GPU FOR
tic
for i=1:T1
    
    for j=1:T2
        
        rho = corr2(gpu_A,gpu_B);
        %a_std = std(gpu_A);
        
    end
    
end
toc

%% GPU PAR
tic
for i=1:T1
    
    parfor j=1:T2
        
        rho = corr2(gpu_A,gpu_B);
        %a_std = std(gpu_A);
        
    end
    
end
toc

matlabpool CLOSE