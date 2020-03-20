
N = 5000;

A = rand(N,N);

%%% TESTE 1

% matlabpool OPEN
% 
% for i=1:N
%     
%     X = A(:,i);
%    
%     gpu_X = gpuArray(X);
%     
%     tic
%     parfor j=1:N
%         
%         Y = A(:,j);
%     
%         gpu_Y = gpuArray(Y);
%         
%         rho = corr2(gpu_X,gpu_Y);
%         
%     end
%     toc
%     
% end
% 
% matlabpool CLOSE

%%% TESTE 2

matlabpool OPEN

gpu_A = gpuArray(A);

for i=1:N
    
    gpu_X = gpu_A(:,i);
   
    %gpu_X = gpuArray(X);
    
    tic
    parfor j=1:N
        
        gpu_Y = gpu_A(:,j);
    
        %gpu_Y = gpuArray(Y);
        
        rho = corr2(gpu_X,gpu_Y);
        
    end
    toc
    
end

matlabpool CLOSE