function quociente = testeGPU

A = rand(3,3);

B = rand(3,3);

gpuA = gpuArray(A);

gpuB = gpuArray(B);
        
parfor i=1:5

    disp(int2str(i));
       
    p(i) = ( gpuA ./ gpuB ) .* i ;
        
    end
    
%end

quocience = gather(p);


end

