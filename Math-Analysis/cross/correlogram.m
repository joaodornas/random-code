function correlation = correlogram(A,B,kind)

bin_size = 1;

min_A = min(A);
max_A = max(A);

min_B = min(B);
max_B = max(B);

min_A = 0;
min_B = 0;

if strcmp(kind,'spike_train')

    nBins_A = round( ( (max_A - min_A) / bin_size ) * 1000 );
    nBins_B = round( ( (max_B - min_B) / bin_size ) * 1000 );

elseif strcmp(kind,'time_vector')

    nBins_A = (max_A - min_A);
    nBins_B = (max_B - min_B);

end

size = max([nBins_A nBins_B]);

time_vector_A = zeros(1,size);
time_vector_B = zeros(1,size);

A(A==0) = [];
B(B==0) = [];

for i=1:length(A)
   
    if isnan(A(i)), A(i) = 0; end
    if isinf(A(i)), A(i) = 0; end
   
    if strcmp(kind,'spike_train')

        spike_time = round( ( A(i) / bin_size ) * 1000 );

    elseif strcmp(kind,'time_vector');

        spike_time = A(i);

    end
    
    if spike_time ~= 0
        
        time_vector_A(uint16(spike_time)) = 1;
        
    end
    
end

for i=1:length(B)
    
    if isnan(B(i)), B(i) = 0; end
    if isinf(B(i)), B(i) = 0; end
    
    if strcmp(kind,'spike_train')

        spike_time = round( ( B(i) / bin_size ) * 1000 );

    elseif strcmp(kind,'time_vector');

        spike_time = B(i);

    end
    
    if spike_time ~= 0
        
        time_vector_B(uint16(spike_time)) = 1;
        
    end
    
end

correlation = zeros(1,2*size + 1);

for i=1:length(time_vector_A)

    if time_vector_A(i) == 1
        
        if i ~= 1
            
            for j=(size + 2 - i):(size)
            
                w = i - (size + 1 - j) + 1; 
            
                correlation(j) = correlation(j) + time_vector_B(w);
            
            end
            
        elseif i == 1
            
            for j=(size + 2 - i):(size)
            
                w = i - (size + 1 - j); 
            
                correlation(j) = correlation(j) + time_vector_B(w);
            
            end
            
        end
            
        
        if i ~= size
            
            for j=(size + 2):(2*size + 1 - i)
            
                w = i + (j - (size + 1)) - 1;
            
                correlation(j) = correlation(j) + time_vector_B(w);
            
            end
            
        elseif i == size
            
              for j=(size + 2):(2*size + 1 - i)
            
                w = i + (j - (size + 1));
            
                correlation(j) = correlation(j) + time_vector_B(w);
            
              end
            
        end
        
        correlation(size + 1) = correlation(size + 1) + time_vector_B(i);
               
    end
    
end

end

