
howManyZeros = zeros(size(data,1),size(data,2),size(data,3));

run = 1;
for xdim=1:size(data,1)

    for ydim=1:size(data,2)

        for zdim=1:size(data,3)

            voxel_time_series = squeeze(data(xdim,ydim,zdim,:));
            
            howManyZeros(xdim,ydim,zdim) = length(find(voxel_time_series==0));
            
        end
        
    end
    
end
            
            