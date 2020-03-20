
function output = mean_4D(input)

howManyZeros = zeros(size(input,1),size(input,2),size(input,3));

nTR = size(input,4);

all_data_points = [];

for xdim=1:size(input,1)

    for ydim=1:size(input,2)

        for zdim=1:size(input,3)

            voxel_time_series = squeeze(input(xdim,ydim,zdim,:));
            
            howManyZeros(xdim,ydim,zdim) = length(find(voxel_time_series>8));
            
        end
        
    end
    
end

index = find(howManyZeros>0);

%tic
output = zeros(length(nTR),1);
for iTR=1:nTR
    
    all_voxels = [];
    selected_voxels = [];
    
    %tic
    all_voxels = squeeze(input(:,:,:,iTR));
    selected_voxels = all_voxels(index);
    %toc
    
    output(iTR) = mean(selected_voxels);
    
end
%toc

end
            
            