
function [ output ] = remmean_whole_brain( input, mask_input )

nTime = size(input,4);

output = zeros(size(input));

for iTime=1:nTime
    
    get_3D_data = squeeze(input(:,:,:,iTime));
    
    get_3D_data_vec = reshape(get_3D_data,[1,size(input,1)*size(input,2)*size(input,3)]);
    
    mask_vec = reshape(mask_input,[1,size(input,1)*size(input,2)*size(input,3)]);
    
    idx_voxels_brain = find(mask_vec == 1);
    
    idx_voxels_space = find(mask_vec == 0);
    
    idx_brain_zero = find(get_3D_data_vec(idx_voxels_brain) == 0);
    idx_voxels_brain(idx_brain_zero) = [];
    
    mean_everything = mean(get_3D_data_vec(idx_voxels_brain));
    
    remmean_3D_data = get_3D_data;
    
    remmean_3D_data(idx_voxels_brain) = get_3D_data(idx_voxels_brain) - mean_everything;
    
    remmean_3D_data(idx_voxels_space) = 0;

    output(:,:,:,iTime) = remmean_3D_data;
    
    clear get_3D_data
    clear get_3D_data_vec
    clear mean_everything
    clear std_everything
    clear remmean_3D_data
    
end

end

