
function [ output_int, output_double ] = zscore_whole_brain( input, mask_input, scale )

%scale = 9;

nTime = size(input,4);

output_int = int32(zeros(size(input)));
output_double = double(zeros(size(input)));

for iTime=1:nTime
    
    %if mod(iTime,2) == 0; disp(int2str(iTime)); end
    
    get_3D_data = squeeze(input(:,:,:,iTime));
    
    get_3D_data_vec = reshape(get_3D_data,[1,size(input,1)*size(input,2)*size(input,3)]);
    
    mask_vec = reshape(mask_input,[1,size(input,1)*size(input,2)*size(input,3)]);
    
    idx_voxels_brain = find(mask_vec == 1);
    
    idx_voxels_space = find(mask_vec == 0);
    
    voxels_data = get_3D_data_vec(idx_voxels_brain);
    idx_non_zero = find(voxels_data ~= 0);
    voxels_data_zscored = zeros(1,length(voxels_data));
    voxels_data_zscored(idx_non_zero) = zscore(voxels_data(idx_non_zero)); 
    
    new_get_3D_data_vec_int = int32(zeros(1,size(input,1)*size(input,2)*size(input,3)));
    new_get_3D_data_vec_int(idx_voxels_brain) = int32(voxels_data_zscored*10^(scale));
    
    new_get_3D_data_vec_double = zeros(1,size(input,1)*size(input,2)*size(input,3));
    new_get_3D_data_vec_double(idx_voxels_brain) = voxels_data_zscored(:);
    
    zscored_3D_data_int = reshape(new_get_3D_data_vec_int,[size(input,1),size(input,2),size(input,3)]);
    zscored_3D_data_double = reshape(new_get_3D_data_vec_double,[size(input,1),size(input,2),size(input,3)]);
    
    output_int(:,:,:,iTime) = zscored_3D_data_int;
    output_double(:,:,:,iTime) = zscored_3D_data_double;
    
    clear get_3D_data
    clear get_3D_data_vec
    clear mean_everything
    clear std_everything
    clear zscored_3D_data
    
end

end

