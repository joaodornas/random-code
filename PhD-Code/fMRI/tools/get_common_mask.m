function common_mask = get_common_mask( masks )


nMasks = length(masks);

current_mask = masks(1).mask;

for iMask=2:nMasks
    
    next_mask = masks(iMask).mask;
    
    idx_current_mask = find(current_mask == 1);
    
    idx_cross_mask = find(next_mask(idx_current_mask) == 1);
    
    all_idx_current_mask = find(current_mask < 2);
    
    all_idx_current_mask(idx_current_mask(idx_cross_mask)) = [];
    
    current_mask(all_idx_current_mask) = 0;
    
    common_mask = current_mask;
    
end
  
end

