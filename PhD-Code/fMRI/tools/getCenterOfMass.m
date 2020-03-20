function [center_of_X, center_of_Y, center_of_Z] = getCenterOfMass(idx,AAL_img)

nIndices = length(idx);

for iidx=1:nIndices
    
    [X(iidx),Y(iidx),Z(iidx)] = ind2sub(size(AAL_img),idx(iidx));
    
end

center_of_X = round(mean(X));
center_of_Y = round(mean(Y));
center_of_Z = round(mean(Z));

end