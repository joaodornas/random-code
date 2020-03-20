function dat = fromImageToSPMdat(imgfile)
    
file = nifti(imgfile);
img = file.dat(:,:,:);
        
img_vec = reshape(img,[1,size(img,1)*size(img,2)*size(img,3)]);

[idxx, idxy, idxz] = ind2sub(size(img),1:size(img,1)*size(img,2)*size(img,3));

xSPM.XYZ = zeros(3,size(img,1)*size(img,2)*size(img,3));

xSPM.XYZ(1,:) = idxx;
xSPM.XYZ(2,:) = idxy;
xSPM.XYZ(3,:) = idxz;

xSPM.Z = zeros(1,size(img,1)*size(img,2)*size(img,3));

xSPM.Z(:) = img_vec(:);

xSPM.M = file.mat;

xSPM.DIM = file.dat.dim.';
   
dat = struct( 'XYZ', xSPM.XYZ,...
    't', xSPM.Z',...
    'mat', xSPM.M,...
    'dim',  xSPM.DIM);
    
end