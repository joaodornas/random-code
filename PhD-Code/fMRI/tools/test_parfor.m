
for iVoxel=1:nVoxels
    
    voxel(iVoxel).ts = zeros(1,nTR);

    clean(iVoxel).ts = zeros(1,nTR);
    noise(iVoxel).ts = zeros(1,nTR);

end

idxx = zeros(1,nTR);
idxy = zeros(1,nTR);
idxz = zeros(1,nTR);

parfor iVoxel=1:nVoxels

    iVoxel
    
    [idxx(iVoxel),idxy(iVoxel),idxz(iVoxel)] = ind2sub(size(maskdata),idx_brain(iVoxel));

    voxel(iVoxel).ts = squeeze(data(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),1:nTR));

    [clean(iVoxel).ts, noise(iVoxel).ts] = WaveletDespike(voxel(iVoxel).ts,'wavelet',wavelet,'threshold',threshold,'verbose',0);

end
