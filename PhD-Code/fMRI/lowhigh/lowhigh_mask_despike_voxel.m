

%settings_jan_0805;
settings_elena_2905;

%get_at_this_preprocessed_step = settings.FSL.folders.warped;
get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.filtered;
mask = settings.FSL.files.mask.custom;

label = file;

settings = lowhigh_concatenate_folders_strings_FSL(settings,get_at_this_preprocessed_step);

run = 1;
filename{1} = strcat(settings.mot4.FSL.run(run).fullpath,'\',file);
filename{2} = strcat(settings.mot2.FSL.run(run).fullpath,'\',file);
filename{3} = strcat(settings.restingstate.FSL.run(run).fullpath,'\',file);
maskname{1} = strcat(settings.mot4.FSL.run(run).fullpath,'\',file);
maskname{2} = strcat(settings.mot2.FSL.run(run).fullpath,'\',file);
maskname{3} = strcat(settings.restingstate.FSL.run(run).fullpath,'\',file);

run = 2;
filename{4} = strcat(settings.mot4.FSL.run(run).fullpath,'\',file);
filename{5} = strcat(settings.mot2.FSL.run(run).fullpath,'\',file);
filename{6} = strcat(settings.restingstate.FSL.run(run).fullpath,'\',file);
maskname{4} = strcat(settings.mot4.FSL.run(run).fullpath,'\',file);
maskname{5} = strcat(settings.mot2.FSL.run(run).fullpath,'\',file);
maskname{6} = strcat(settings.restingstate.FSL.run(run).fullpath,'\',file);

condition_label{1} = 'MOT4-Run1';
condition_label{2} = 'MOT2-Run1';
condition_label{3} = 'RestingState-Run1';
condition_label{4} = 'MOT4-Run2';
condition_label{5} = 'MOT2-Run2';
condition_label{6} = 'RestingState-Run2';

for ifile=1:length(filename)
%for ifile=4
    
    disp(strcat('ifile:',int2str(ifile)));
    
    datafile = nifti(strcat(filename{ifile},'.nii'));
    datafile.dat.fname = strcat(filename{ifile},'.nii');
    data = datafile.dat(:,:,:,:);
    
    maskfile = nifti(strcat(maskname{ifile},'.nii'));
    maskfile.dat.fname = strcat(maskname{ifile},'.nii');
    maskdata = maskfile.dat(:,:,:);
    
    nTR = size(data,4);
    idx_brain = find(maskdata);
    nVoxels = length(idx_brain);

    threshold = 100;
    wavelet = 'd4';
    
    clean_four_D = zeros(size(data));
    noise_four_D = zeros(size(data));
    
    voxel(1:nVoxels) = struct('ts',zeros(1,nTR));
    
    clean(1:nVoxels) = struct('ts',zeros(1,nTR));
    noise(1:nVoxels) = struct('ts',zeros(1,nTR));
    
    idxx = zeros(1,nTR);
    idxy = zeros(1,nTR);
    idxz = zeros(1,nTR);
    
    k = 0;
    
    parfor_progress(nVoxels);
    
    parfor iVoxel=1:nVoxels
        
        [idxx(iVoxel),idxy(iVoxel),idxz(iVoxel)] = ind2sub(size(maskdata),idx_brain(iVoxel));
        
        voxel(iVoxel).ts = squeeze(data(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),1:nTR));
        
        [clean(iVoxel).ts, noise(iVoxel).ts] = WaveletDespike(voxel(iVoxel).ts,'wavelet',wavelet,'threshold',threshold,'verbose',0);
        
        parfor_progress;
        
    end
    
    for iVoxel=1:nVoxels
        
        clean_four_D(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),1:nTR) = clean(iVoxel).ts(1:nTR);
        
        noise_four_D(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),1:nTR) = noise(iVoxel).ts(1:nTR);
        
    end
    
    dim = size(data);
    dtype = 'FLOAT32';
    offset = 0;
    scl_slope = 1;
    scl_inter = 0;
    nifti_file = datafile;
    descrip = 'Wavelet';
    
    fname = strcat(filename{ifile},'-clean-voxel.nii');
    input_data = clean_four_D;
    lowhigh_save_image;
    
    fname = strcat(filename{ifile},'-noise-voxel.nii');
    input_data = noise_four_D;
    lowhigh_save_image;
    
end

    