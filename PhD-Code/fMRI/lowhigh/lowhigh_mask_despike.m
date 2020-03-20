

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

threshold = 10;
wavelet = 'd4';
RAM = 100;

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
    
    condition = zeros(nVoxels,nTR);
    
    for iVoxel=1:nVoxels
        
        [idxx,idxy,idxz] = ind2sub(size(maskdata),idx_brain(iVoxel));
        
        condition(iVoxel,1:nTR) = data(idxx,idxy,idxz,1:nTR);
        
    end
    
    %threshold = [5 10 15 20 25];
    threshold = 10;
    wavelet = 'd4';
    
%     for ithreshold=1:length(threshold)
%         
%         disp(strcat('ithreshold:',int2str(ithreshold)));
%         
%         [clean(ithreshold).ts noise(ithreshold).ts] = WaveletDespike(condition,'wavelet',wavelet,'threshold',threshold(ithreshold));
% 
%     end
    
%     save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-',condition_label{ifile},'-Wavelet','.mat'),'clean','noise','condition','idx_brain','maskdata','label','-v7.3');
    
    [clean, noise] = WaveletDespike(condition,'wavelet',wavelet,'threshold',threshold,'LimitRAM',RAM);

    clean_four_D = zeros(size(data));
    noise_four_D = zeros(size(data));
    
    for iVoxel=1:nVoxels
        
        [idxx,idxy,idxz] = ind2sub(size(maskdata),idx_brain(iVoxel));
        
        clean_four_D(idxx,idxy,idxz,1:nTR) = clean(iVoxel,1:nTR);
        
        noise_four_D(idxx,idxy,idxz,1:nTR) = noise(iVoxel,1:nTR);
        
    end
    
    dim = size(data);
    dtype = 'FLOAT32';
    offset = 0;
    scl_slope = 1;
    scl_inter = 0;
    nifti_file = datafile;
    descrip = 'Wavelet';
    
    fname = strcat(filename{ifile},'-clean.nii');
    input_data = clean_four_D;
    lowhigh_save_image;
    
    fname = strcat(filename{ifile},'-noise.nii');
    input_data = noise_four_D;
    lowhigh_save_image;
    
    clear condition
    clear clean
    clear noise
    
end

    