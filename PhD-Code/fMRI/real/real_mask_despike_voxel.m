function real_mask_despike_voxel

% all_settings = getAllSettingsPetra;
% 
% nSettings = length(all_settings);
% 
% for iSet=1:nSettings
% 
%     settings = all_settings(iSet).settings;
% 
%     doTheMath(settings);
% 
% end

% OnlyOneSubjectAndRun;

% settings_martin;
% doTheMathMARTIN(settings);

% doTheMathMARTIN2;

% settings_subj5_again;
% doTheMathSUBJ5Again(settings);

% doTheMathMARTINCustom;

doTheMathSUBJ5AgainCustomPrefiltered;

end

function doTheMath(settings)

get_at_this_preprocessed_step = settings.FSL.folders.melodic;
file = settings.FSL.files.functional.warped;
mask = settings.FSL.files.mask.warped;

label = file;

settings = real_concatenate_folders_strings_FSL_Petra(settings,get_at_this_preprocessed_step);

run = 1;

filename{1} = strcat(settings.restingstate.FSL.run(run).fullpath,'\',file);

maskname{1} = strcat(settings.restingstate.FSL.run(run).fullpath,'\',mask);

condition_label{1} = strcat('RestingState-Run',int2str(run));

for ifile=1:length(filename)
    
    disp(strcat('ifile:',int2str(ifile)));
    start = now;
    disp(datestr(start));
    
    datafile = nifti(strcat(filename{ifile},'.nii'));
    datafile.dat.fname = strcat(filename{ifile},'.nii');
    data = datafile.dat(:,:,:,:);
    
    maskfile = nifti(strcat(maskname{ifile},'.nii'));
    maskfile.dat.fname = strcat(maskname{ifile},'.nii');
    maskdata = maskfile.dat(:,:,:);
    
    nTR = size(data,4);
    idx_brain = find(maskdata);
    nVoxels = length(idx_brain);
    
    disp(strcat('nVoxels:',int2str(nVoxels)));

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
    
    %parfor_progress(nVoxels);
    
    parfor iVoxel=1:nVoxels
        
        %if mod(iVoxel,1000) == 0; disp(strcat('iVoxel:',int2str(iVoxel))); end
        
        [idxx(iVoxel),idxy(iVoxel),idxz(iVoxel)] = ind2sub(size(maskdata),idx_brain(iVoxel));
        
        voxel(iVoxel).ts = squeeze(data(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),1:nTR));
        
        [clean(iVoxel).ts, noise(iVoxel).ts] = WaveletDespike(voxel(iVoxel).ts,'wavelet',wavelet,'threshold',threshold,'verbose',0);
        
        %parfor_progress;
        
    end
    
    %parfor_progress(0);
    
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
    real_save_image;
    
    fname = strcat(filename{ifile},'-noise-voxel.nii');
    input_data = noise_four_D;
    real_save_image;
    
    end_ = now;
    disp(strcat('Lasted:',datestr(end_-start)));
    
end

end

function OnlyOneSubjectAndRun

fullpath = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-1-22-10-2015\preprocessed\T2-Stimulus-RestingState\Run-1-4-1\FSL\Melodic.ica';

file = 'filtered_func_data2standard.nii';
mask = 'mask2standard.nii';

filename = strcat(fullpath,'\',file);

maskname = strcat(fullpath,'\',mask);

start = now;
disp(datestr(start));
    
datafile = nifti(filename);
datafile.dat.fname = filename;
data = datafile.dat(:,:,:,:);
    
maskfile = nifti(maskname);
maskfile.dat.fname = maskname;
maskdata = maskfile.dat(:,:,:);
    
nTR = size(data,4);
idx_brain = find(maskdata);
nVoxels = length(idx_brain);
    
disp(strcat('nVoxels:',int2str(nVoxels)));

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
    
parfor iVoxel=1:nVoxels

    [idxx(iVoxel),idxy(iVoxel),idxz(iVoxel)] = ind2sub(size(maskdata),idx_brain(iVoxel));

    voxel(iVoxel).ts = squeeze(data(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),1:nTR));

    [clean(iVoxel).ts, noise(iVoxel).ts] = WaveletDespike(voxel(iVoxel).ts,'wavelet',wavelet,'threshold',threshold,'verbose',0);

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

fname = strcat(filename,'-clean-voxel.nii');
input_data = clean_four_D;
real_save_image;

fname = strcat(filename,'-noise-voxel.nii');
input_data = noise_four_D;
real_save_image;

end_ = now;
disp(strcat('Lasted:',datestr(end_-start)));
    
end

function doTheMathMARTIN(settings)

get_at_this_preprocessed_step = settings.FSL.folders.melodic;
file = settings.FSL.files.functional.warped;
mask = settings.FSL.files.mask.warped;

label = file;

settings = real_concatenate_folders_strings_FSL_Martin(settings,get_at_this_preprocessed_step);

nRuns = 4;

for iRun=1:nRuns;

    filename{iRun} = strcat(settings.restingstate.FSL.run(iRun).fullpath,'\',file);

    maskname{iRun} = strcat(settings.restingstate.FSL.run(iRun).fullpath,'\',mask);

    condition_label{iRun} = strcat('RestingState-Run',int2str(iRun));

    for ifile=1:length(filename)

        disp(strcat('ifile:',int2str(ifile)));
        start = now;
        disp(datestr(start));

        datafile = nifti(strcat(filename{ifile},'.nii'));
        datafile.dat.fname = strcat(filename{ifile},'.nii');
        data = datafile.dat(:,:,:,:);

        maskfile = nifti(strcat(maskname{ifile},'.nii'));
        maskfile.dat.fname = strcat(maskname{ifile},'.nii');
        maskdata = maskfile.dat(:,:,:);

        nTR = size(data,4);
        idx_brain = find(maskdata);
        nVoxels = length(idx_brain);

        disp(strcat('nVoxels:',int2str(nVoxels)));

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

        %parfor_progress(nVoxels);

        parfor iVoxel=1:nVoxels

            %if mod(iVoxel,1000) == 0; disp(strcat('iVoxel:',int2str(iVoxel))); end

            [idxx(iVoxel),idxy(iVoxel),idxz(iVoxel)] = ind2sub(size(maskdata),idx_brain(iVoxel));

            voxel(iVoxel).ts = squeeze(data(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),1:nTR));

            [clean(iVoxel).ts, noise(iVoxel).ts] = WaveletDespike(voxel(iVoxel).ts,'wavelet',wavelet,'threshold',threshold,'verbose',0);

            %parfor_progress;

        end

        %parfor_progress(0);

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
        real_save_image;

        fname = strcat(filename{ifile},'-noise-voxel.nii');
        input_data = noise_four_D;
        real_save_image;

        end_ = now;
        disp(strcat('Lasted:',datestr(end_-start)));

    end

end

end

function doTheMathSUBJ5Again(settings)

get_at_this_preprocessed_step = settings.FSL.folders.melodic;
file = settings.FSL.files.functional.warped;
mask = settings.FSL.files.mask.warped;

label = file;

settings = real_concatenate_folders_strings_FSL_SUBJ5_Again(settings,get_at_this_preprocessed_step);

nRuns = 4;

for iRun=1:nRuns;

    filename{iRun} = strcat(settings.restingstate.FSL.run(iRun).fullpath,'\',file);

    maskname{iRun} = strcat(settings.restingstate.FSL.run(iRun).fullpath,'\',mask);

    condition_label{iRun} = strcat('RestingState-Run',int2str(iRun));
    
end

    for ifile=3:length(filename)

        disp(strcat('ifile:',int2str(ifile)));
        start = now;
        disp(datestr(start));

        datafile = nifti(strcat(filename{ifile},'.nii'));
        datafile.dat.fname = strcat(filename{ifile},'.nii');
        data = datafile.dat(:,:,:,:);

        maskfile = nifti(strcat(maskname{ifile},'.nii'));
        maskfile.dat.fname = strcat(maskname{ifile},'.nii');
        maskdata = maskfile.dat(:,:,:);

        nTR = size(data,4);
        idx_brain = find(maskdata);
        nVoxels = length(idx_brain);

        disp(strcat('nVoxels:',int2str(nVoxels)));

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

        %parfor_progress(nVoxels);

        parfor iVoxel=1:nVoxels

            %if mod(iVoxel,1000) == 0; disp(strcat('iVoxel:',int2str(iVoxel))); end

            [idxx(iVoxel),idxy(iVoxel),idxz(iVoxel)] = ind2sub(size(maskdata),idx_brain(iVoxel));

            voxel(iVoxel).ts = squeeze(data(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),1:nTR));

            [clean(iVoxel).ts, noise(iVoxel).ts] = WaveletDespike(voxel(iVoxel).ts,'wavelet',wavelet,'threshold',threshold,'verbose',0);

            %parfor_progress;

        end

        %parfor_progress(0);

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
        real_save_image;

        fname = strcat(filename{ifile},'-noise-voxel.nii');
        input_data = noise_four_D;
        real_save_image;

        end_ = now;
        disp(strcat('Lasted:',datestr(end_-start)));

    end

end

function doTheMathMARTIN2

nRuns = 4;

file = 'filtered_func_data2standard';
mask = 'mask2standard';

for iRun=1:nRuns;

    filename{iRun} = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\MARTIN\preprocessed\T2-Stimulus-RestingState\dn20_1622\Run-',int2str(iRun),'\FSL\Melodic-Fieldmap.ica','\',file);

    maskname{iRun} = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\MARTIN\preprocessed\T2-Stimulus-RestingState\dn20_1622\Run-',int2str(iRun),'\FSL\Melodic-Fieldmap.ica','\',mask);

    condition_label{iRun} = strcat('RestingState-Run',int2str(iRun));
    
end

    for ifile=1:length(filename)

        disp(strcat('ifile:',int2str(ifile)));
        start = now;
        disp(datestr(start));

        datafile = nifti(strcat(filename{ifile},'.nii'));
        datafile.dat.fname = strcat(filename{ifile},'.nii');
        data = datafile.dat(:,:,:,:);

        maskfile = nifti(strcat(maskname{ifile},'.nii'));
        maskfile.dat.fname = strcat(maskname{ifile},'.nii');
        maskdata = maskfile.dat(:,:,:);

        nTR = size(data,4);
        idx_brain = find(maskdata);
        nVoxels = length(idx_brain);

        disp(strcat('nVoxels:',int2str(nVoxels)));

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

        %parfor_progress(nVoxels);

        parfor iVoxel=1:nVoxels

            %if mod(iVoxel,1000) == 0; disp(strcat('iVoxel:',int2str(iVoxel))); end

            [idxx(iVoxel),idxy(iVoxel),idxz(iVoxel)] = ind2sub(size(maskdata),idx_brain(iVoxel));

            voxel(iVoxel).ts = squeeze(data(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),1:nTR));

            [clean(iVoxel).ts, noise(iVoxel).ts] = WaveletDespike(voxel(iVoxel).ts,'wavelet',wavelet,'threshold',threshold,'verbose',0);

            %parfor_progress;

        end

        %parfor_progress(0);

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
        real_save_image;

        fname = strcat(filename{ifile},'-noise-voxel.nii');
        input_data = noise_four_D;
        real_save_image;

        end_ = now;
        disp(strcat('Lasted:',datestr(end_-start)));

    end

end

function doTheMathMARTINCustom

nRuns = 4;

file = 'filtered_func_data_mcf_unwarp2standard';
mask = 'mask2standard';

for iRun=1:nRuns;

    filename{iRun} = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\MARTIN\preprocessed\T2-Stimulus-RestingState\dn20_1590\Run-',int2str(iRun),'\FSL\custom','\',file);

    maskname{iRun} = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\MARTIN\preprocessed\T2-Stimulus-RestingState\dn20_1590\Run-',int2str(iRun),'\FSL\custom','\',mask);

    condition_label{iRun} = strcat('RestingState-Run',int2str(iRun));
    
end

doTheWavelet(filename,maskname);

clear filename maskname

for iRun=1:nRuns;

    filename{iRun} = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\MARTIN\preprocessed\T2-Stimulus-RestingState\dn20_1622\Run-',int2str(iRun),'\FSL\custom','\',file);

    maskname{iRun} = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\MARTIN\preprocessed\T2-Stimulus-RestingState\dn20_1622\Run-',int2str(iRun),'\FSL\custom','\',mask);

    condition_label{iRun} = strcat('RestingState-Run',int2str(iRun));
    
end

doTheWavelet(filename,maskname);

end

function doTheMathSUBJ5AgainCustomPrefiltered

nRuns = 4;

file = 'filtered_func_data_mcf_unwarp2standard';
mask = 'maskprefiltered2standard';

filename{1} = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-5-preprocess-test\preprocessed\T2-Stimulus-RestingState\Run-1-4-1\FSL\custom','\',file);
maskname{1} = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-5-preprocess-test\preprocessed\T2-Stimulus-RestingState\Run-1-4-1\FSL\custom','\',mask);

filename{2} = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-5-preprocess-test\preprocessed\T2-Stimulus-RestingState\Run-2-8-1\FSL\custom','\',file);
maskname{2} = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-5-preprocess-test\preprocessed\T2-Stimulus-RestingState\Run-2-8-1\FSL\custom','\',mask);

filename{3} = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-5-preprocess-test\preprocessed\T2-Stimulus-RestingState\Run-3-4-2\FSL\custom','\',file);
maskname{3} = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-5-preprocess-test\preprocessed\T2-Stimulus-RestingState\Run-3-4-2\FSL\custom','\',mask);

filename{4} = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-5-preprocess-test\preprocessed\T2-Stimulus-RestingState\Run-4-8-2\FSL\custom','\',file);
maskname{4} = strcat('Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\SUBJECT-5-preprocess-test\preprocessed\T2-Stimulus-RestingState\Run-4-8-2\FSL\custom','\',mask);

doTheWaveletPre(filename,maskname);

end

function doTheWaveletPre(filename,maskname)

    for ifile=1:length(filename)

        disp(strcat('ifile:',int2str(ifile)));
        start = now;
        disp(datestr(start));

        datafile = nifti(strcat(filename{ifile},'.nii'));
        datafile.dat.fname = strcat(filename{ifile},'.nii');
        data = datafile.dat(:,:,:,:);

        maskfile = nifti(strcat(maskname{ifile},'.nii'));
        maskfile.dat.fname = strcat(maskname{ifile},'.nii');
        maskdata = maskfile.dat(:,:,:);

        nTR = size(data,4);
        idx_brain = find(maskdata);
        nVoxels = length(idx_brain);

        disp(strcat('nVoxels:',int2str(nVoxels)));

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

        %parfor_progress(nVoxels);

        parfor iVoxel=1:nVoxels

            %if mod(iVoxel,1000) == 0; disp(strcat('iVoxel:',int2str(iVoxel))); end

            [idxx(iVoxel),idxy(iVoxel),idxz(iVoxel)] = ind2sub(size(maskdata),idx_brain(iVoxel));

            voxel(iVoxel).ts = squeeze(data(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),1:nTR));

            [clean(iVoxel).ts, noise(iVoxel).ts] = WaveletDespike(voxel(iVoxel).ts,'wavelet',wavelet,'threshold',threshold,'verbose',0);

            %parfor_progress;

        end

        %parfor_progress(0);

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

        fname = strcat(filename{ifile},'-clean-voxel-pre.nii');
        input_data = clean_four_D;
        real_save_image;

        fname = strcat(filename{ifile},'-noise-voxel-pre.nii');
        input_data = noise_four_D;
        real_save_image;

        end_ = now;
        disp(strcat('Lasted:',datestr(end_-start)));

    end

end

function doTheWavelet(filename,maskname)

    for ifile=1:length(filename)

        disp(strcat('ifile:',int2str(ifile)));
        start = now;
        disp(datestr(start));

        datafile = nifti(strcat(filename{ifile},'.nii'));
        datafile.dat.fname = strcat(filename{ifile},'.nii');
        data = datafile.dat(:,:,:,:);

        maskfile = nifti(strcat(maskname{ifile},'.nii'));
        maskfile.dat.fname = strcat(maskname{ifile},'.nii');
        maskdata = maskfile.dat(:,:,:);

        nTR = size(data,4);
        idx_brain = find(maskdata);
        nVoxels = length(idx_brain);

        disp(strcat('nVoxels:',int2str(nVoxels)));

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

        %parfor_progress(nVoxels);

        parfor iVoxel=1:nVoxels

            %if mod(iVoxel,1000) == 0; disp(strcat('iVoxel:',int2str(iVoxel))); end

            [idxx(iVoxel),idxy(iVoxel),idxz(iVoxel)] = ind2sub(size(maskdata),idx_brain(iVoxel));

            voxel(iVoxel).ts = squeeze(data(idxx(iVoxel),idxy(iVoxel),idxz(iVoxel),1:nTR));

            [clean(iVoxel).ts, noise(iVoxel).ts] = WaveletDespike(voxel(iVoxel).ts,'wavelet',wavelet,'threshold',threshold,'verbose',0);

            %parfor_progress;

        end

        %parfor_progress(0);

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
        real_save_image;

        fname = strcat(filename{ifile},'-noise-voxel.nii');
        input_data = noise_four_D;
        real_save_image;

        end_ = now;
        disp(strcat('Lasted:',datestr(end_-start)));

    end

end


    