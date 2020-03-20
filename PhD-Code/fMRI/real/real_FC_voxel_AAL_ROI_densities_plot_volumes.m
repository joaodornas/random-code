
function real_FC_voxel_AAL_ROI_densities_plot_volumes

%plotOnlyResultsVolumes;
plotResultsWithNets;

end

function plotOnlyResultsVolumes

    results_folder = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\densities\rendering\both-flipped\clusters-dlh=1.5-both-fl-2mm-FS_spm_CanonicalBrain_norecon-unzip';


    filename{1} = 'LHR-Correlation-Contrast-Passive-Neg-z-clu-both-fl-2-mm';
    filename{2} = 'LHR-Correlation-Contrast-Passive-Pos-z-clu-both-fl-2-mm';
    filename{3} = 'LHR-Correlation-Contrast-Track-Neg-z-clu-both-fl-2-mm';
    filename{4} = 'LHR-Correlation-Contrast-Track-Pos-z-clu-both-fl-2-mm';

    for ifile=1:length(filename)

        datafile = nifti(strcat(results_folder,'\',filename{ifile},'.nii'));
        datafile.dat.fname = strcat(results_folder,'\',filename{ifile},'.nii');

        img = datafile.dat(:,:,:);

        idx_negative = find(img < 0);
        idx_positive = find(img > 0);

        new_img = zeros(size(img));

        new_img(idx_negative) = 1;
        new_img(idx_positive) = 2;

        nifti_file = datafile;
        offset = datafile.dat.offset;
        scl_slope = datafile.dat.scl_slope;
        scl_inter = datafile.dat.scl_inter;

        dtype = 'FLOAT32';
        offset = 0;

        dim = datafile.dat.dim;

        descrip = 'paint';

        fname = strcat(results_folder,'\',filename{ifile},'-paint','.nii');
        input_data = new_img; 
        real_save_image;

    end

end

function plotResultsWithNets

    results_folder = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION\all-subjects\Real\FC_Voxels_AAL_ROI\densities\rendering\both-flipped\clusters-dlh=1.5-both-fl-2mm-FS_spm_CanonicalBrain_norecon-unzip';

    nets_folder = 'Z:\Dropbox (Uni Magdeburg)\_DATA\Parcellation\Functional_Parcellation\v4\Final_Parcellation\net_by_net_individually';
    
    filename{1} = 'LHR-Correlation-Contrast-Passive-Neg-z-clu-both-fl-2-mm';
    filename{2} = 'LHR-Correlation-Contrast-Passive-Pos-z-clu-both-fl-2-mm';
    filename{3} = 'LHR-Correlation-Contrast-Track-Neg-z-clu-both-fl-2-mm';
    filename{4} = 'LHR-Correlation-Contrast-Track-Pos-z-clu-both-fl-2-mm';
    
    DAN_file = strcat(nets_folder,'\','DAN-bin','.nii');
    VAN_file = strcat(nets_folder,'\','VAN-bin','.nii');
    FPC_file = strcat(nets_folder,'\','FPC-bin','.nii');

    for ifile=1:length(filename)

        datafile = nifti(strcat(results_folder,'\',filename{ifile},'.nii'));
        datafile.dat.fname = strcat(results_folder,'\',filename{ifile},'.nii');

        DAN_datafile = nifti(DAN_file);
        VAN_datafile = nifti(VAN_file);
        FPC_datafile = nifti(FPC_file);
        
        DAN_datafile.dat.fname = DAN_file;
        VAN_datafile.dat.fname = VAN_file;
        FPC_datafile.dat.fname = FPC_file;
        
        img = datafile.dat(:,:,:);

        DAN_img = DAN_datafile.dat(:,:,:);
        VAN_img = VAN_datafile.dat(:,:,:);
        FPC_img = FPC_datafile.dat(:,:,:);
        
        idx_DAN = find(DAN_img == 1);
        idx_VAN = find(VAN_img == 1);
        idx_FPC = find(FPC_img == 1);
        
        idx_negative = find(img < 0);
        idx_positive = find(img > 0);

        %%% DAN
        
        new_img_DAN = zeros(size(img));

        new_img_DAN(idx_DAN) = 3;
        new_img_DAN(idx_negative) = 1;
        new_img_DAN(idx_positive) = 2;

        nifti_file = datafile;
        offset = datafile.dat.offset;
        scl_slope = datafile.dat.scl_slope;
        scl_inter = datafile.dat.scl_inter;

        dtype = 'FLOAT32';
        offset = 0;

        dim = datafile.dat.dim;

        descrip = 'paint';

        fname = strcat(results_folder,'\',filename{ifile},'-paint-DAN','.nii');
        input_data = new_img_DAN; 
        real_save_image;

        %%% VAN
        
        new_img_VAN = zeros(size(img));

        new_img_VAN(idx_VAN) = 3;
        new_img_VAN(idx_negative) = 1;
        new_img_VAN(idx_positive) = 2;

        nifti_file = datafile;
        offset = datafile.dat.offset;
        scl_slope = datafile.dat.scl_slope;
        scl_inter = datafile.dat.scl_inter;

        dtype = 'FLOAT32';
        offset = 0;

        dim = datafile.dat.dim;

        descrip = 'paint';

        fname = strcat(results_folder,'\',filename{ifile},'-paint-VAN','.nii');
        input_data = new_img_VAN; 
        real_save_image;
        
        %%% FPC
        
        new_img_FPC = zeros(size(img));

        new_img_FPC(idx_FPC) = 3;
        new_img_FPC(idx_negative) = 1;
        new_img_FPC(idx_positive) = 2;

        nifti_file = datafile;
        offset = datafile.dat.offset;
        scl_slope = datafile.dat.scl_slope;
        scl_inter = datafile.dat.scl_inter;

        dtype = 'FLOAT32';
        offset = 0;

        dim = datafile.dat.dim;

        descrip = 'paint';

        fname = strcat(results_folder,'\',filename{ifile},'-paint-FPC','.nii');
        input_data = new_img_FPC; 
        real_save_image;

    end

end

