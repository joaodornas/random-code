

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

net_folder = 'Z:\Dropbox (Uni Magdeburg)\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually';

for iNet=1:length(all_networks_labels)
   
    filename = strcat(all_networks_labels{iNet},'-','bin','.nii');
    
    img_file = nifti(filename);
    
    img_file.dat.fname = strcat(net_folder,'\',filename);
    
    img = img_file.dat(:,:,:);
    
    nVoxels = length(find(img == 1));
    
    all_count{iNet,1} = all_networks_labels{iNet};
    all_count{iNet,2} = nVoxels;
    
end

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

nVoxels = length(find(AAL_img));

all_count{iNet+1,1} = 'AAL';
all_count{iNet+1,2} = nVoxels;


