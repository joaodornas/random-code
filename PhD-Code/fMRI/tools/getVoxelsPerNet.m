
function all_nets = getVoxelsPerNet

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','DMN','AUD'};

nNet = 8;

for iNet=1:nNet
    
   filename = strcat(all_networks_labels{iNet},'-bin.nii');
   
   foldername = 'Z:\Dropbox (Uni Magdeburg)\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually';
    
   load_parcellation = nifti(filename);
   load_parcellation.dat.fname = strcat(foldername,'\',load_parcellation.dat.fname);

   all_nets(iNet).img = load_parcellation.dat(:,:,:);
   
   all_nets(iNet).idx_voxels = find(Net(iNet).img);
   
end


end

