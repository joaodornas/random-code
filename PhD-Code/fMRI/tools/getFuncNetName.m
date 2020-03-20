function [nets_names, seeds_names] = getFuncNetName(X,Y,Z)

%%% LOAD AAL

%load_funcnets = nifti('All-8-Nets-merged.nii');
load_funcnets = nifti('All-8-Nets-seeds-merged.nii');
load_funcnets.dat.fname = strcat('/Volumes/dropbox/_DATA/Parcellation/Functional_Parcellation/v4/Final_Parcellation/net_by_net_individually/',load_funcnets.dat.fname);

Func_Net_img = load_funcnets.dat(:,:,:,:);

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

MNI_x_center = 45;
MNI_y_center = 63;
MNI_z_center = 36;

size_voxels_mm = 2;

x_center = MNI_x_center - round(X/size_voxels_mm);
y_center = MNI_y_center + round(Y/size_voxels_mm);
z_center = MNI_z_center + round(Z/size_voxels_mm);

idx_NET = find(squeeze(Func_Net_img(x_center + 1,y_center + 1,z_center + 1,:)));

if ~isempty(idx_NET)

    nets_names = all_networks_labels{idx_NET(1)};
    
    idx_seed = Func_Net_img(x_center + 1,y_center + 1,z_center + 1,idx_NET(1)) - idx_NET(1)*100;
    
    seeds = getFunctionalSeeds_v4(all_networks_labels{idx_NET(1)});
    
    seeds_names = seeds.ROI(idx_seed).label;
    
    for i=2:length(idx_NET)
        
        new_name = all_networks_labels{idx_NET(i)};

        nets_names = strcat(nets_names,'-',new_name);
        
        idx_seed = Func_Net_img(x_center + 1,y_center + 1,z_center + 1,idx_NET(i)) - idx_NET(i)*100;
    
        seeds = getFunctionalSeeds_v4(all_networks_labels{idx_NET(i)});
    
        new_seeds_names = seeds.ROI(idx_seed).label;
        
        seeds_names = strcat(seeds_names,'-',new_seeds_names);
    
    end

else
    
    nets_names = 'no nets';
    seeds_names = 'no seeds'; 

end


end

