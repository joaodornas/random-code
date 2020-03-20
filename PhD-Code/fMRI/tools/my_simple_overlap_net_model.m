
%%% SIMPLE METHOD TO REMOVE OVERLAPS ON FUNCTIONAL NETWORKS

folder = 'Z:\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually';

nNets = 8;

all_networks_labels = {'DAN','VAN','VIS','LAN','FPC','SMN','AUD','DMN'};

file{1} = 'LHR-All-Subjects-Functional-Parcels-Pop-Map-DAN-zstat.nii';
file{2} = 'LHR-All-Subjects-Functional-Parcels-Pop-Map-VAN-zstat.nii';
file{3} = 'LHR-All-Subjects-Functional-Parcels-Pop-Map-VIS-zstat.nii';
file{4} = 'LHR-All-Subjects-Functional-Parcels-Pop-Map-LAN-zstat.nii';
file{5} = 'LHR-All-Subjects-Functional-Parcels-Pop-Map-FPC-zstat.nii';
file{6} = 'LHR-All-Subjects-Functional-Parcels-Pop-Map-SMN-zstat.nii';
file{7} = 'LHR-All-Subjects-Functional-Parcels-Pop-Map-AUD-zstat.nii';
file{8} = 'LHR-All-Subjects-Functional-Parcels-Pop-Map-DMN-zstat.nii';

MNI_size = [91 109 91];

allNets_img = zeros(nNets,91*109*91);

for iFile=1:nNets
    
    this = nifti(file{iFile});
    
    this.dat.fname = strcat(folder,'\',file{iFile});
    
    this_img = this.dat(:,:,:);
    
    allNets_img(iFile,:) = this_img(:);
    
end

s_allNets_img = squeeze(sum(allNets_img,1));

s_allNets_img = reshape(s_allNets_img,MNI_size);
idx_zeros = find(s_allNets_img==0);

[Y,I] = max(allNets_img,[],1);

I_img = reshape(I,MNI_size);

I_img(idx_zeros) = 0;

nifti_file = this;
offset = this.dat.offset;
scl_slope = this.dat.scl_slope;
scl_inter = this.dat.scl_inter;
dtype = 'FLOAT32';
offset = 0;
dim = this.dat.dim;

descrip = 'Functional-Network';
fname = strcat('Functional-Network-no-overlap-expected','.nii');
input_data = I_img; 
real_save_image;

