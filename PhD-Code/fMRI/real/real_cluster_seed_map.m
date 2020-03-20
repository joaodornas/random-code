
%settings_subj1_2210;
%settings_subj2_2610;
%settings_subj3_0311;
%settings_subj4_0211;
%settings_subj5_0211;
%settings_subj6_2411;

subject_label = 'All-Subjects';

%all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','DMN','AUD'};
all_networks_labels = {'LAN','DMN'};

analysis_label = 'Functional-Parcels';
analysis_step_label = 'Pop-Map';

nNet = length(all_networks_labels);

for iNet=1:nNet
   
    network_seeds = getFunctionalSeeds_v5(all_networks_labels{iNet});
    
    nROI = length(network_seeds.ROI);
    
    for iROI=1:nROI
    
        seed_map_label = strcat('LHR','-',subject_label,'-',analysis_label,'-',analysis_step_label,'-',network_seeds.network_label,'-',strrep(network_seeds.ROI(iROI).label,'/','-'));
    
        system(strcat('/usr/local/fsl/bin/cluster --in=',seed_map_label,'.nii --zthresh=3 --pthresh=0.05 --dlh=0.15 --volume=228483 --othresh=',seed_map_label,'-clu',' >',seed_map_label,'-cluster_info.txt'));
        
    end
    
end