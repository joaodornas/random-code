%%% generate MD758 - Funcitonal Parcellation (Seeds and Population Maps)
%%% Overlaps

%%% GET INFORMATION ABOUT SEEDS OVERLAP WITH CLUSTERS

clear all

load('/Volumes/INDIREA/_DATA/Parcellation/758-Cluster/Corbetta/ROI_MD758Corbetta.mat');
% load('ROI_MD758Corbetta.mat');

networks = {'DAN' 'VAN' 'VIS' 'LAN' 'FPC' 'SMN' 'AUD' 'DMN'};
nNets = 8;
nROIs = 90;

for iNet=1:nNets
    
    clusters_per_net(iNet).net_label = networks{iNet};
    
end

for iNet=1:nNets
    
    net_label = networks{iNet};
    
    %%% GET CLUSTERS WITHIN NETWORK
    
    iiselected = 0;
    iiCluster = 0;
    
    for iROI=1:nROIs
       
        nClusters = ROI_MD758Corbetta(iROI).nClusters;
        
        for iCluster=1:nClusters
            
            iiCluster = iiCluster + 1;
            
            if isfield(ROI_MD758Corbetta(iROI).clusters(iCluster),'Nets')
           
                this_cluster_nets = ROI_MD758Corbetta(iROI).clusters(iCluster).Nets;

                if iscell(this_cluster_nets)

                    for iCell=1:length(this_cluster_nets)

                        this_net = this_cluster_nets{iCell};

                        if strfind(this_net,net_label)

                            iiselected = iiselected + 1;
                            
                            seed_label = this_net(5:end);
                            
                            idx_slash = strfind(seed_label,'-');
                            
                            if length(idx_slash) > 1
                                
                                seed_label = seed_label(idx_slash(1)+1:idx_slash(2)-1);
                                
                            elseif ~isempty(idx_slash)
                                
                                seed_label = seed_label(idx_slash(1)+1:end);
                                
                            end

                            clusters_per_net(iNet).clusters(iiselected) = iiCluster;
                            clusters_per_net(iNet).seed_label{iiselected} = seed_label;
                        
                        end

                    end

            end
            
            end
            
        end
        
    end
    
end


%%% GETS INFORMATION ABOUT POPULATION MAP OVERLAP WITH CLUSTER

% folder_nets = 'Z:\_DATA\Parcellation\Functional_Parcellation\v5\Final_Parcellation\net_by_net_individually\';
folder_nets = '/Volumes/INDIREA/_DATA/Parcellation/Functional_Parcellation/v5/Final_Parcellation/net_by_net_individually/';

DAN = nifti(strcat(folder_nets,'DAN-bin.nii'));
DAN.dat.fname = strcat(folder_nets,'DAN-bin.nii');
DAN_img = DAN.dat(:,:,:);

VAN = nifti(strcat(folder_nets,'VAN-bin.nii'));
VAN.dat.fname = strcat(folder_nets,'VAN-bin.nii');
VAN_img = VAN.dat(:,:,:);

VIS = nifti(strcat(folder_nets,'VIS-bin.nii'));
VIS.dat.fname = strcat(folder_nets,'VIS-bin.nii');
VIS_img = VIS.dat(:,:,:);

AUD = nifti(strcat(folder_nets,'AUD-bin.nii'));
AUD.dat.fname = strcat(folder_nets,'AUD-bin.nii');
AUD_img = AUD.dat(:,:,:);

LAN = nifti(strcat(folder_nets,'LAN-bin.nii'));
LAN.dat.fname = strcat(folder_nets,'LAN-bin.nii');
LAN_img = LAN.dat(:,:,:);

FPC = nifti(strcat(folder_nets,'FPC-bin.nii'));
FPC.dat.fname = strcat(folder_nets,'FPC-bin.nii');
FPC_img = FPC.dat(:,:,:);

SMN = nifti(strcat(folder_nets,'SMN-bin.nii'));
SMN.dat.fname = strcat(folder_nets,'SMN-bin.nii');
SMN_img = SMN.dat(:,:,:);

DMN = nifti(strcat(folder_nets,'DMN-bin.nii'));
DMN.dat.fname = strcat(folder_nets,'DMN-bin.nii');
DMN_img = DMN.dat(:,:,:);

% folder_fun = 'Z:\_DATA\Parcellation\758-Cluster\';
folder_fun = '/Volumes/INDIREA/_DATA/Parcellation/758-Cluster/';

MD = nifti(strcat(folder_fun,'LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii'));
MD.dat.fname = strcat(folder_fun,'LHR-All-Subjects-FC-Voxel-AAL-ROI-KMeans-Parcellation.nii');
MD_img = MD.dat(:,:,:);

nNets = 8;
nROIs = 90;
nTotalRuns = 32;
nTotalVoxels = 160990;
nTotalClusters = 758;
MNI_size = [91 109 91];

for iNet=1:nNets
    
    clusters_per_net(iNet).population_map = zeros(1,nTotalClusters);
    
end

iCluster = 0;
iiCluster = 0;
for iROI=1:nROIs
    
    nClusters = ROI_MD758Corbetta(iROI).nClusters;
    
    for iClusters=1:nClusters
        
        iiCluster = iiCluster + 1;

        idx_voxels = ROI_MD758Corbetta(iROI).clusters(iClusters).idx_voxels;
        
        nVoxels = length(idx_voxels);
        
        for iVoxel=1:nVoxels
            
            [idxx,idxy,idxz] = ind2sub(MNI_size,idx_voxels(iVoxel));
            
            if DAN_img(idxx,idxy,idxz)
                
                iNet = 1;
                clusters_per_net(iNet).population_map(iiCluster) = 1;
            
            end
            
            if VAN_img(idxx,idxy,idxz)
                
                iNet = 2;
                clusters_per_net(iNet).population_map(iiCluster) = 1;
            
            end
            
            if VIS_img(idxx,idxy,idxz)
                
                iNet = 3;
                clusters_per_net(iNet).population_map(iiCluster) = 1;
            
            end
            
            if LAN_img(idxx,idxy,idxz)
                
                iNet = 4;
                clusters_per_net(iNet).population_map(iiCluster) = 1;
            
            end
            
            if FPC_img(idxx,idxy,idxz)
                
                iNet = 5;
                clusters_per_net(iNet).population_map(iiCluster) = 1;
            
            end
            
            if SMN_img(idxx,idxy,idxz)
                
                iNet = 6;
                clusters_per_net(iNet).population_map(iiCluster) = 1;
            
            end
            
            if AUD_img(idxx,idxy,idxz)
                
                iNet = 7;
                clusters_per_net(iNet).population_map(iiCluster) = 1;
            
            end
            
            if DMN_img(idxx,idxy,idxz)
                
                iNet = 8;
                clusters_per_net(iNet).population_map(iiCluster) = 1;
            
            end
            
        end
        
    end
    
end

%%% PLOT POPULATION MAP

Net_Population_Map_Color(1,:) = [0,0,128]./255;
Net_Population_Map_Color(2,:) = [160,32,240]./255;
Net_Population_Map_Color(3,:) = [0,100,0]./255;
Net_Population_Map_Color(6,:) = [0,191,255]./255;
Net_Population_Map_Color(4,:) = [255,165,0]./255;
Net_Population_Map_Color(5,:) = [255,255,0]./255;
Net_Population_Map_Color(7,:) = [192,192,192]./255;
Net_Population_Map_Color(8,:) = [255,0,0]./255;


for iNet=1:nNets
    
    f = figure;
    set(gcf, 'color', 'w');
    colormap(jet);
   
    population_map = clusters_per_net(iNet).population_map;
        
    bar(population_map,'FaceColor',Net_Population_Map_Color(iNet,:));
    hold on
    
%     if iNet == 7
%         
%         marker_color = 'w*';
%         
%     else
%         
%         marker_color = 'k*';
%         
%     end

    marker_color = 'k*';
    
    for iCluster=1:length(clusters_per_net(iNet).clusters)
        
        plot(clusters_per_net(iNet).clusters(iCluster),0.5,marker_color,'MarkerSize',12);
        
    end
    
    xlim([1 nTotalClusters]);
    axis off
    print(f,strcat('MD758-Networks-Overlaps-',networks{iNet},'.eps'),'-depsc'); 
    
end




