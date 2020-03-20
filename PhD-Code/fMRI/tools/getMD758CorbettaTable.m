
clear all

load('/Volumes/INDIREA/_DATA/Parcellation/758-Cluster/Corbetta/ROI_MD758Corbetta.mat');

networks = {'DAN' 'VAN' 'SMN' 'VIS' 'FPC' 'LAN' 'DMN' 'AUD'};

nROI = 90;
iiCluster = 0;
iiCor = 0;
for iROI=1:nROI
    
    nClusters = ROI_MD758Corbetta(iROI).nClusters;
    
    for iCluster=1:nClusters
        
        iiCluster = iiCluster + 1;
        
        if isfield(ROI_MD758Corbetta(iROI).clusters(iCluster),'Nets')
        
            if ~isempty(ROI_MD758Corbetta(iROI).clusters(iCluster).Nets)

                iiCor = iiCor + 1;

                nets = ROI_MD758Corbetta(iROI).clusters(iCluster).Nets;

                nNets = length(nets);

                correspondence{iiCor,1} = ROI_MD758Corbetta(iROI).label;
                correspondence{iiCor,2} = iiCluster;
                
                % correspondence{iiCor,3} = strjoin(nets);
                
                for iNet=1:nNets
                    
                    this_net = nets{iNet};
                    
                    the_net = this_net(1:3);
                    the_seed = this_net(5:end);
                    
                    correspondence{iiCor,2+find(strcmp(networks,the_net))} = the_seed;
                    
                end

            end
        
        end
        
    end
    
end

