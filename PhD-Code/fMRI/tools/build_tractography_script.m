
load('FC-Voxels-AAL-ROI-corr-KMeans-MeaningFul-FC-GC-Att-Stim-Only.mat');

Hubs = [264 266 511 358 573];

for iHub=Hubs
    
    all_clusters = [];

    for iMean=2:length(MeaningFul_Attention)
    
        if MeaningFul_Attention{iMean,1} == iHub
           
            clusters = MeaningFul_Attention{iMean,7};
            all_clusters = [all_clusters; clusters(:)];
            
            clusters = MeaningFul_Attention{iMean,8};
            all_clusters = [all_clusters; clusters(:)];
            
            clusters = MeaningFul_Attention{iMean,11};
            all_clusters = [all_clusters; clusters(:)];
            
            clusters = MeaningFul_Attention{iMean,12};
            all_clusters = [all_clusters; clusters(:)];
            
            clusters = MeaningFul_Attention{iMean,15};
            all_clusters = [all_clusters; clusters(:)];
            
            clusters = MeaningFul_Attention{iMean,16};
            all_clusters = [all_clusters; clusters(:)];
            
        end
        
    end
    
    prefix = strcat('fslmaths cluster-',int2str(iHub));

    command = sprintf('%s -add cluster-%s',prefix,int2str(all_clusters(1)));

    for iClusters=2:length(all_clusters)

        command = sprintf('%s -add cluster-%s',command,int2str(all_clusters(iClusters)));

    end
    
    command = sprintf('%s cluster-%s-complete',command,int2str(iHub));
        
    fileID = fopen(strcat('fslmaths-hub-',int2str(iHub),'.txt'),'w');
    fprintf(fileID,'%s',command);
    fclose(fileID);
    
    clear command
    
end