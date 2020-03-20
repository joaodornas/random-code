
zthreshold = 2.5;
pvalue = 0.05;
MNI = 'MNI Structural Atlas';

% prefix = strcat('LHR','-','Correlation-Contrast');
% prefix = strcat('LHR');

% correlation_contrast{1} = 'PR-All-z';
% correlation_contrast{2} = 'PR-Neg-z';
% correlation_contrast{3} = 'PR-Pos-z';
% correlation_contrast{4} = 'PT-All-z';
% correlation_contrast{5} = 'PT-Neg-z';
% correlation_contrast{6} = 'PT-Pos-z';
% correlation_contrast{7} = 'RP-All-z';
% correlation_contrast{8} = 'RP-Neg-z';
% correlation_contrast{9} = 'RP-Pos-z';
% correlation_contrast{10} = 'TP-All-z';
% correlation_contrast{11} = 'TP-Neg-z';
% correlation_contrast{12} = 'TP-Pos-z';

% correlation_contrast{1} = 'PR-All-z-clu';
% correlation_contrast{2} = 'PR-Neg-z-clu';
% correlation_contrast{3} = 'PR-Pos-z-clu';
% correlation_contrast{4} = 'PT-All-z-clu';
% correlation_contrast{5} = 'PT-Neg-z-clu';
% correlation_contrast{6} = 'PT-Pos-z-clu';
% correlation_contrast{7} = 'RP-All-z-clu';
% correlation_contrast{8} = 'RP-Neg-z-clu';
% correlation_contrast{9} = 'RP-Pos-z-clu';
% correlation_contrast{10} = 'TP-All-z-clu';
% correlation_contrast{11} = 'TP-Neg-z-clu';
% correlation_contrast{12} = 'TP-Pos-z-clu';

% correlation_contrast{1} = 'PR-All-z-clu-mask';
% correlation_contrast{2} = 'PR-Neg-z-clu-mask';
% correlation_contrast{3} = 'PR-Pos-z-clu-mask';
% correlation_contrast{4} = 'PT-All-z-clu-mask';
% correlation_contrast{5} = 'PT-Neg-z-clu-mask';
% correlation_contrast{6} = 'PT-Pos-z-clu-mask';
% correlation_contrast{7} = 'RP-All-z-clu-mask';
% correlation_contrast{8} = 'RP-Neg-z-clu-mask';
% correlation_contrast{9} = 'RP-Pos-z-clu-mask';
% correlation_contrast{10} = 'TP-All-z-clu-mask';
% correlation_contrast{11} = 'TP-Neg-z-clu-mask';
% correlation_contrast{12} = 'TP-Pos-z-clu-mask';

% correlation_contrast{1} = 'PR-All-z_cluster_info_data';
% correlation_contrast{2} = 'PR-Neg-z_cluster_info_data';
% correlation_contrast{3} = 'PR-Pos-z_cluster_info_data';
% correlation_contrast{4} = 'PT-All-z_cluster_info_data';
% correlation_contrast{5} = 'PT-Neg-z_cluster_info_data';
% correlation_contrast{6} = 'PT-Pos-z_cluster_info_data';
% correlation_contrast{7} = 'RP-All-z_cluster_info_data';
% correlation_contrast{8} = 'RP-Neg-z_cluster_info_data';
% correlation_contrast{9} = 'RP-Pos-z_cluster_info_data';
% correlation_contrast{10} = 'TP-All-z_cluster_info_data';
% correlation_contrast{11} = 'TP-Neg-z_cluster_info_data';
% correlation_contrast{12} = 'TP-Pos-z_cluster_info_data';

correlation_contrast{1} = 'Global-Attention-Neg-Decrease-GN-PerRun-fl-1mm';
correlation_contrast{2} = 'Global-Attention-Neg-Increase-GN-PerRun-fl-1mm';
correlation_contrast{3} = 'Global-Attention-Pos-Decrease-GN-PerRun-fl-1mm';
correlation_contrast{4} = 'Global-Attention-Pos-Increase-GN-PerRun-fl-1mm';
correlation_contrast{5} = 'Global-Stimulus-Neg-Decrease-GN-PerRun-fl-1mm';
correlation_contrast{6} = 'Global-Stimulus-Neg-Increase-GN-PerRun-fl-1mm';
correlation_contrast{7} = 'Global-Stimulus-Pos-Decrease-GN-PerRun-fl-1mm';
correlation_contrast{8} = 'Global-Stimulus-Pos-Increase-GN-PerRun-fl-1mm';

for iseed=1:length(correlation_contrast)
    
     inputimg = strcat(correlation_contrast{iseed});
     outputimg = strcat(inputimg,'-clu');
    %%% --volume=228483
    %%% --dlh=0.15
     system (sprintf('/usr/local/fsl/bin/cluster --in=%s --zthresh=2.5 --pthresh=0.05 --dlh=1.5 --volume=902629 --mm --othresh=%s >%s_cluster_info.txt',inputimg,outputimg,inputimg));

%     inputimg = strcat(prefix,'-',correlation_contrast{iseed});
%     system(sprintf('/usr/local/fsl/bin/fslmaths %s -ptoz %s-z',inputimg,inputimg));

%     inputimg = strcat(prefix,'-',correlation_contrast{iseed});
%     system(sprintf('/usr/local/fsl/bin/fslmaths %s -bin %s-mask',inputimg,inputimg));
    
%     inputimg = strcat(prefix,'-',correlation_contrast{iseed});
%     system(sprintf('/usr/local/fsl/bin/atlasquery -a "Harvard-Oxford Subcortical Structural Atlas" -m %s >%s_query.txt',inputimg,inputimg));

%if iseed ~= 8, inputimg = strcat(prefix,'-',correlation_contrast{iseed},'.txt'); end

% if iseed ~= 8 && iseed ~= 2 && iseed ~= 5, inputimg = strcat(prefix,'-',correlation_contrast{iseed},'.txt'); end
% 

% inputimg = strcat(prefix,'-',correlation_contrast{iseed},'_cluster_info_data.txt');
% 
% if ~isempty(inputimg) 
% 
%     data = csvread(inputimg);
%     nClusters = size(data,1);
%     clusters = cell(nClusters,2);
%     for iCluster=1:nClusters
% 
%         X = data(end-iCluster+1,6);
%         Y = data(end-iCluster+1,7);
%         Z = data(end-iCluster+1,8);
% 
%         AAL_name = getROIName(X,Y,Z);
% 
%         clusters{iCluster,1} = int2str(iCluster);
%         clusters{iCluster,2} = AAL_name;
% 
% 
% 
%     end
% 
%     if ~isempty(clusters)
%         
%         clusters = flipud(clusters);
% 
%         formatSpec = '%s %s\n';
%         fileID = fopen(strcat(inputimg,'_AAL.txt'),'w');
%         for row = 1:nClusters
%             fprintf(fileID,formatSpec,clusters{row,:});
%         end
%         fclose(fileID);
% 
%     end
% 
% end

%     inputimg = strcat(prefix,'-',correlation_contrast{iseed});
%     outputimg = strcat(inputimg,'-clu');
%     %% --volume=228483
%     %% --dlh=0.15
%     system (sprintf('/usr/local/fsl/bin/cluster --in=%s --zthresh=2.5 --pthresh=0.05 --dlh=1.5 --volume=902629 --mm --othresh=%s >%s_cluster_info.txt',inputimg,outputimg,inputimg));


end
