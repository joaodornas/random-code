function real_FC_voxel_AAL_ROI_center_of_mass(settings,idx_ROI)

%%% LOAD AAL

disp('Load AAL');

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nRuns = 4;
nClusters = 10;

for iROI=1:length(idx_ROI)

    label_ROI = AAL_ROI(idx_ROI(iROI)).Nom_L;
    idx_AAL = AAL_ROI(idx_ROI(iROI)).ID;

    area_label = strrep(label_ROI,'_','-');       

    disp(area_label);

    idx_voxels = find(AAL_img==idx_AAL);
    [center_of_X_ROI, center_of_Y_ROI, center_of_Z_ROI] = getCenterOfMass(idx_voxels,AAL_img);

    TrackKmeans = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Track','-',area_label,'-','KMeans','.mat'));
    PassiveKmeans = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','Passive','-',area_label,'-','KMeans','.mat'));
    RestingStateKmeans = load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-','RestingState','-',area_label,'-','KMeans','.mat'));

    for irun=1:nRuns

       track_clusters = TrackKmeans.FC_KMeans.run(irun).IdxClusters;
       passive_clusters = PassiveKmeans.FC_KMeans.run(irun).IdxClusters;
       restingstate_clusters = RestingStateKmeans.FC_KMeans.run(irun).IdxClusters;

       for iCluster=1:nClusters

            idx_voxels_track = find(track_clusters==iCluster);
            idx_voxels_track = idx_voxels(idx_voxels_track);
            [center_of_X_track(irun,iCluster), center_of_Y_track(irun,iCluster), center_of_Z_track(irun,iCluster)] = getCenterOfMass(idx_voxels_track,AAL_img);

            idx_voxels_passive = find(passive_clusters==iCluster);
            idx_voxels_passive = idx_voxels(idx_voxels_passive);
            [center_of_X_passive(irun,iCluster), center_of_Y_passive(irun,iCluster), center_of_Z_passive(irun,iCluster)] = getCenterOfMass(idx_voxels_passive,AAL_img);

            idx_voxels_rest = find(restingstate_clusters==iCluster);
            idx_voxels_rest = idx_voxels(idx_voxels_rest);
            [center_of_X_rest(irun,iCluster), center_of_Y_rest(irun,iCluster), center_of_Z_rest(irun,iCluster)] = getCenterOfMass(idx_voxels_rest,AAL_img);

       end

    end

    distance_cluster_track = getDistanceFromCenterOfMass(center_of_X_ROI,center_of_Y_ROI,center_of_Z_ROI,center_of_X_track,center_of_Y_track,center_of_Z_track); 
    distance_cluster_passive = getDistanceFromCenterOfMass(center_of_X_ROI,center_of_Y_ROI,center_of_Z_ROI,center_of_X_passive,center_of_Y_passive,center_of_Z_passive); 
    distance_cluster_rest = getDistanceFromCenterOfMass(center_of_X_ROI,center_of_Y_ROI,center_of_Z_ROI,center_of_X_rest,center_of_Y_rest,center_of_Z_rest); 


    f = figure;

    for irun=1:nRuns

       plot(sort(distance_cluster_track(irun,:)),'b');
       hold on
       plot(sort(distance_cluster_passive(irun,:)),'r');
       plot(sort(distance_cluster_rest(irun,:)),'k');

    end

    print(f,'-djpeg',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-',area_label,'-Distance-Center-Of-Mass','.jpeg'));
    print(f,'-depsc',strcat(settings.codes.experiment,'-',settings.codes.subject,'-','FC-Voxels-AAL-ROI','-',area_label,'-Distance-Center-Of-Mass','.eps'));

end

end

function [center_of_X, center_of_Y, center_of_Z] = getCenterOfMass(idx,AAL_img)

nIndices = length(idx);

for iidx=1:nIndices
    
    [X(iidx),Y(iidx),Z(iidx)] = ind2sub(size(AAL_img),idx(iidx));
    
end

center_of_X = round(mean(X));
center_of_Y = round(mean(Y));
center_of_Z = round(mean(Z));

end

function distance_cluster = getDistanceFromCenterOfMass(X_ROI,Y_ROI,Z_ROI,X_cluster,Y_cluster,Z_cluster)

nCluster = 10;
nRun = 4;

for irun=1:nRun
    
   for iCluster=1:nCluster
       
      distance_cluster(irun,iCluster) = sqrt((X_cluster(irun,iCluster) - X_ROI)^2 + (Y_cluster(irun,iCluster) - Y_ROI)^2 + (Z_cluster(irun,iCluster) - Z_ROI)^2); 
       
   end
    
end


end
