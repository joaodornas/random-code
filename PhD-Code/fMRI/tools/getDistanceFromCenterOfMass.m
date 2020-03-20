function distance_cluster = getDistanceFromCenterOfMass(X_ROI,Y_ROI,Z_ROI,X_cluster,Y_cluster,Z_cluster)

nCluster = 10;
nRun = 4;

for irun=1:nRun
    
   for iCluster=1:nCluster
       
      distance_cluster(irun,iCluster) = sqrt((X_cluster(irun,iCluster) - X_ROI)^2 + (Y_cluster(irun,iCluster) - Y_ROI)^2 + (Z_cluster(irun,iCluster) - Z_ROI)^2); 
       
   end
    
end


end