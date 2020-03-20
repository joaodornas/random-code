nclusters = size (rho_resp,1);
recordSiteDist = zeros (nclusters,1);

for i=1:nclusters
    cell1=eval(['cluster' int2str(i) '{1}']);
    cell2=eval(['cluster' int2str(i) '{2}']);
    site1=[cell1(1:7) cell1(10:11)];
    site2=[cell2(1:7) cell2(10:11)];
    a1=strmatch(site1,recordingSites);
    a2=strmatch(site2,recordingSites);
    recordSiteDist(i)=abs(recordingDepths(a1)-recordingDepths(a2));
end

[r_distPerRho,p_distPerRho]=corr(rho_resp,recordSiteDist,'rows','pairwise','type','spearman');