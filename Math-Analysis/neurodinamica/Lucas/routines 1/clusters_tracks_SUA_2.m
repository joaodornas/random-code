load ('population_data','cells','Pref_SF','Pref_TF','SF_BW','TF_BW');
load ('track_names');
nclusters=size(track_names,1)/2;
sf_x=zeros(nclusters,1);
sf_y=zeros(nclusters,1);
tf_x=zeros(nclusters,1);
tf_y=zeros(nclusters,1);
sfbw_x=zeros(nclusters,1);
sfbw_y=zeros(nclusters,1);
tfbw_x=zeros(nclusters,1);
tfbw_y=zeros(nclusters,1);
rho_resp=zeros(nclusters,1);
p_resp=zeros(nclusters,1);

for i=1:2:101
    cluster=[track_names(i,:);track_names(i+1,:)];
    assignin('base',strcat('cluster',int2str((i+1)/2)),cluster);
    for j=1:size(cluster,1)
        unit=cluster(j,:);
        ind=strmatch(unit,cells);
        vector=strcat('responseVector',int2str(ind));
        load('partial_corr_bydata',vector);
        vector=eval(vector);
        vector_norm=vector./max(vector);
        assignin('base',strcat('norm_response_',int2str(i),'_',int2str(j)),vector_norm)
        if j==1
            sf_x((i+1)/2)=Pref_SF(ind);
            tf_x((i+1)/2)=Pref_TF(ind);
            sfbw_x((i+1)/2)=SF_BW(ind);
            tfbw_x((i+1)/2)=TF_BW(ind);
        elseif j==2
            sf_y((i+1)/2)=Pref_SF(ind);
            tf_y((i+1)/2)=Pref_TF(ind);
            sfbw_y((i+1)/2)=SF_BW(ind);
            tfbw_y((i+1)/2)=TF_BW(ind);
        end
    end
    [r,p]=corr(eval(strcat('norm_response_',int2str(i),'_1')),eval(strcat('norm_response_',int2str(i),'_2')),'rows','pairwise','type','Spearman');
    assignin('base',strcat('rho_cluster',int2str(i),'_1x2'),r)
    assignin('base',strcat('p_cluster',int2str(i),'_1x2'),p)
    rho_resp((i+1)/2)=r;
    p_resp((i+1)/2)=p;
end

[r_sf,p_sf]=corr(sf_x,sf_y,'rows','pairwise','type','Spearman');
figure;
plot(sf_x,sf_y);
title 'Pref SF';
[r_tf,p_tf]=corr(tf_x,tf_y,'rows','pairwise','type','Spearman');
figure;
plot(tf_x,tf_y);
title 'Pref TF';
[r_sfbw,p_sfbw]=corr(sfbw_x,sfbw_y,'rows','pairwise','type','Spearman');
figure;
plot(sfbw_x,sfbw_y);
title 'SF-BW';
[r_tfbw,p_tfbw]=corr(tfbw_x,tfbw_y,'rows','pairwise','type','Spearman');
figure;
plot(tfbw_x,tfbw_y);
title 'TF-BW';

mean_rho_resp=mean(rho_resp);
sem_rho_resp=std(rho_resp)/sqrt(51);

%save('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\clusters\clusters_tracks_SUA');