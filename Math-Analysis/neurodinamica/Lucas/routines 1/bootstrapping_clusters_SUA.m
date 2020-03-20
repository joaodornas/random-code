sf_x_SUA_differentSite=zeros(32,1);
sf_y_SUA_differentSite=zeros(32,1);
sfbw_x_SUA_differentSite=zeros(32,1);
sfbw_y_SUA_differentSite=zeros(32,1);
tf_x_SUA_differentSite=zeros(32,1);
tf_y_SUA_differentSite=zeros(32,1);
tfbw_x_SUA_differentSite=zeros(32,1);
tfbw_y_SUA_differentSite=zeros(32,1);

for i = 1:32
    cellname_x=eval(['cellpairs_SUA_differentSites(' num2str(i) ',1);']);
    cellname_y=eval(['cellpairs_SUA_differentSites(' num2str(i) ',2);']);
    i_x=strmatch(cellname_x,cells);
    i_y=strmatch(cellname_y,cells);
    sf_x_SUA_differentSite(i)=Pref_SF(i_x);
    sf_y_SUA_differentSite(i)=Pref_SF(i_y);
    sfbw_x_SUA_differentSite(i)=SF_BW(i_x);
    sfbw_y_SUA_differentSite(i)=SF_BW(i_y);
    tf_x_SUA_differentSite(i)=Pref_TF(i_x);
    tf_y_SUA_differentSite(i)=Pref_TF(i_y);
    tfbw_x_SUA_differentSite(i)=TF_BW(i_x);
    tfbw_y_SUA_differentSite(i)=TF_BW(i_y);
end

boot_i_differentSite=ceil(32.*rand(20,1000));
boot_i_sameSite=ceil(48.*rand(20,1000));

r_Vector_sf_differentSite=nan(1000,1);
p_Vector_sf_differentSite=nan(1000,1);
r_Vector_sf_sameSite=nan(1000,1);
p_Vector_sf_sameSite=nan(1000,1);
r_Vector_tf_differentSite=nan(1000,1);
p_Vector_tf_differentSite=nan(1000,1);
r_Vector_tf_sameSite=nan(1000,1);
p_Vector_tf_sameSite=nan(1000,1);
r_Vector_sfbw_differentSite=nan(1000,1);
p_Vector_sfbw_differentSite=nan(1000,1);
r_Vector_sfbw_sameSite=nan(1000,1);
p_Vector_sfbw_sameSite=nan(1000,1);
r_Vector_tfbw_differentSite=nan(1000,1);
p_Vector_tfbw_differentSite=nan(1000,1);
r_Vector_tfbw_sameSite=nan(1000,1);
p_Vector_tfbw_sameSite=nan(1000,1);

for k=1:1000
    
    sf_x_same=sf_x_SUA_sameSite(boot_i_sameSite(:,k));
    sf_y_same=sf_y_SUA_sameSite(boot_i_sameSite(:,k));
    sfbw_x_same=sfbw_x_SUA_sameSite(boot_i_sameSite(:,k));
    sfbw_y_same=sfbw_y_SUA_sameSite(boot_i_sameSite(:,k));
    tf_x_same=tf_x_SUA_sameSite(boot_i_sameSite(:,k));
    tf_y_same=tf_y_SUA_sameSite(boot_i_sameSite(:,k));
    tfbw_x_same=tfbw_x_SUA_sameSite(boot_i_sameSite(:,k));
    tfbw_y_same=tfbw_y_SUA_sameSite(boot_i_sameSite(:,k));
    
    [r_sf_same,p_sf_same]=corr(sf_x_same,sf_y_same,'rows','pairwise','type','Spearman');
    [r_sfbw_same,p_sfbw_same]=corr(sfbw_x_same,sfbw_y_same,'rows','pairwise','type','Spearman');
    [r_tf_same,p_tf_same]=corr(tf_x_same,tf_y_same,'rows','pairwise','type','Spearman');
    [r_tfbw_same,p_tfbw_same]=corr(tfbw_x_same,tfbw_y_same,'rows','pairwise','type','Spearman');
    
    r_Vector_sf_sameSite(k)=r_sf_same;
    p_Vector_sf_sameSite(k)=p_sf_same;
    r_Vector_sfbw_sameSite(k)=r_sfbw_same;
    p_Vector_sfbw_sameSite(k)=p_sfbw_same;
    r_Vector_tf_sameSite(k)=r_tf_same;
    p_Vector_tf_sameSite(k)=p_tf_same;
    r_Vector_tfbw_sameSite(k)=r_tfbw_same;
    p_Vector_tfbw_sameSite(k)=p_tfbw_same;
    
    sf_x_different=sf_x_SUA_differentSite(boot_i_differentSite(:,k));
    sf_y_different=sf_y_SUA_differentSite(boot_i_differentSite(:,k));
    sfbw_x_different=sfbw_x_SUA_differentSite(boot_i_differentSite(:,k));
    sfbw_y_different=sfbw_y_SUA_differentSite(boot_i_differentSite(:,k));
    tf_x_different=tf_x_SUA_differentSite(boot_i_differentSite(:,k));
    tf_y_different=tf_y_SUA_differentSite(boot_i_differentSite(:,k));
    tfbw_x_different=tfbw_x_SUA_differentSite(boot_i_differentSite(:,k));
    tfbw_y_different=tfbw_y_SUA_differentSite(boot_i_differentSite(:,k));
    
    [r_sf_different,p_sf_different]=corr(sf_x_different,sf_y_different,'rows','pairwise','type','Spearman');
    [r_sfbw_different,p_sfbw_different]=corr(sfbw_x_different,sfbw_y_different,'rows','pairwise','type','Spearman');
    [r_tf_different,p_tf_different]=corr(tf_x_different,tf_y_different,'rows','pairwise','type','Spearman');
    [r_tfbw_different,p_tfbw_different]=corr(tfbw_x_different,tfbw_y_different,'rows','pairwise','type','Spearman');
    
    r_Vector_sf_differentSite(k)=r_sf_different;
    p_Vector_sf_differentSite(k)=p_sf_different;
    r_Vector_sfbw_differentSite(k)=r_sfbw_different;
    p_Vector_sfbw_differentSite(k)=p_sfbw_different;
    r_Vector_tf_differentSite(k)=r_tf_different;
    p_Vector_tf_differentSite(k)=p_tf_different;
    r_Vector_tfbw_differentSite(k)=r_tfbw_different;
    p_Vector_tfbw_differentSite(k)=p_tfbw_different;
    
end

for i=1:105
 
        load('partial_corr_bydata',['responseVector' num2str(i)]);
        vector=eval(['responseVector' int2str(i)]);
        vector=vector./max(vector);
        eval(['resp_vector_norm.cell' num2str(i) '=vector;']);
        
end

rho_resp_SUA_differentSite=nan(32,1);
p_resp_SUA_differentSite=nan(32,1);

for i = 1:32
    cellname_x=eval(['cellpairs_SUA_differentSites(' num2str(i) ',1);']);
    cellname_y=eval(['cellpairs_SUA_differentSites(' num2str(i) ',2);']);
    i_x=strmatch(cellname_x,cells);
    i_y=strmatch(cellname_y,cells);
    vector_x=eval(['resp_vector_norm.cell' num2str(i_x) ';']);
    vector_y=eval(['resp_vector_norm.cell' num2str(i_y) ';']);
    [r,p]=corr(vector_x,vector_y,'rows','pairwise','type','Spearman');
    rho_resp_SUA_differentSite(i)=r;
    p_resp_SUA_differentSite(i)=p;
    clear cellname_x cellname_y i_x i_y vector_x vector_y r p
end

mean_rho_vector_resp_SUA_sameSite=nan(1000,1);
mean_p_vector_resp_SUA_sameSite=nan(1000,1);
percent_signif_mean_rho_vector_resp_SUA_sameSite=nan(1000,1);
mean_rho_vector_resp_SUA_differentSite=nan(1000,1);
mean_p_vector_resp_SUA_differentSite=nan(1000,1);
percent_signif_mean_rho_vector_resp_SUA_differentSite=nan(1000,1);

for k = 1:1000
    rho_vector_same=rho_resp_SUA_sameSite(boot_i_sameSite(:,k));
    p_vector_same=rho_resp_SUA_sameSite(boot_i_sameSite(:,k));
    rho_vector_different=rho_resp_SUA_differentSite(boot_i_differentSite(:,k));
    p_vector_different=rho_resp_SUA_differentSite(boot_i_differentSite(:,k));
    mean_rho_vector_resp_SUA_sameSite(k)=mean(rho_vector_same);
    mean_p_vector_resp_SUA_sameSite(k)=mean(p_vector_same);
    percent_signif_mean_rho_vector_resp_SUA_sameSite(k)=size(find(p_vector_same<0.5),1)*5;
    mean_rho_vector_resp_SUA_differentSite(k)=mean(rho_vector_different);
    mean_p_vector_resp_SUA_differentSite(k)=mean(p_vector_different);
    percent_signif_mean_rho_vector_resp_SUA_differentSite(k)=size(find(p_vector_different<0.5),1)*5;
end

mean_mean_rho_vector_resp_SUA_sameSite=nanmean(mean_rho_vector_resp_SUA_sameSite)
sem_mean_rho_vector_resp_SUA_sameSite=nanstd(mean_rho_vector_resp_SUA_sameSite)/sqrt(1000)
mean_mean_p_vector_resp_SUA_sameSite=nanmean(mean_p_vector_resp_SUA_sameSite)
sem_mean_p_vector_resp_SUA_sameSite=nanstd(mean_p_vector_resp_SUA_sameSite)/sqrt(1000)
mean_percent_signif_mean_rho_vector_resp_SUA_sameSite=nanmean(percent_signif_mean_rho_vector_resp_SUA_sameSite)
sem_percent_signif_mean_rho_vector_resp_SUA_sameSite=nanstd(percent_signif_mean_rho_vector_resp_SUA_sameSite)/sqrt(1000)
mean_mean_rho_vector_resp_SUA_differentSite=nanmean(mean_rho_vector_resp_SUA_differentSite)
sem_mean_rho_vector_resp_SUA_differentSite=nanstd(mean_rho_vector_resp_SUA_differentSite)/sqrt(1000)
mean_mean_p_vector_resp_SUA_differentSite=nanmean(mean_p_vector_resp_SUA_differentSite)
sem_mean_p_vector_resp_SUA_differentSite=nanstd(mean_p_vector_resp_SUA_differentSite)/sqrt(1000)
mean_percent_signif_mean_rho_vector_resp_SUA_differentSite=nanmean(percent_signif_mean_rho_vector_resp_SUA_differentSite)
sem_percent_signif_mean_rho_vector_resp_SUA_differentSite=nanstd(percent_signif_mean_rho_vector_resp_SUA_differentSite)/sqrt(1000)

pvalue_mean_mean_rho_vector_resp_SUA=ranksum(mean_mean_rho_vector_resp_SUA_sameSite,mean_mean_rho_vector_resp_SUA_differentSite)
pvalue_mean_percent_signif_mean_rho_vector_resp_SUA=ranksum(mean_percent_signif_mean_rho_vector_resp_SUA_sameSite,mean_percent_signif_mean_rho_vector_resp_SUA_differentSite)