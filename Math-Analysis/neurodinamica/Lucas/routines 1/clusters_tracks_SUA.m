load ('population_data','cells','Pref_SF','Pref_TF','SF_BW','TF_BW');
sf_x=[];
sf_y=[];
tf_x=[];
tf_y=[];
sfbw_x=[];
sfbw_y=[];
tfbw_x=[];
tfbw_y=[];
rho_resp=[];
p_resp=[];

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
        if size(cluster,1)==2
            if j==1
                a=find(sf_x,1,'last');
                t=isempty(a);
                if t==1
                sf_x(1,1)=Pref_SF(ind);
                tf_x(1,1)=Pref_TF(ind);
                sfbw_x(1,1)=SF_BW(ind);
                tfbw_x(1,1)=TF_BW(ind);
                else sf_x(a+1,1)=Pref_SF(ind);
                    tf_x(a+1,1)=Pref_TF(ind);
                    sfbw_x(a+1,1)=SF_BW(ind);
                    tfbw_x(a+1,1)=TF_BW(ind);
                end
            elseif j==2
                a=find(sf_y,1,'last');
                t=isempty(a);
                if t==1
                sf_y(1,1)=Pref_SF(ind);
                tf_y(1,1)=Pref_TF(ind);
                sfbw_y(1,1)=SF_BW(ind);
                tfbw_y(1,1)=TF_BW(ind);
                else sf_y(a+1,1)=Pref_SF(ind);
                    tf_y(a+1,1)=Pref_TF(ind);
                    sfbw_y(a+1,1)=SF_BW(ind);
                    tfbw_y(a+1,1)=TF_BW(ind);
                end
            end
        elseif size(cluster,1)==3
            if j==1
                a=find(sf_x,1,'last');
                t=isempty(a);
                if t==1
                sf_x(1:2,1)=Pref_SF(ind);
                tf_x(1:2,1)=Pref_TF(ind);
                sfbw_x(1:2,1)=SF_BW(ind);
                tfbw_x(1:2,1)=TF_BW(ind);
                else sf_x(a+1:a+2,1)=Pref_SF(ind);
                    tf_x(a+1:a+2,1)=Pref_TF(ind);
                    sfbw_x(a+1:a+2,1)=SF_BW(ind);
                    tfbw_x(a+1:a+2,1)=TF_BW(ind);
                end
            elseif j==2
                a=find(sf_x,1,'last');
                t=isempty(a);
                if t==1
                sf_x(1,1)=Pref_SF(ind);
                tf_x(1,1)=Pref_TF(ind);
                sfbw_x(1,1)=SF_BW(ind);
                tfbw_x(1,1)=SF_BW(ind);
                else sf_x(a+1,1)=Pref_SF(ind);
                    tf_x(a+1,1)=Pref_TF(ind);
                    sfbw_x(a+1,1)=SF_BW(ind);
                    tfbw_x(a+1,1)=TF_BW(ind);
                end
                a=find(sf_y,1,'last');
                t=isempty(a);
                if t==1
                sf_y(1,1)=Pref_SF(ind);
                tf_y(1,1)=Pref_TF(ind);
                sfbw_y(1,1)=SF_BW(ind);
                tfbw_y(1,1)=TF_BW(ind);
                else sf_y(a+1,1)=Pref_SF(ind);
                    tf_y(a+1,1)=Pref_TF(ind);
                    sfbw_y(a+1,1)=SF_BW(ind);
                    tfbw_y(1,1)=TF_BW(ind);
                end
            elseif j==3
                a=find(sf_y,1,'last');
                sf_y(a+1:a+2,1)=Pref_SF(ind);
                tf_y(a+1:a+2,1)=Pref_TF(ind);
                sfbw_y(a+1:a+2,1)=SF_BW(ind);
                tfbw_y(a+1:a+2,1)=TF_BW(ind);
            end
        end
    end
    if size(cluster,1)==2
        [r,p]=corr(eval(strcat('norm_response_',int2str(i),'_1')),eval(strcat('norm_response_',int2str(i),'_2')),'rows','pairwise','type','Spearman');
        assignin('base',strcat('rho_cluster',int2str(i),'_1x2'),r)
        assignin('base',strcat('p_cluster',int2str(i),'_1x2'),p)
        a=find(rho_resp,1,'last');
        t=isempty(a);
        if t==1
            rho_resp(1,1)=r;
            p_resp(1,1)=p;
        else rho_resp(a+1,1)=r;
            p_resp(a+1,1)=p;
        end
    elseif size(cluster,1)==3
        [r,p]=corr(eval(strcat('norm_response_',int2str(i),'_1')),eval(strcat('norm_response_',int2str(i),'_2')),'rows','pairwise','type','Spearman');
        assignin('base',strcat('rho_cluster',int2str(i),'_1x2'),r)
        assignin('base',strcat('p_cluster',int2str(i),'_1x2'),p)
        a=find(rho_resp,1,'last');
        t=isempty(a);
        if t==1
            rho_resp(1,1)=r;
            p_resp(1,1)=p;
        else rho_resp(a+1,1)=r;
            p_resp(a+1,1)=p;
        end
        [r,p]=corr(eval(strcat('norm_response_',int2str(i),'_1')),eval(strcat('norm_response_',int2str(i),'_3')),'rows','pairwise','type','Spearman');
        assignin('base',strcat('rho_cluster',int2str(i),'_1x3'),r)
        assignin('base',strcat('p_cluster',int2str(i),'_1x3'),p)
        a=find(rho_resp,1,'last');
        t=isempty(a);
        if t==1
            rho_resp(1,1)=r;
            p_resp(1,1)=p;
        else rho_resp(a+1,1)=r;
            p_resp(a+1,1)=p;
        end
        [r,p]=corr(eval(strcat('norm_response_',int2str(i),'_2')),eval(strcat('norm_response_',int2str(i),'_3')),'rows','pairwise','type','Spearman');
        assignin('base',strcat('rho_cluster',int2str(i),'_2x3'),r)
        assignin('base',strcat('p_cluster',int2str(i),'_2x3'),p)
        a=find(rho_resp,1,'last');
        t=isempty(a);
        if t==1
            rho_resp(1,1)=r;
            p_resp(1,1)=p;
        else rho_resp(a+1,1)=r;
            p_resp(a+1,1)=p;
        end
    end
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