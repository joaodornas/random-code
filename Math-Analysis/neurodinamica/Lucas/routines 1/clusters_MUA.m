load('population_data','cells','Pref_SF','Pref_TF','SF_BW','TF_BW');
sf_sua=zeros(48,1);
sf_mua=zeros(48,1);
tf_sua=zeros(48,1);
tf_mua=zeros(48,1);
sfbw_sua=zeros(48,1);
sfbw_mua=zeros(48,1);
tfbw_sua=zeros(48,1);
tfbw_mua=zeros(48,1);
rho_resp=zeros(48,1);
p_resp=zeros(48,1);

for i=1:48
    cluster=eval(strcat('cluster',int2str(i)));
    sua=cluster(1,:);
    ind=strmatch(sua,cells);
    sf_sua(i)=Pref_SF(ind);
    tf_sua(i)=Pref_TF(ind);
    sfbw_sua(i)=SF_BW(ind);
    tfbw_sua(i)=TF_BW(ind);
    vector=strcat('responseVector',int2str(ind));
    load('partial_corr_bydata',vector);
    vector=eval(vector);
    sua_vector_norm=vector./max(vector); 
    assignin('base',strcat('norm_response_sua',int2str(i)),sua_vector_norm);
    recording=cluster(1,1:10);
    channel=cluster(1,11);
    unit=cluster(1,12);
    mua=strcat(recording,channel,'m_min_',unit);
    ind=strmatch(mua,MUA_names);
    sf_mua(i)=MUA_Pref_SF(ind);
    tf_mua(i)=MUA_Pref_TF(ind);
    sfbw_mua(i)=MUA_SF_BW(ind);
    tfbw_mua(i)=MUA_TF_BW(ind);
    mua_path=strcat('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\MUA\',mua);
    load(mua_path,'STTCMeanResponse','Sf_octave','Tf_octave');
    vector=zeros(size(Sf_octave,2)*size(Tf_octave,2),1);
    for l=1:size(Sf_octave,2)
        vector((l-1)*size(Tf_octave,2)+1:(l-1)*size(Tf_octave,2)+size(Tf_octave,2))=STTCMeanResponse(l,:);
    end
    mua_vector_norm=vector./max(vector);
    assignin('base',strcat('norm_response_mua',int2str(i)),mua_vector_norm);
    [r,p]=corr(sua_vector_norm,mua_vector_norm,'rows','pairwise','type','Spearman');
    assignin('base',strcat('rho_cluster',int2str(i)),r)
    assignin('base',strcat('p_cluster',int2str(i)),p)
    rho_resp(i)=r;
    p_resp(i)=p;
    clear(strcat('responseVector',int2str(ind)),'vector','STTCMeanResponse','Sf_octave','Tf_octave','recording','mua_path','cluster','r','p','unit','channel','sua','mua','norm_response_mua','norm_response_sua','ind');
end
    
[r_sf,p_sf]=corr(sf_sua,sf_mua,'rows','pairwise','type','Spearman');
figure;
plot(sf_sua,sf_mua);
title 'SUA x MUA: PREF SF'
[r_sfbw,p_sfbw]=corr(sfbw_sua,sfbw_mua,'rows','pairwise','type','Spearman');
figure;
plot(sfbw_sua,sfbw_mua);
title 'SUA x MUA: SF-BW'
[r_tf,p_tf]=corr(tf_sua,tf_mua,'rows','pairwise','type','Spearman');
figure;
plot(tf_sua,tf_mua);
title 'SUA x MUA: PREF TF'
[r_tfbw,p_tfbw]=corr(tfbw_sua,tfbw_mua,'rows','pairwise','type','Spearman');
figure;
plot(tfbw_sua,tfbw_mua);
title 'SUA x MUA: TF-BW'

save('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\clusters\clusters_MUA');
