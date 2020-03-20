% mean_rates_250_sf=nan(size(lat_SFs_atBestTF));
% mean_rates_250_tf=nan(size(lat_TFs_atBestSF));
% 
% for i=[1:4,6:12,14:92]
% 
%     Sf_octave=eval(['SFvalues_octave.complexCell' num2str(i) ';']);
%     Tf_octave=eval(['TFvalues_octave.complexCell' num2str(i) ';']);
%     
%     SFs_i=eval(['SFs_i_atBestTF.complexCell' num2str(i) ';']);
%     TFs_i=eval(['TFs_i_atBestSF.complexCell' num2str(i) ';']);
% 
%     sf_column_i=find(SF_columnLabels_octv==Sf_octave(1));
%     tf_column_i=find(TF_columnLabels_octv==Tf_octave(1));
%     
%     nSFs=numel(Sf_octave); nTFs=numel(Tf_octave);
%     
%     load(['psth_' char(complex_cells(i)) '_convolved_sigma3.mat'],'psth');
% 
%     mean_rates_sf_i=nan(1,nSFs); mean_rates_tf_i=nan(1,nTFs);
%     
%     for k=1:nSFs
%        psth_i=eval(['psth.condition' num2str(SFs_i(k)) ';']);
%        mean_rates_sf_i(k)=mean(psth_i(1001:1250));
%     end
%     
%     for k=1:nTFs
%        psth_i=eval(['psth.condition' num2str(TFs_i(k)) ';']);
%        mean_rates_tf_i(k)=mean(psth_i(1001:1250));
%     end
%     
%     mean_rates_250_sf(i,sf_column_i:sf_column_i+nSFs-1)=mean_rates_sf_i;
%     mean_rates_250_tf(i,tf_column_i:tf_column_i+nTFs-1)=mean_rates_tf_i;
% 
% end
% 
% clear psth Sf_octvae Tf_octave sf_column_i tf_column_i i k nTFs nSFs ...
%     mean_rates_sf_i mean_rates_tf_i psth_i
% 
% average_meanRates_SF=nanmean(mean_rates_250_sf);
% average_meanRates_TF=nanmean(mean_rates_250_tf);
% 
% figure; plot(average_meanRates_SF); figure; plot(average_meanRates_TF);
% 
% lat_vector=nan(92*8,1);
% rate_vector=nan(92*8,1);
% for i=1:92
%     lat_vector((i-1)*8+1:i*8,1)=lat_SFs_atBestTF(i,:);
%     rate_vector((i-1)*8+1:i*8,1)=mean_rates_250_sf(i,:);
% end
% [r_rateXlat,p_rateVSlat]=corr(lat_vector,rate_vector,'rows','pairwise','type','spearman')

peak_amps_sf=nan(size(lat_SFs_atBestTF));
peak_amps_tf=nan(size(lat_TFs_atBestSF));

for i=[1:4,6:12,14:92]

    Sf_octave=eval(['SFvalues_octave.complexCell' num2str(i) ';']);
    Tf_octave=eval(['TFvalues_octave.complexCell' num2str(i) ';']);
    
    SFs_i=eval(['SFs_i_atBestTF.complexCell' num2str(i) ';']);
    TFs_i=eval(['TFs_i_atBestSF.complexCell' num2str(i) ';']);

    sf_column_i=find(SF_columnLabels_octv==Sf_octave(1));
    tf_column_i=find(TF_columnLabels_octv==Tf_octave(1));
    
    nSFs=numel(Sf_octave); nTFs=numel(Tf_octave);
    
    load(['psth_' char(complex_cells(i)) '_convolved_sigma3.mat'],'psth');

    peak_amps_sf_i=nan(1,nSFs); peak_amps_tf_i=nan(1,nTFs);
    
    for k=1:nSFs
       psth_i=eval(['psth.condition' num2str(SFs_i(k)) ';']);
       peak_amps_sf_i(k)=psth_i(find(psth_i(1001:1200)==max(psth_i(1001:1200)),1,'first'));
    end
    
    for k=1:nTFs
       psth_i=eval(['psth.condition' num2str(TFs_i(k)) ';']);
       peak_amps_tf_i(k)=psth_i(find(psth_i(1001:1200)==max(psth_i(1001:1200)),1,'first'));
    end
    
    peak_amps_sf(i,sf_column_i:sf_column_i+nSFs-1)=peak_amps_sf_i;
    peak_amps_tf(i,tf_column_i:tf_column_i+nTFs-1)=peak_amps_tf_i;

end

clear psth Sf_octave Tf_octave sf_column_i tf_column_i i k nTFs nSFs ...
    mean_rates_sf_i mean_rates_tf_i psth_i