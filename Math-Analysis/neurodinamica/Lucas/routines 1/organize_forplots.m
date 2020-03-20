load('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\1D Fits\Spatial Frequency\Matlab workspaces\SF_stt047a02_1b');
SF_fitx1b=analysisresults1.xi;
SF_fity1b=analysisresults1.yfit;
load('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\1D Fits\Temporal Frequency\Matlab workspaces\TF_stt047a02_1b','analysisresults1');
TF_fitx1b=analysisresults1.xi;
TF_fity1b=analysisresults1.yfit;
load('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\Multiple1DFits\workspaces\multiple_stt047a02_1b','STTCMeanResponse');
SF_SEM1b=zeros(1,6);
TF_SEM1b=zeros(1,6);
meanRates_TemporalFreq_minBase1b=meanRates_TemporalFreq_minBase;
meanRates_SpatialFreq_minBase1b=meanRates_SpatialFreq_minBase;
mean_SF_values1b=mean_SF_Values;
mean_TF_values1b=mean_TF_Values;
twoSEM1b=twoSEM;
sfwonan=trialRates_SpatialFreq_minBase;
tfwonan=trialRates_TemporalFreq_minBase;
for i=1:66
    t=isnan(trialRates_SpatialFreq_minBase(i));
    if t==1
        sfwonan(i)=-1;
    end
    t=isnan(trialRates_TemporalFreq_minBase(i));
    if t==1
        tfwonan(i)=-1;
    end
end
for i=1:6
    a=find(sfwonan((i-1)*11+1:(i-1)*11+11)>=0);
    ntrials=size(a,2);
    SF_SEM1b(i)=nanstd(trialRates_SpatialFreq_minBase((i-1)*11+1:(i-1)*11+11))/sqrt(ntrials);
    a=find(tfwonan((i-1)*11+1:(i-1)*11+11)>=0);
    ntrials=size(a,2);
    TF_SEM1b(i)=nanstd(trialRates_TemporalFreq_minBase((i-1)*11+1:(i-1)*11+11))/sqrt(ntrials);
end