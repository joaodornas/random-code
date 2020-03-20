load ('latForEachFreq','complex_cells','SFs_i_atBestTF','SFvalues_octave','TFs_i_atBestSF','TFvalues_octave','lat_SFs_atBestTF','SF_columnLabels_octv','lat_TFs_atBestSF','TF_columnLabels_octv');
load('exp_fit_psth','tau_vector','i_peakPsth');
ncells=numel(complex_cells);
baseline_duration=1000;

pref_SF_tillTau=nan(ncells,1);
pref_SF_afterTau=nan(ncells,1);
pref_SF_200=nan(ncells,1);
pref_SF_800=nan(ncells,1);
pref_TF_tillTau=nan(ncells,1);
pref_TF_afterTau=nan(ncells,1);
pref_TF_200=nan(ncells,1);
pref_TF_800=nan(ncells,1);
sigma_SF_tillTau=nan(ncells,1);
sigma_SF_afterTau=nan(ncells,1);
sigma_SF_200=nan(ncells,1);
sigma_SF_800=nan(ncells,1);
sigma_TF_tillTau=nan(ncells,1);
sigma_TF_afterTau=nan(ncells,1);
sigma_TF_200=nan(ncells,1);
sigma_TF_800=nan(ncells,1);
fitp_SF_tillTau=nan(ncells,1);
fitp_SF_afterTau=nan(ncells,1);
fitp_SF_200=nan(ncells,1);
fitp_SF_800=nan(ncells,1);
fitp_TF_tillTau=nan(ncells,1);
fitp_TF_afterTau=nan(ncells,1);
fitp_TF_200=nan(ncells,1);
fitp_TF_800=nan(ncells,1);
rsquare_SF_tillTau=nan(ncells,1);
rsquare_SF_afterTau=nan(ncells,1);
rsquare_SF_200=nan(ncells,1);
rsquare_SF_800=nan(ncells,1);
rsquare_TF_tillTau=nan(ncells,1);
rsquare_TF_afterTau=nan(ncells,1);
rsquare_TF_200=nan(ncells,1);
rsquare_TF_800=nan(ncells,1);

priebe_1D=fittype('A.*((exp((-((log2(x)-log2(ps))^2)./(2.*(tuningWidth+skew.*(log2(x)-log2(ps))).^2))))-(exp(-1./((skew).^2))))');

for i=[1:4,6:12,14:ncells]
    
    tic
    
    disp(['fitting cell #' num2str(i)])
    
    %%%%%%%%%% SF %%%%%%%%%%%%%%%%%%
    
    %get SF conditions and values for cell i
    eval(['SFs_i=SFs_i_atBestTF.complexCell' num2str(i) ';']);
    eval(['SFvalues_octv=SFvalues_octave.complexCell' num2str(i) ';']);

    SFvalues_lin=2.^SFvalues_octv;
    eval(['SFvalues_linScale.complexCell' num2str(i) '=SFvalues_lin;']);
    
    %load and shape psths in matrices where each row is a condition
    load(['psth_' char(complex_cells(i)) '_convolved_sigma3.mat'],'psth');

    SF_psths=zeros(numel(SFs_i),numel(psth.condition1));
    for k=1:numel(SFs_i)
        SF_psths(k,:)=eval(['psth.condition' num2str(SFs_i(k))]);
        mean_baseline=mean(SF_psths(k,1:baseline_duration));
        SF_psths(k,:)=SF_psths(k,:)-mean_baseline;
    end
    
    %latencies for each condition, if nan, then lat=55 (mean latency)
    sf_column_i=find(SF_columnLabels_octv==SFvalues_octv(1));
    lat=lat_SFs_atBestTF(i,sf_column_i:sf_column_i+numel(SFvalues_octv)-1)';
    lat(find(isnan(lat)))=55;
    sf_responses_tillTau=zeros(size(lat));
    sf_responses_afterTau=zeros(size(lat));
    sf_responses_200=zeros(size(lat));
    sf_responses_800=zeros(size(lat));
    
    for j=1:numel(lat)
        if isnan(i_peakPsth(i))
            i_peak=nanmean(i_peakPsth(i));
        else i_peak=i_peakPsth(i);
        end
        if isnan(tau_vector(i)) || tau_vector(i)>1800
            tau=50;
        else tau=tau_vector(i);
        end
        sf_responses_tillTau(j)=mean(SF_psths(j,baseline_duration+lat(j):baseline_duration+lat(j)+i_peak+tau));
        sf_responses_afterTau(j)=mean(SF_psths(j,baseline_duration+2000+lat(j):baseline_duration+2000+lat(j)+i_peak+tau));
        sf_responses_200(j)=mean(SF_psths(j,baseline_duration+50:baseline_duration+150));
        sf_responses_200(j)=mean(SF_psths(j,baseline_duration+1050:baseline_duration+1150));
    end
    
    initial_values=[max(sf_responses_tillTau) SFvalues_lin(find(sf_responses_tillTau==max(sf_responses_tillTau))) 0 1.5];
    lower_limits=[0.1,0.05,-1,0.1];
    upper_limits=[max(sf_responses_tillTau)+20,10,1,8];
    [fittedmodel,gof,output]=fit(SFvalues_lin',sf_responses_tillTau,priebe_1D,'maxIter',1000,'startpoint',initial_values,'lower',lower_limits,'upper',upper_limits);
    [fitp_SF_tillTau(i)]=fitp(gof,output);
    rsquare_SF_tillTau(i)=gof.rsquare;
    pref_SF_tillTau(i)=fittedmodel.ps;
    sigma_SF_tillTau(i)=fittedmodel.tuningWidth;
    eval(['sf_fittedmodels.cell' num2str(i) '.tillTau=fittedmodel;']);
    eval(['sf_gofs.cell' num2str(i) '.tillTau=gof;']);
    eval(['sf_outputs.cell' num2str(i) '.tillTau=output;']);
    
    initial_values=[max(sf_responses_afterTau) SFvalues_lin(find(sf_responses_afterTau==max(sf_responses_afterTau))) 0 1.5];
    lower_limits=[0.1,0.05,-1,0.1];
    upper_limits=[max(sf_responses_afterTau)+20,10,1,8];
    [fittedmodel,gof,output]=fit(SFvalues_lin',sf_responses_afterTau,priebe_1D,'maxIter',1000,'startpoint',initial_values,'lower',lower_limits,'upper',upper_limits);
    [fitp_SF_afterTau(i)]=fitp(gof,output);
    rsquare_SF_afterTau(i)=gof.rsquare;
    pref_SF_afterTau(i)=fittedmodel.ps;
    sigma_SF_afterTau(i)=fittedmodel.tuningWidth;
    eval(['sf_fittedmodels.cell' num2str(i) '.afterTau=fittedmodel;']);
    eval(['sf_gofs.cell' num2str(i) '.afterTau=gof;']);
    eval(['sf_outputs.cell' num2str(i) '.afterTau=output;']);
    
    initial_values=[max(sf_responses_200) SFvalues_lin(find(sf_responses_200==max(sf_responses_200))) 0 1.5];
    lower_limits=[0.1,0.05,-1,0.1];
    upper_limits=[max(sf_responses_200)+20,10,1,8];
    [fittedmodel,gof,output]=fit(SFvalues_lin',sf_responses_200,priebe_1D,'maxIter',1000,'startpoint',initial_values,'lower',lower_limits,'upper',upper_limits);
    [fitp_SF_200(i)]=fitp(gof,output);
    rsquare_SF_200(i)=gof.rsquare;
    pref_SF_200(i)=fittedmodel.ps;
    sigma_SF_200(i)=fittedmodel.tuningWidth;
    eval(['sf_fittedmodels.cell' num2str(i) '.till_150=fittedmodel;']);
    eval(['sf_gofs.cell' num2str(i) '.till_150=gof;']);
    eval(['sf_outputs.cell' num2str(i) '.till_150=output;']);
    
    initial_values=[max(sf_responses_800) SFvalues_lin(find(sf_responses_800==max(sf_responses_800))) 0 1.5];
    lower_limits=[0.1,0.05,-1,0.1];
    upper_limits=[max(sf_responses_800)+20,10,1,8];
    [fittedmodel,gof,output]=fit(SFvalues_lin',sf_responses_800,priebe_1D,'maxIter',1000,'startpoint',initial_values,'lower',lower_limits,'upper',upper_limits);
    [fitp_SF_800(i)]=fitp(gof,output);
    rsquare_SF_800(i)=gof.rsquare;
    pref_SF_800(i)=fittedmodel.ps;
    sigma_SF_800(i)=fittedmodel.tuningWidth;
    eval(['sf_fittedmodels.cell' num2str(i) '.after_1000=fittedmodel;']);
    eval(['sf_gofs.cell' num2str(i) '.after_1000=gof;']);
    eval(['sf_outputs.cell' num2str(i) '.after_1000=output;']);
    
    %%%%%%%% TF %%%%%%%%%%
    
    %get TF conditions and values for cell i
    eval(['TFs_i=TFs_i_atBestSF.complexCell' num2str(i) ';']);
    eval(['TFvalues_octv=TFvalues_octave.complexCell' num2str(i) ';']);

    TFvalues_lin=2.^TFvalues_octv;
    eval(['TFvalues_linScale.complexCell' num2str(i) '=TFvalues_lin;']);

    TF_psths=zeros(numel(TFs_i),numel(psth.condition1));
    for k=1:numel(TFs_i)
        TF_psths(k,:)=eval(['psth.condition' num2str(TFs_i(k))]);
        mean_baseline=mean(TF_psths(k,1:baseline_duration));
        TF_psths(k,:)=TF_psths(k,:)-mean_baseline;
    end
    
    %latencies for each condition, if nan, then lat=55 (mean latency)
    tf_column_i=find(TF_columnLabels_octv==TFvalues_octv(1));
    lat=lat_TFs_atBestSF(i,tf_column_i:tf_column_i+numel(TFvalues_octv)-1)';
    lat(find(isnan(lat)))=55;
    tf_responses_tillTau=zeros(size(lat));
    tf_responses_afterTau=zeros(size(lat));
    tf_responses_200=zeros(size(lat));
    tf_responses_800=zeros(size(lat));
    
    for j=1:numel(lat)
        if isnan(i_peakPsth(i))
            i_peak=nanmean(i_peakPsth(i));
        else i_peak=i_peakPsth(i);
        end
        if isnan(tau_vector(i)) || tau_vector(i)>1800
            tau=50;
        else tau=tau_vector(i);
        end
        tf_responses_tillTau(j)=mean(TF_psths(j,baseline_duration+lat(j):baseline_duration+lat(j)+i_peak+tau));
        tf_responses_afterTau(j)=mean(TF_psths(j,baseline_duration+2000+lat(j):baseline_duration+2000+lat(j)+i_peak+tau));
        tf_responses_200(j)=mean(TF_psths(j,baseline_duration+50:baseline_duration+150));
        tf_responses_200(j)=mean(TF_psths(j,baseline_duration+1050:baseline_duration+1150));
    end
    
    initial_values=[max(tf_responses_tillTau) TFvalues_lin(find(tf_responses_tillTau==max(tf_responses_tillTau))) 0 1.5];
    lower_limits=[0.1,0.05,-1,0.1];
    upper_limits=[max(tf_responses_tillTau)+20,10,1,8];
    [fittedmodel,gof,output]=fit(TFvalues_lin',tf_responses_tillTau,priebe_1D,'maxIter',1000,'startpoint',initial_values,'lower',lower_limits,'upper',upper_limits);
    [fitp_TF_tillTau(i)]=fitp(gof,output);
    rsquare_TF_tillTau(i)=gof.rsquare;
    pref_TF_tillTau(i)=fittedmodel.ps;
    sigma_TF_tillTau(i)=fittedmodel.tuningWidth;
    eval(['tf_fittedmodels.cell' num2str(i) '.tillTau=fittedmodel;']);
    eval(['tf_gofs.cell' num2str(i) '.tillTau=gof;']);
    eval(['TF_outputs.cell' num2str(i) '.tillTau=output;']);
    
    initial_values=[max(tf_responses_afterTau) TFvalues_lin(find(tf_responses_afterTau==max(tf_responses_afterTau))) 0 1.5];
    lower_limits=[0.1,0.05,-1,0.1];
    upper_limits=[max(tf_responses_afterTau)+20,10,1,8];
    [fittedmodel,gof,output]=fit(TFvalues_lin',tf_responses_afterTau,priebe_1D,'maxIter',1000,'startpoint',initial_values,'lower',lower_limits,'upper',upper_limits);
    [fitp_TF_afterTau(i)]=fitp(gof,output);
    rsquare_TF_afterTau(i)=gof.rsquare;
    pref_TF_afterTau(i)=fittedmodel.ps;
    sigma_TF_afterTau(i)=fittedmodel.tuningWidth;
    eval(['tf_fittedmodels.cell' num2str(i) '.afterTau=fittedmodel;']);
    eval(['tf_gofs.cell' num2str(i) '.afterTau=gof;']);
    eval(['tf_outputs.cell' num2str(i) '.afterTau=output;']);
    
    initial_values=[max(tf_responses_200) TFvalues_lin(find(tf_responses_200==max(tf_responses_200))) 0 1.5];
    lower_limits=[0.1,0.05,-1,0.1];
    upper_limits=[max(tf_responses_200)+20,10,1,8];
    [fittedmodel,gof,output]=fit(TFvalues_lin',tf_responses_200,priebe_1D,'maxIter',1000,'startpoint',initial_values,'lower',lower_limits,'upper',upper_limits);
    [fitp_TF_200(i)]=fitp(gof,output);
    rsquare_TF_200(i)=gof.rsquare;
    pref_TF_200(i)=fittedmodel.ps;
    sigma_TF_200(i)=fittedmodel.tuningWidth;
    eval(['tf_fittedmodels.cell' num2str(i) '.till_150=fittedmodel;']);
    eval(['tf_gofs.cell' num2str(i) '.till_150=gof;']);
    eval(['tf_outputs.cell' num2str(i) '.till_150=output;']);
    
    initial_values=[max(tf_responses_800) TFvalues_lin(find(tf_responses_800==max(tf_responses_800))) 0 1.5];
    lower_limits=[0.1,0.05,-1,0.1];
    upper_limits=[max(tf_responses_800)+20,10,1,8];
    [fittedmodel,gof,output]=fit(TFvalues_lin',tf_responses_800,priebe_1D,'maxIter',1000,'startpoint',initial_values,'lower',lower_limits,'upper',upper_limits);
    [fitp_TF_800(i)]=fitp(gof,output);
    rsquare_TF_800(i)=gof.rsquare;
    pref_TF_800(i)=fittedmodel.ps;
    sigma_TF_800(i)=fittedmodel.tuningWidth;
    eval(['tf_fittedmodels.cell' num2str(i) '.after_1000=fittedmodel;']);
    eval(['tf_gofs.cell' num2str(i) '.after_1000=gof;']);
    eval(['tf_outputs.cell' num2str(i) '.after_1000=output;']);

end
