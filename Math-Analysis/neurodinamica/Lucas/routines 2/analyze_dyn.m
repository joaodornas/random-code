non_empty_cols_lsfs=nan(1,analysis_period);
non_empty_cols_hsfs=nan(1,analysis_period);
non_empty_cols_bw=nan(1,analysis_period);
non_empty_cols_prefSF=nan(1,analysis_period);
non_empty_cols_sigma=nan(1,analysis_period);
non_empty_cols_amp=nan(1,analysis_period);

for k=1:analysis_period
    non_empty_cols_lsfs(k)=numel(find(~isnan(LSFSs_dyn(:,k))));
    non_empty_cols_hsfs(k)=numel(find(~isnan(HSFSs_dyn(:,k))));
    non_empty_cols_bw(k)=numel(find(~isnan(BWs_dyn(:,k))));
    non_empty_cols_prefSF(k)=numel(find(~isnan(pref_SFs_dyn(:,k))));
    non_empty_cols_sigma(k)=numel(find(~isnan(sigmas_dyn(:,k))));
    non_empty_cols_amp(k)=numel(find(~isnan(response_amp_dyn(:,k))));
end

mean_lsfs_dyn=nanmean(LSFSs_dyn);
sem_lsfs_dyn=nanstd(LSFSs_dyn)./sqrt(non_empty_cols_lsfs);
mean_hsfs_dyn=nanmean(HSFSs_dyn);
sem_hsfs_dyn=nanstd(HSFSs_dyn)./sqrt(non_empty_cols_hsfs);
mean_bw_dyn=nanmean(BWs_dyn);
sem_bw_dyn=nanstd(BWs_dyn)./sqrt(non_empty_cols_bw);
mean_prefSF_dyn=nanmean(pref_SFs_dyn);
sem_preSF_dyn=nanstd(pref_SFs_dyn)./sqrt(non_empty_cols_prefSF);
mean_sigma_dyn=nanmean(sigmas_dyn);
sem_sigma_dyn=nanstd(sigmas_dyn)./sqrt(non_empty_cols_sigma);
mean_amp_dyn=nanmean(response_amp_dyn);
sem_amp_dyn=nanstd(response_amp_dyn)./sqrt(non_empty_cols_amp);

delta_sf=nanmean(pref_SFs_dyn(:,analysis_period-9:analysis_period),2)-nanmean(pref_SFs_dyn(:,1:10),2);
mean_delta_sf=nanmean(delta_sf)
sem_delta_sf=nanstd(delta_sf)/sqrt(numel(find(~isnan(delta_sf))));

delta_sigma=nanmean(sigmas_dyn(:,analysis_period-9:analysis_period),2)-nanmean(sigmas_dyn(:,1:10),2);
mean_delta_sigma=nanmean(delta_sigma)
sem_delta_sigma=nanstd(delta_sigma)/sqrt(numel(find(~isnan(delta_sigma))));

delta_lsfs=nanmean(LSFSs_dyn(:,analysis_period-9:analysis_period),2)-nanmean(LSFSs_dyn(:,1:10),2);
mean_delta_lsfs=nanmean(delta_lsfs)
sem_delta_lsfs=nanstd(delta_lsfs)/sqrt(numel(find(~isnan(delta_lsfs))));

delta_hsfs=nanmean(HSFSs_dyn(:,analysis_period-9:analysis_period),2)-nanmean(HSFSs_dyn(:,1:10),2);
mean_delta_hsfs=nanmean(delta_hsfs)
sem_delta_hsfs=nanstd(delta_hsfs)/sqrt(numel(find(~isnan(delta_hsfs))));

delta_amp=nanmean(response_amp_dyn(:,analysis_period-9:analysis_period),2)-nanmean(response_amp_dyn(:,1:10),2);
mean_delta_amp=nanmean(delta_amp)
sem_delta_amp=nanstd(delta_amp)/sqrt(numel(find(~isnan(delta_amp))));

delta_sf_octv=log2(nanmean(pref_SFs_dyn(:,analysis_period-9:analysis_period),2))-log2(nanmean(pref_SFs_dyn(:,1:10),2));
mean_delta_sf_octv=nanmean(delta_sf_octv);
sem_delta_sf_octv=nanstd(delta_sf_octv)/sqrt(numel(find(~isnan(delta_sf_octv))));

[h,p_runs_sigma]=runstest(mean_sigma_dyn)
[h,p_runs_prefSF]=runstest(mean_prefSF_dyn)
[h,p_runs_lsfs]=runstest(mean_lsfs_dyn)
[h,p_runs_hsfs]=runstest(mean_hsfs_dyn)
[h,p_runs_amp]=runstest(mean_amp_dyn)

clear h

p_signrank_delta_sf_octv=signrank(delta_sf_octv)
p_signrank_delta_sigma=signrank(delta_sigma)
p_signrank_delta_lsfs=signrank(delta_lsfs)
p_signrank_delta_hsfs=signrank(delta_hsfs)
p_signrank_delta_amp=signrank(delta_amp)

pref_SFs_dyn_octv=log2(pref_SFs_dyn);
mean_prefSF_dyn_octv=nanmean(pref_SFs_dyn_octv);
sem_prefSF_dyn_octv=nanstd(pref_SFs_dyn_octv)./sqrt(non_empty_cols_prefSF);

slopes_sf_octv=nan(92,1);
slopes_sigma=nan(92,1);
slopes_lsfs=nan(92,1);
slopes_hsfs=nan(92,1);
slopes_amp=nan(92,1);
x_curve=1:analysis_period;

for i=[1:4,6:12,14:92]
    
    try
        curve=pref_SFs_dyn_octv(i,:);
        curve=curve(find(~isnan(curve)));
        x=x_curve(find(~isnan(curve)));
        fitResult=fit(x',curve','poly1');
        slopes_sf_octv(i)=fitResult.p1;
    catch
        slopes_sf_octv(i)=nan;
    end

    curve=sigmas_dyn(i,:);
    curve=curve(find(~isnan(curve)));
    x=x_curve(find(~isnan(curve)));
    fitResult=fit(x',curve','poly1');
    slopes_sigma(i)=fitResult.p1;

    curve=LSFSs_dyn(i,:);
    curve=curve(find(~isnan(curve)));
    x=x_curve(find(~isnan(curve)));
    fitResult=fit(x',curve','poly1');
    slopes_lsfs(i)=fitResult.p1;
    
    curve=HSFSs_dyn(i,:);
    curve=curve(find(~isnan(curve)));
    x=x_curve(find(~isnan(curve)));
    fitResult=fit(x',curve','poly1');
    slopes_hsfs(i)=fitResult.p1;
    
    curve=response_amp_dyn(i,:);
    curve=curve(find(~isnan(curve)));
    x=x_curve(find(~isnan(curve)));
    fitResult=fit(x',curve','poly1');
    slopes_amp(i)=fitResult.p1;

    clear curve fitResult x
end

clear x_curve

mean_slopes_sf_octv=nanmean(slopes_sf_octv)
sem_slopes_sf=nanstd(slopes_sf_octv)/sqrt(numel(find(~isnan(slopes_sf_octv))));

mean_slopes_sigma=nanmean(slopes_sigma)
sem_slopes_sigma=nanstd(slopes_sigma)/sqrt(numel(find(~isnan(slopes_sigma))));

mean_slopes_lsfs=nanmean(slopes_lsfs)
sem_slopes_lsfs=nanstd(slopes_lsfs)/sqrt(numel(find(~isnan(slopes_lsfs))));

mean_slopes_hsfs=nanmean(slopes_hsfs)
sem_slopes_hsfs=nanstd(slopes_hsfs)/sqrt(numel(find(~isnan(slopes_hsfs))));

mean_slopes_amp=nanmean(slopes_amp)
sem_slopes_amp=nanstd(slopes_amp)/sqrt(numel(find(~isnan(slopes_amp))));

p_signrank_slopes_sf_octv=signrank(slopes_sf_octv)
p_signrank_slopes_sigma=signrank(slopes_sigma)
p_signrank_slopes_lsfs=signrank(slopes_lsfs)
p_signrank_slopes_hsfs=signrank(slopes_hsfs)
p_signrank_slopes_amp=signrank(slopes_amp)


h=lillietest(delta_sf_octv);
if h==1
    p_signrank_delta_sf_octv=signrank(delta_sf_octv)
elseif h==0
    [a,p_ttest_delta_sf_octv]=ttest(delta_sf_octv)
    clear a
end

h=lillietest(delta_hsfs);
if h==1
    p_signrank_delta_hsfs=signrank(delta_hsfs)
elseif h==0
    [a,p_ttest_delta_hsfs]=ttest(delta_hsfs)
    clear a p_signrank_delta_hsfs
end

h=lillietest(delta_lsfs);
if h==1
    p_signrank_delta_lsfs=signrank(delta_lsfs)
elseif h==0
    [a,p_ttest_delta_lsfs]=ttest(delta_lsfs)
    clear a p_signrank_delta_lsfs
end

h=lillietest(delta_sigma);
if h==1
    p_signrank_delta_sigma=signrank(delta_sigma)
elseif h==0
    [a,p_ttest_delta_sigma]=ttest(delta_sigma)
    clear a p_signrank_delta_sigma
end

clear h

prefSF_octv_plusSEM=mean_prefSF_dyn_octv+sem_prefSF_dyn_octv;
prefSF_octv_minusSEM=mean_prefSF_dyn_octv-sem_prefSF_dyn_octv;
sigma_plusSEM=mean_sigma_dyn+sem_sigma_dyn;
sigma_minusSEM=mean_sigma_dyn-sem_sigma_dyn;
lsfs_plusSEM=mean_lsfs_dyn+sem_lsfs_dyn;
lsfs_minusSEM=mean_lsfs_dyn-sem_lsfs_dyn;
hsfs_plusSEM=mean_hsfs_dyn+sem_hsfs_dyn;
hsfs_minusSEM=mean_hsfs_dyn-sem_hsfs_dyn;

[y_hist_delta_hsfs,x_hist_delta_hsfs]=makeHistogram(delta_hsfs,0.25);
[y_hist_delta_lsfs,x_hist_delta_lsfs]=makeHistogram(delta_lsfs,0.25);
[y_hist_delta_sigma,x_hist_delta_sigma]=makeHistogram(delta_sigma,1);
[y_hist_delta_prefSF_octv,x_hist_delta_prefSF_octv]=makeHistogram(delta_sf_octv,1);

sign_index_amp=nan(ncells,1);
sign_index_lsfs50=nan(ncells,1);
sign_index_hsfs50=nan(ncells,1);
sign_index_sf50=nan(ncells,1);
sign_index_sigma50=nan(ncells,1);

for i=[1:4,6:12,14:92]
%     larger=find(response_amp_dyn(i,:)>response_amp_dyn(i,1));
%     smaller=find(response_amp_dyn(i,:)<response_amp_dyn(i,1));
%     sign_index_amp(i)=(numel(larger)-numel(smaller))/(numel(larger)+numel(smaller));
%     
    larger=find(LSFSs_dyn(i,2:51)>LSFSs_dyn(i,1));
    smaller=find(LSFSs_dyn(i,2:51)<LSFSs_dyn(i,1));
    sign_index_lsfs50(i)=(numel(larger)-numel(smaller))/(numel(larger)+numel(smaller));
    
    larger=find(HSFSs_dyn(i,2:51)>HSFSs_dyn(i,1));
    smaller=find(HSFSs_dyn(i,2:51)<HSFSs_dyn(i,1));
    sign_index_hsfs50(i)=(numel(larger)-numel(smaller))/(numel(larger)+numel(smaller));
    
    larger=find(sigmas_dyn(i,2:51)>sigmas_dyn(i,1));
    smaller=find(sigmas_dyn(i,:)<sigmas_dyn(i,1));
    sign_index_sigma50(i)=(numel(larger)-numel(smaller))/(numel(larger)+numel(smaller));
    
    larger=find(pref_SFs_dyn(i,2:51)>pref_SFs_dyn(i,1));
    smaller=find(pref_SFs_dyn(i,2:51)<pref_SFs_dyn(i,1));
    sign_index_sf59(i)=(numel(larger)-numel(smaller))/(numel(larger)+numel(smaller));
end

clear larger smaller i

mean_sign_index_amp=nanmean(sign_index_amp) 
sem_sign_index_amp=nanstd(sign_index_amp)/sqrt(numel(find(~isnan(sign_index_amp))));
mean_sign_index_lsfs=nanmean(sign_index_lsfs) 
sem_sign_index_lsfs=nanstd(sign_index_lsfs)/sqrt(numel(find(~isnan(sign_index_lsfs))));
mean_sign_index_hsfs=nanmean(sign_index_hsfs) 
sem_sign_index_hsfs=nanstd(sign_index_hsfs)/sqrt(numel(find(~isnan(sign_index_hsfs))));
mean_sign_index_sigma=nanmean(sign_index_sigma)
sem_sign_index_sigma=nanstd(sign_index_sigma)/sqrt(numel(find(~isnan(sign_index_sigma))));
mean_sign_index_sf=nanmean(sign_index_sf)
sem_sign_index_sf=nanstd(sign_index_sf)/sqrt(numel(find(~isnan(sign_index_sf))));

h=lillietest(sign_index_amp);
if h==1
    p_signrank_sign_index_amp=signrank(sign_index_amp)
elseif h==0
    [a,p_ttest_sign_index_amp]=ttest(sign_index_amp)
    clear a
end

h=lillietest(sign_index_lsfs);
if h==1
    p_signrank_sign_index_lsfs=signrank(sign_index_lsfs)
elseif h==0
    [a,p_ttest_sign_index_lsfs]=ttest(sign_index_lsfs)
    clear a 
end

h=lillietest(sign_index_hsfs);
if h==1
    p_signrank_sign_index_hsfs=signrank(sign_index_hsfs)
elseif h==0
    [a,p_ttest_sign_index_hsfs]=ttest(sign_index_hsfs)
    clear a
end

h=lillietest(sign_index_sigma);
if h==1
    p_signrank_sign_index_sigma=signrank(sign_index_sigma)
elseif h==0
    [a,p_ttest_sign_index_sigma]=ttest(sign_index_sigma)
    clear a 
end

h=lillietest(sign_index_sf);
if h==1
    p_signrank_sign_index_sf=signrank(sign_index_sf)
elseif h==0
    [a,p_ttest_sign_index_sf]=ttest(sign_index_sf)
    clear a h
end

[y_hist_signInd_hsfs,x_hist_signInd_hsfs]=makeHistogram(sign_index_hsfs,0.2);
[y_hist_signInd_lsfs,x_hist_signInd_lsfs]=makeHistogram(sign_index_lsfs,0.2);
[y_hist_signInd_sigma,x_hist_signInd_sigma]=makeHistogram(sign_index_sigma,0.2);
[y_hist_signInd_prefSF,x_hist_signInd_prefSF]=makeHistogram(sign_index_sf,0.2);

save([savepath 'SF_tuningDynamics_analysis_windowSize' num2str(window_size) '_anaPeriod' num2str(analysis_period) '.mat']);

for i=20:20:980
runningDelta_sf(:,i/20)=log2(nanmean(pref_SFs_dyn(:,i:i+10),2))-log2(nanmean(pref_SFs_dyn(:,1:10),2));
runningDelta_sigma(:,i/20)=nanmean(sigmas_dyn(:,i:i+10),2)-nanmean(sigmas_dyn(:,1:10),2);
runningDelta_lsfs(:,i/20)=nanmean(LSFSs_dyn(:,i:i+10),2)-nanmean(LSFSs_dyn(:,1:10),2);
runningDelta_hsfs(:,i/20)=nanmean(HSFSs_dyn(:,i:i+10),2)-nanmean(HSFSs_dyn(:,1:10),2);
end