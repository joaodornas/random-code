non_empty_cols_lsfs=nan(1,analysis_period);
non_empty_cols_bw=nan(1,analysis_period);
non_empty_cols_prefTF=nan(1,analysis_period);
non_empty_cols_sigma=nan(1,analysis_period);
non_empty_cols_hsfs=nan(1,analysis_period);
non_empty_cols_amp=nan(1,analysis_period);

for k=1:analysis_period
    non_empty_cols_lsfs(k)=numel(find(~isnan(LSFSs_dyn(:,k))));
    non_empty_cols_bw(k)=numel(find(~isnan(BWs_dyn(:,k))));
    non_empty_cols_prefTF(k)=numel(find(~isnan(pref_TFs_dyn(:,k))));
    non_empty_cols_sigma(k)=numel(find(~isnan(sigmas_dyn(:,k))));
    non_empty_cols_hsfs(k)=numel(find(~isnan(HSFSs_dyn(:,k))));
    non_empty_cols_amp(k)=numel(find(~isnan(response_amp_dyn(:,k))));
end

mean_lsfs_dyn=nanmean(LSFSs_dyn);
sem_lsfs_dyn=nanstd(LSFSs_dyn)./sqrt(non_empty_cols_lsfs);
mean_bw_dyn=nanmean(BWs_dyn);
sem_bw_dyn=nanstd(BWs_dyn)./sqrt(non_empty_cols_bw);
mean_prefTF_dyn=nanmean(pref_TFs_dyn);
sem_preTF_dyn=nanstd(pref_TFs_dyn)./sqrt(non_empty_cols_prefTF);
mean_sigma_dyn=nanmean(sigmas_dyn);
sem_sigma_dyn=nanstd(sigmas_dyn)./sqrt(non_empty_cols_sigma);
mean_hsfs_dyn=nanmean(HSFSs_dyn);
sem_hsfs_dyn=nanstd(HSFSs_dyn)./sqrt(non_empty_cols_hsfs);
mean_amp_dyn=nanmean(response_amp_dyn);
sem_amp_dyn=nanstd(response_amp_dyn)./sqrt(non_empty_cols_amp);

delta_tf=nanmean(pref_TFs_dyn(:,analysis_period-9:analysis_period),2)-nanmean(pref_TFs_dyn(:,1:10),2);
mean_delta_tf=nanmean(delta_tf)
sem_delta_tf=nanstd(delta_tf)/sqrt(numel(find(~isnan(delta_tf))));

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

delta_tf_octv=log2(nanmean(pref_TFs_dyn(:,analysis_period-9:analysis_period),2))-log2(nanmean(pref_TFs_dyn(:,1:10),2));
mean_delta_tf_octv=nanmean(delta_tf_octv);
sem_delta_tf_octv=nanstd(delta_tf_octv)/sqrt(numel(find(~isnan(delta_tf_octv))));


h=lillietest(delta_tf_octv);
if h==1
    p_signrank_delta_tf_octv=signrank(delta_tf_octv)
elseif h==0
    [a,p_ttest_delta_tf_octv]=ttest(delta_tf_octv)
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

[h,p_runs_sigma]=runstest(mean_sigma_dyn)
[h,p_runs_prefSF]=runstest(mean_prefTF_dyn)
[h,p_runs_lsfs]=runstest(mean_lsfs_dyn)
[h,p_runs_hsfs]=runstest(mean_hsfs_dyn)
[h,p_runs_amp]=runstest(mean_amp_dyn)

clear h

slopes_tf=nan(92,1);
slopes_sigma=nan(92,1);
slopes_lsfs=nan(92,1);
slopes_hsfs=nan(92,1);
slopes_amp=nan(92,1);
x_curve=1:analysis_period;

for i=[1:4,6:12,14:92]
    curve=pref_TFs_dyn(i,:);
    curve=curve(find(~isnan(curve)));
    x=x_curve(find(~isnan(curve)));
    fitResult=fit(x',curve','poly1');
    slopes_tf(i)=fitResult.p1;

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
    
    curve=response_amp_dyn(i,:);
    curve=curve(find(~isnan(curve)));
    x=x_curve(find(~isnan(curve)));
    fitResult=fit(x',curve','poly1');
    slopes_amp(i)=fitResult.p1;
    
    curve=HSFSs_dyn(i,:);
    curve=curve(find(~isnan(curve)));
    x=x_curve(find(~isnan(curve)));
    fitResult=fit(x',curve','poly1');
    slopes_hsfs(i)=fitResult.p1;
    
    clear curve fitResult x
end

clear x_curve i k

mean_slopes_amp=nanmean(slopes_amp)
sem_slopes_amp=nanstd(slopes_amp)/sqrt(numel(find(~isnan(slopes_amp))));

mean_slopes_tf=nanmean(slopes_tf)
sem_slopes_tf=nanstd(slopes_tf)/sqrt(numel(find(~isnan(slopes_tf))));

mean_slopes_sigma=nanmean(slopes_sigma)
sem_slopes_sigma=nanstd(slopes_sigma)/sqrt(numel(find(~isnan(slopes_sigma))));

mean_slopes_lsfs=nanmean(slopes_lsfs)
sem_slopes_lsfs=nanstd(slopes_lsfs)/sqrt(numel(find(~isnan(slopes_lsfs))));

mean_slopes_hsfs=nanmean(slopes_hsfs)
sem_slopes_hsfs=nanstd(slopes_hsfs)/sqrt(numel(find(~isnan(slopes_hsfs))));

p_signrank_slopes_tf=signrank(slopes_tf)
p_signrank_slopes_sigma=signrank(slopes_sigma)
p_signrank_slopes_lsfs=signrank(slopes_lsfs)
p_signrank_slopes_amp=signrank(slopes_amp)
p_signrank_slopes_hsfs=signrank(slopes_hsfs)

diff_amp=nan(ncells,analysis_period-1); sign_diff_amp=nan(size(diff_amp)); 
diff_lsfs=nan(ncells,analysis_period-1); sign_diff_lsfs=nan(size(diff_lsfs)); 
diff_hsfs=nan(ncells,analysis_period-1); sign_diff_hsfs=nan(size(diff_hsfs)); 
diff_tf=nan(ncells,analysis_period-1); sign_diff_tf=nan(size(diff_tf)); 
diff_sigma=nan(ncells,analysis_period-1); sign_diff_sigma=nan(size(diff_sigma)); 

for i=[1:4,6:12,14:92]
% diff_amp(i,:)=diff(response_amp_dyn(i,:)); sign_diff_amp(i,:)=sign(diff_amp(i,:)); 
% sign_index_amp(i)=(numel(find(sign_diff_amp(i,:)==1))*abs(sum(diff_amp(i,find(sign_diff_amp(i,:)==1))))-numel(find(sign_diff_amp(i,:)==-1))*abs(sum(diff_amp(i,find(sign_diff_amp(i,:)==-1)))))/(numel(find(sign_diff_amp(i,:)==1))*abs(sum(diff_amp(i,find(sign_diff_amp(i,:)==1))))+numel(find(sign_diff_amp(i,:)==-1))*abs(sum(diff_amp(i,find(sign_diff_amp(i,:)==-1)))));
% 
% diff_lsfs(i,:)=diff(LSFSs_dyn(i,:)); sign_diff_lsfs(i,:)=sign(diff_lsfs(i,:)); 
% sign_index_lsfs(i)=(abs(sum(diff_lsfs(i,find(sign_diff_lsfs(i,:)==1))))-abs(sum(diff_lsfs(i,find(sign_diff_lsfs(i,:)==-1)))))/(abs(sum(diff_lsfs(i,find(sign_diff_lsfs(i,:)==1))))+abs(sum(diff_lsfs(i,find(sign_diff_lsfs(i,:)==-1)))));

% diff_hsfs(i,:)=diff(HSFSs_dyn(i,:)); sign_diff_hsfs(i,:)=sign(diff_hsfs(i,:)); 
% sign_index_hsfs(i)=(numel(find(sign_diff_hsfs(i,:)==1))*abs(sum(diff_hsfs(i,find(sign_diff_hsfs(i,:)==1))))-numel(find(sign_diff_hsfs(i,:)==-1))*abs(sum(diff_hsfs(i,find(sign_diff_hsfs(i,:)==-1)))))/(numel(find(sign_diff_hsfs(i,:)==1))*abs(sum(diff_hsfs(i,find(sign_diff_hsfs(i,:)==1))))+numel(find(sign_diff_hsfs(i,:)==-1))*abs(sum(diff_hsfs(i,find(sign_diff_hsfs(i,:)==-1)))));
% 
% diff_tf(i,:)=diff(pref_TFs_dyn(i,:)); sign_diff_tf(i,:)=sign(diff_tf(i,:)); 
% sign_index_tf(i)=(numel(find(sign_diff_tf(i,:)==1))*abs(sum(diff_tf(i,find(sign_diff_tf(i,:)==1))))-numel(find(sign_diff_tf(i,:)==-1))*abs(sum(diff_tf(i,find(sign_diff_tf(i,:)==-1)))))/(numel(find(sign_diff_tf(i,:)==1))*abs(sum(diff_tf(i,find(sign_diff_tf(i,:)==1))))+numel(find(sign_diff_tf(i,:)==-1))*abs(sum(diff_tf(i,find(sign_diff_tf(i,:)==-1)))));
% 
% diff_sigma(i,:)=diff(sigmas_dyn(i,:)); sign_diff_sigma(i,:)=sign(diff_sigma(i,:)); 
% sign_index_sigma(i)=(numel(find(sign_diff_sigma(i,:)==1))*abs(sum(diff_sigma(i,find(sign_diff_sigma(i,:)==1))))-numel(find(sign_diff_sigma(i,:)==-1))*abs(sum(diff_sigma(i,find(sign_diff_sigma(i,:)==-1)))))/(numel(find(sign_diff_sigma(i,:)==1))*abs(sum(diff_sigma(i,find(sign_diff_sigma(i,:)==1))))+numel(find(sign_diff_sigma(i,:)==-1))*abs(sum(diff_sigma(i,find(sign_diff_sigma(i,:)==-1)))));
end

sign_index_amp=nan(ncells,1);
sign_index_lsfs=nan(ncells,1);
sign_index_hsfs=nan(ncells,1);
sign_index_tf=nan(ncells,1);
sign_index_sigma=nan(ncells,1);

for i=[1:4,6:12,14:92]

    a=find(~isnan(response_amp_dyn(i,:)),1,'first');
    if a < 10
        larger=find(response_amp_dyn(i,:)>response_amp_dyn(i,a));
        smaller=find(response_amp_dyn(i,:)<response_amp_dyn(i,a));
        sign_index_amp(i)=(numel(larger)-numel(smaller))/(numel(larger)+numel(smaller));
    end

    a=find(~isnan(LSFSs_dyn(i,:)),1,'first');
    if a < 10
        larger=find(LSFSs_dyn(i,:)>LSFSs_dyn(i,a));
        smaller=find(LSFSs_dyn(i,:)<LSFSs_dyn(i,a));
        sign_index_lsfs(i)=(numel(larger)-numel(smaller))/(numel(larger)+numel(smaller));
    end

    a=find(~isnan(HSFSs_dyn(i,:)),1,'first');
    if a < 10
        larger=find(HSFSs_dyn(i,:)>HSFSs_dyn(i,a));
        smaller=find(HSFSs_dyn(i,:)<HSFSs_dyn(i,a));
        sign_index_hsfs(i)=(numel(larger)-numel(smaller))/(numel(larger)+numel(smaller));
    end

    a=find(~isnan(sigmas_dyn(i,:)),1,'first');
    if a < 10
        larger=find(sigmas_dyn(i,:)>sigmas_dyn(i,a));
        smaller=find(sigmas_dyn(i,:)<sigmas_dyn(i,a));
        sign_index_sigma(i)=(numel(larger)-numel(smaller))/(numel(larger)+numel(smaller));
    end

    a=find(~isnan(pref_TFs_dyn(i,:)),1,'first');
    if a < 10
        larger=find(pref_TFs_dyn(i,:)>pref_TFs_dyn(i,a));
        smaller=find(pref_TFs_dyn(i,:)<pref_TFs_dyn(i,a));
        sign_index_tf(i)=(numel(larger)-numel(smaller))/(numel(larger)+numel(smaller));
    end
end

clear larger smaller i a

mean_sign_index_amp=nanmean(sign_index_amp) 
sem_sign_index_amp=nanstd(sign_index_amp)/sqrt(numel(find(~isnan(sign_index_amp))));
mean_sign_index_lsfs=nanmean(sign_index_lsfs) 
sem_sign_index_lsfs=nanstd(sign_index_lsfs)/sqrt(numel(find(~isnan(sign_index_lsfs))));
mean_sign_index_hsfs=nanmean(sign_index_hsfs) 
sem_sign_index_hsfs=nanstd(sign_index_hsfs)/sqrt(numel(find(~isnan(sign_index_hsfs))));
mean_sign_index_sigma=nanmean(sign_index_sigma)
sem_sign_index_sigma=nanstd(sign_index_sigma)/sqrt(numel(find(~isnan(sign_index_sigma))));
mean_sign_index_tf=nanmean(sign_index_tf)
sem_sign_index_tf=nanstd(sign_index_tf)/sqrt(numel(find(~isnan(sign_index_tf))));

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

h=lillietest(sign_index_tf);
if h==1
    p_signrank_sign_index_tf=signrank(sign_index_tf)
elseif h==0
    [a,p_ttest_sign_index_tf]=ttest(sign_index_tf)
    clear a h
end

[y_hist_signInd_hsfs,x_hist_signInd_hsfs]=makeHistogram(sign_index_hsfs,0.2);
[y_hist_signInd_lsfs,x_hist_signInd_lsfs]=makeHistogram(sign_index_lsfs,0.2);
[y_hist_signInd_sigma,x_hist_signInd_sigma]=makeHistogram(sign_index_sigma,0.2);
[y_hist_signInd_prefTF,x_hist_signInd_prefTF]=makeHistogram(sign_index_tf,0.2);

% [y_hist_delta_hsfs,x_hist_delta_hsfs]=makeHistogram(delta_hsfs,0.25);
% [y_hist_delta_lsfs,x_hist_delta_lsfs]=makeHistogram(delta_lsfs,0.25);
% [y_hist_delta_sigma,x_hist_delta_sigma]=makeHistogram(delta_sigma,1);
% [y_hist_delta_prefTF_octv,x_hist_delta_prefTF_octv]=makeHistogram(delta_tf_octv,1);

save([savepath 'TF_tuningDynamics_analysis_windowSize' num2str(window_size) '_anaPeriod' num2str(analysis_period) '.mat']);