load('TF_tuningDynamics_analysis_windowSize10_anaPeriod1000','sigmas_dyn','LSFSs_dyn','HSFSs_dyn','pref_TFs_dyn_octv', ...
    'delta_tf','delta_sigma','delta_lsfs','delta_hsfs','sign_index_tf','sign_index_sigma','sign_index_lsfs','sign_index_hsfs');
load('bootstrap_sf_dyn_analysis','boot_i')

% % bootstrapping indexes, 20 cells, 1000 iterations
% boot_i=ceil(92.*rand(20,1000)); 
% 
% % exclude 5's and 13's (all nan)
% a=nan;
% for i=1:1000
%     while ~isempty(a)
%         a=find(boot_i(:,i)==5 | boot_i(:,i)==13);
%         boot_i(a,i)=ceil(92.*rand(numel(a),1));
%     end
% end

%initialize matrices
pref_TF_boot=nan(1000,1000);
sigma_boot=nan(1000,1000);
lsfs_boot=nan(1000,1000);
hsfs_boot=nan(1000,1000);

%bootstrap, saving mean values per iteration
for i=1:1000
    pref_TF_boot(i,:)=nanmean(pref_TFs_dyn_octv(boot_i(:,i),:));
    sigma_boot(i,:)=nanmean(sigmas_dyn(boot_i(:,i),:));
    lsfs_boot(i,:)=nanmean(LSFSs_dyn(boot_i(:,i),:));
    hsfs_boot(i,:)=nanmean(HSFSs_dyn(boot_i(:,i),:));
end

mean_pref_TF_boot=nanmean(pref_TF_boot);
mean_sigma_boot=nanmean(sigma_boot);
mean_lsfs_boot=nanmean(lsfs_boot);
mean_hsfs_boot=nanmean(hsfs_boot);

mean_pref_TF=nanmean(pref_TFs_dyn_octv);
mean_sigma=nanmean(sigmas_dyn);
mean_lsfs=nanmean(LSFSs_dyn);
mean_hsfs=nanmean(HSFSs_dyn);

sign_index_lsfs_boot=nan(1000,1);
sign_index_hsfs_boot=nan(1000,1);
sign_index_tf_boot=nan(1000,1);
sign_index_sigma_boot=nan(1000,1);

for i=1:1000
    
    larger=find(lsfs_boot(i,:)>lsfs_boot(i,1));
    smaller=find(lsfs_boot(i,:)<lsfs_boot(i,1));
    sign_index_lsfs_boot(i)=(numel(larger)-numel(smaller))/(numel(larger)+numel(smaller));
    
    larger=find(hsfs_boot(i,:)>hsfs_boot(i,1));
    smaller=find(hsfs_boot(i,:)<hsfs_boot(i,1));
    sign_index_hsfs_boot(i)=(numel(larger)-numel(smaller))/(numel(larger)+numel(smaller));
    
    larger=find(sigma_boot(i,:)>sigma_boot(i,1));
    smaller=find(sigma_boot(i,:)<sigma_boot(i,1));
    sign_index_sigma_boot(i)=(numel(larger)-numel(smaller))/(numel(larger)+numel(smaller));
    
    larger=find(pref_TF_boot(i,1:100)>pref_TF_boot(i,1));
    smaller=find(pref_TF_boot(i,1:100)<pref_TF_boot(i,1));
    sign_index_tf_boot(i)=(numel(larger)-numel(smaller))/(numel(larger)+numel(smaller));
end

delta_tf_boot=nanmean(pref_TF_boot(:,901:1000),2)-nanmean(pref_TF_boot(:,1:10),2);
delta_sigma_boot=nanmean(sigma_boot(:,901:1000),2)-nanmean(sigma_boot(:,1:10),2);
delta_lsfs_boot=nanmean(lsfs_boot(:,901:1000),2)-nanmean(lsfs_boot(:,1:10),2);
delta_hsfs_boot=nanmean(hsfs_boot(:,901:1000),2)-nanmean(hsfs_boot(:,1:10),2);

p_deltaTF=ranksum(delta_tf(find(~isnan(delta_tf))),delta_tf_boot(find(~isnan(delta_tf_boot))))
p_deltaSigma=ranksum(delta_sigma(find(~isnan(delta_sigma))),delta_sigma_boot(find(~isnan(delta_sigma_boot))))
p_deltaLSFS=ranksum(delta_lsfs(find(~isnan(delta_lsfs))),delta_lsfs_boot(find(~isnan(delta_lsfs_boot))))
p_deltaHSFS=ranksum(delta_hsfs(find(~isnan(delta_hsfs))),delta_hsfs_boot(find(~isnan(delta_hsfs_boot))))
p_signIndTF=ranksum(sign_index_tf(find(~isnan(sign_index_tf))),sign_index_tf_boot(find(~isnan(sign_index_tf_boot))))
p_signIndSigma=ranksum(sign_index_sigma(find(~isnan(sign_index_sigma))),sign_index_sigma_boot(find(~isnan(sign_index_sigma_boot))))
p_signIndLSFS=ranksum(sign_index_lsfs(find(~isnan(sign_index_lsfs))),sign_index_lsfs_boot(find(~isnan(sign_index_lsfs_boot))))
p_signIndHSFS=ranksum(sign_index_hsfs(find(~isnan(sign_index_hsfs))),sign_index_hsfs_boot(find(~isnan(sign_index_hsfs_boot))))
[r_TFdyn,p_TFdyn]=corr(mean_pref_TF',mean_pref_TF_boot','type','spearman','rows','pairwise')
[r_SigmaDyn,p_SigmaDyn]=corr(mean_sigma',mean_sigma_boot','type','spearman','rows','pairwise')
[r_LSFSdyn,p_LSFSdyn]=corr(mean_lsfs',mean_lsfs_boot','type','spearman','rows','pairwise')
[r_HSFSdyn,p_HSFSdyn]=corr(mean_hsfs',mean_hsfs_boot','type','spearman','rows','pairwise')

r_TFdyn_iterations=nan(1000,1); p_TFdyn_iterations=nan(1000,1);
r_SigmaDyn_iterations=nan(1000,1); p_SigmaDyn_iterations=nan(1000,1);
r_LSFSdyn_iterations=nan(1000,1); p_LSFSdyn_iterations=nan(1000,1);
r_HSFSdyn_iterations=nan(1000,1); p_HSFSdyn_iterations=nan(1000,1);

for i=1:1000
    [r_TFdyn_iterations(i),p_TFdyn_iterations(i)]=corr(mean_pref_TF',pref_TF_boot(i,:)','type','spearman','rows','pairwise');
    [r_SigmaDyn_iterations(i),p_SigmaDyn_iterations(i)]=corr(mean_sigma',sigma_boot(i,:)','type','spearman','rows','pairwise');
    [r_LSFSdyn_iterations(i),p_LSFSdyn_iterations(i)]=corr(mean_lsfs',lsfs_boot(i,:)','type','spearman','rows','pairwise');
    [r_HSFSdyn_iterations(i),p_HSFSdyn_iterations(i)]=corr(mean_hsfs',hsfs_boot(i,:)','type','spearman','rows','pairwise');
end
clear i
mean_r_iterations_sf=nanmean(r_TFdyn_iterations)
mean_p_iterations_sf=nanmean(p_TFdyn_iterations)
nSiginif_iterations_sf=numel(find(p_TFdyn_iterations<0.05))
mean_r_iterations_sigma=nanmean(r_SigmaDyn_iterations)
mean_p_iterations_sigma=nanmean(p_SigmaDyn_iterations)
nSiginif_iterations_sigma=numel(find(p_SigmaDyn_iterations<0.05))
mean_r_iterations_lsfs=nanmean(r_LSFSdyn_iterations)
mean_p_iterations_lsfs=nanmean(p_LSFSdyn_iterations)
nSiginif_iterations_lsfs=numel(find(p_LSFSdyn_iterations<0.05))
mean_r_iterations_hsfs=nanmean(r_HSFSdyn_iterations)
mean_p_iterations_hsfs=nanmean(p_HSFSdyn_iterations)
nSiginif_iterations_hsfs=numel(find(p_HSFSdyn_iterations<0.05))