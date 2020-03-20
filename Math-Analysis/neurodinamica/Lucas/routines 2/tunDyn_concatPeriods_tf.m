% %response_amp_dyn=nan(92,analysis_period);
% HSFSs_dyn=nan(92,analysis_period);
% 
% for i=[1:4,6:12,14:92]
%     for k=1:analysis_period
%         try
%             TFvalues_lin=eval(['TFvalues_linScale.complexCell' num2str(i)]);
%             fittedmodel=eval(['fittedmodels.cell' num2str(i) '.window' num2str(k) ';']);
%             %response_amp_dyn(i,k)=fittedmodel.A;
%             A=fittedmodel.A; tuningWidth=fittedmodel.tuningWidth; bestSF=fittedmodel.ps; skew=fittedmodel.skew;
%             sf=bestSF;
%             best_SF_resp=A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))));
%             sf=TFvalues_lin(numel(TFvalues_lin));
%             highest_SF_resp=A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))));
%             HSFSs_dyn(i,k)=highest_SF_resp/best_SF_resp;
%             clear A tuningWidth bestSF skew sf highest_SF_resp best_SF_resp TFvalues_lin fittedmodel
%             %response_amp_dyn(i,k)=nan;
%         catch
%             HSFSs_dyn(i,k)=nan;
%         end
%     end
% end
% clear i k

analysis_period=1000;
savepath='/Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/tuningDynamics/';

resp_800=load('TF_tuningDynamics_windowSize10_800.mat','response_amp_dyn');
response_amp_dyn=[response_amp_dyn resp_800.response_amp_dyn];
clear resp_800
a_800=load('TF_tuningDynamics_windowSize10_800.mat','sigmas_dyn');
sigmas_dyn=[sigmas_dyn a_800.sigmas_dyn];
clear a_800
a_800=load('TF_tuningDynamics_windowSize10_800.mat','HSFSs_dyn');
HSFSs_dyn=[HSFSs_dyn a_800.HSFSs_dyn];
clear a_800
a_800=load('TF_tuningDynamics_windowSize10_800.mat','LSFSs_dyn');
LSFSs_dyn=[LSFSs_dyn a_800.LSFSs_dyn];
clear a_800
a_800=load('TF_tuningDynamics_windowSize10_800.mat','LSFVs_dyn');
LSFVs_dyn=[LSFVs_dyn a_800.LSFVs_dyn];
clear a_800
a_800=load('TF_tuningDynamics_windowSize10_800.mat','pref_TFs_dyn');
pref_TFs_dyn=[pref_TFs_dyn a_800.pref_TFs_dyn];
clear a_800
a_800=load('SF_tuningDynamics_windowSize10_800.mat','BWs_dyn');
BWs_dyn=[BWs_dyn a_800.BWs_dyn];
clear a_800

save([savepath 'TF_tuningDynamics_analysis_windowSize' num2str(window_size) '_anaPeriod' num2str(analysis_period) '.mat']);