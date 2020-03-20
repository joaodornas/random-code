% response_amp_dyn=nan(92,analysis_period);
% HSFSs_dyn=nan(92,analysis_period);
% 
% for i=[1:4,6:12,14:92]
%     for k=1:analysis_period
%         SFvalues_lin=eval(['SFvalues_linScale.complexCell' num2str(i)]);
%         try
%             fittedmodel=eval(['fittedmodels.cell' num2str(i) '.window' num2str(k) ';']);
%             response_amp_dyn(i,k)=fittedmodel.A;
%             A=fittedmodel.A; tuningWidth=fittedmodel.tuningWidth; bestSF=fittedmodel.ps; skew=fittedmodel.skew;
%             sf=bestSF;
%             best_SF_resp=A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))));
%             sf=SFvalues_lin(numel(SFvalues_lin));
%             highest_SF_resp=A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))));
%             HSFSs_dyn(i,k)=highest_SF_resp/best_SF_resp;
%         catch
%             response_amp_dyn(i,k)=nan;
%             HSFSs_dyn(i,k)=nan;
%         end
%     end
% end

ncells=92;
analysis_period=1000;
savepath='/Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/tuningDynamics/';

resp_800=load('SF_tuningDynamics_windowSize10_800.mat','response_amp_dyn');
response_amp_dyn=[response_amp_dyn resp_800.response_amp_dyn];
clear resp_800
a_800=load('SF_tuningDynamics_windowSize10_800.mat','sigmas_dyn');
sigmas_dyn=[sigmas_dyn a_800.sigmas_dyn];
clear a_800
a_800=load('SF_tuningDynamics_windowSize10_800.mat','HSFSs_dyn');
HSFSs_dyn=[HSFSs_dyn a_800.HSFSs_dyn];
clear a_800
a_800=load('SF_tuningDynamics_windowSize10_800.mat','LSFSs_dyn');
LSFSs_dyn=[LSFSs_dyn a_800.LSFSs_dyn];
clear a_800
a_800=load('SF_tuningDynamics_windowSize10_800.mat','LSFVs_dyn');
LSFVs_dyn=[LSFVs_dyn a_800.LSFVs_dyn];
clear a_800
a_800=load('SF_tuningDynamics_windowSize10_800.mat','pref_SFs_dyn');
pref_SFs_dyn=[pref_SFs_dyn a_800.pref_SFs_dyn];
clear a_800
a_800=load('SF_tuningDynamics_windowSize10_800.mat','BWs_dyn');
BWs_dyn=[BWs_dyn a_800.BWs_dyn];
clear a_800

save([savepath 'SF_tuningDynamics_analysis_windowSize' num2str(window_size) '_anaPeriod' num2str(analysis_period) '.mat']);