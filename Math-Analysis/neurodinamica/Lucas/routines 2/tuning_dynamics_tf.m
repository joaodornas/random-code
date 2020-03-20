% load ('latForEachFreq','complex_cells','TFs_i_atBestSF','TFvalues_octave','lat_TFs_atBestSF','TF_columnLabels_octv');
% ncells=numel(complex_cells);
% analysis_period=800;
% window_size=10;
% baseline_duration=1000;
% analysis_startTime=201;
% 
% pref_TFs_dyn=nan(ncells,analysis_period);
% BWs_dyn=nan(ncells,analysis_period);
% sigmas_dyn=nan(ncells,analysis_period);
% LSFVs_dyn=nan(ncells,analysis_period);
% LSFSs_dyn=nan(ncells,analysis_period);
% HSFSs_dyn=nan(ncells,analysis_period);
% response_amp_dyn=nan(ncells,analysis_period);
% fit_pvalue=nan(ncells,analysis_period);
% rsquare=nan(ncells,analysis_period);
% x_axis=(0.05:0.05:64)';

priebe_1D=fittype('A.*((exp((-((log2(x)-log2(ps))^2)./(2.*(tuningWidth+skew.*(log2(x)-log2(ps))).^2))))-(exp(-1./((skew).^2))))');

for i=82:ncells%[1:4,6:12,14:ncells]
    
    tic
    
    disp(['fitting and calculating BWs & LSFSs of cell #' num2str(i)])
    
    fitted_curves_i=nan(numel(x_axis),analysis_period);

    %get SF conditions and values for cell i
    eval(['TFs_i=TFs_i_atBestSF.complexCell' num2str(i) ';']);
    eval(['TFvalues_octv=TFvalues_octave.complexCell' num2str(i) ';']);

    TFvalues_lin=2.^TFvalues_octv;
    eval(['TFvalues_linScale.complexCell' num2str(i) '=TFvalues_lin;']);
    
    %load and shape psths in matrices where each row is a condition
    load(['psth_' char(complex_cells(i)) '_convolved_sigma3.mat'],'psth');

    TF_psths=zeros(numel(TFs_i),numel(psth.condition1));
    for k=1:numel(TFs_i)
        TF_psths(k,:)=eval(['psth.condition' num2str(TFs_i(k))]);
        mean_baseline=mean(TF_psths(k,1:baseline_duration));
        TF_psths(k,:)=TF_psths(k,:)-mean_baseline;
    end

    clear k psth TFs_i
    
    %latencies for each condition, if nan, then lat=55 (mean latency)
    tf_column_i=find(TF_columnLabels_octv==TFvalues_octv(1));
    lat=lat_TFs_atBestSF(i,tf_column_i:tf_column_i+numel(TFvalues_octv)-1)';
    lat(find(isnan(lat)))=55;
    tf_responses=zeros(size(lat));
    tf_responses_cell=nan(numel(TFvalues_octv),analysis_period);
    tf_response_sem=zeros(size(lat));
    tf_response_sem_cell=nan(numel(TFvalues_octv),analysis_period);

    for k=476:analysis_period+analysis_startTime-1 %slide windows, fit and extract parameters
        
        col_i=k-analysis_startTime+1;
        
        for j=1:numel(lat)
        tf_responses(j)=mean(TF_psths(j,baseline_duration+k+lat(j)-(window_size/2):baseline_duration+k+lat(j)+(window_size/2)));
        tf_response_sem(j)=std(TF_psths(j,baseline_duration+k+lat(j)-(window_size/2):baseline_duration+k+lat(j)+(window_size/2)))/sqrt(window_size);
        end
        
        if max(tf_responses)<0.1
            
            eval(['fittedmodels.cell' num2str(i) '.window' num2str(k) '=nan;']);
            eval(['gofs.cell' num2str(i) '.window' num2str(k) '=nan;']);
            eval(['outputs.cell' num2str(i) '.window' num2str(k) '=nan;']);
            
        else
            
        initial_values=[max(tf_responses) TFvalues_lin(find(tf_responses==max(tf_responses))) 0 1.5];
        lower_limits=[0.1,0.05,-1,0.1]; 
        upper_limits=[max(tf_responses)+20,10,1,8];
        [fittedmodel,gof,output]=fit(TFvalues_lin',tf_responses,priebe_1D,'maxIter',1000,'startpoint',initial_values,'lower',lower_limits,'upper',upper_limits);

        [fit_pvalue(i,col_i)]=fitp(gof,output); %evaluate goodness of fit with f stats
        rsquare(i,col_i)=gof.rsquare;

        if rsquare(i,col_i)>0.7 %accept only fits with r2>0.7
            
            pref_TFs_dyn(i,col_i)=fittedmodel.ps;
            sigmas_dyn(i,col_i)=fittedmodel.tuningWidth;
            response_amp_dyn(i,col_i)=fittedmodel.A;
            fitted_curves_i(:,col_i)=feval(fittedmodel,x_axis);
            BWs_dyn(i,col_i)=fullWidthAtHalfHeight(x_axis,fitted_curves_i(:,col_i),'octave');
            
            tf=TFvalues_lin(1); A=fittedmodel.A; tuningWidth=fittedmodel.tuningWidth; bestTF=fittedmodel.ps; skew=fittedmodel.skew;
            lowest_TF_resp=A*((exp((-((log2(tf)-log2(bestTF)).^2)./(2*(tuningWidth+skew.*(log2(tf)-log2(bestTF))).^2))))-(exp(-1/((skew)^2))));
            tf=bestTF;
            best_TF_resp=A*((exp((-((log2(tf)-log2(bestTF)).^2)./(2*(tuningWidth+skew.*(log2(tf)-log2(bestTF))).^2))))-(exp(-1/((skew)^2))));
            tf=TFvalues_lin(numel(TFvalues_lin));
            highest_TF_resp=A*((exp((-((log2(tf)-log2(bestTF)).^2)./(2*(tuningWidth+skew.*(log2(tf)-log2(bestTF))).^2))))-(exp(-1/((skew)^2))));
            HSFSs_dyn(i,col_i)=highest_TF_resp/best_TF_resp;
            LSFSs_dyn(i,col_i)=lowest_TF_resp/best_TF_resp;
            clear tf lowest_TF_resp best_TF_resp A tuningWidth bestTF skew

        end

        tf_responses_cell(:,col_i)=tf_responses;
        tf_response_sem_cell(:,col_i)=tf_response_sem;
        eval(['fittedmodels.cell' num2str(i) '.window' num2str(k) '=fittedmodel;']);
        eval(['gofs.cell' num2str(i) '.window' num2str(k) '=gof;']);
        eval(['outputs.cell' num2str(i) '.window' num2str(k) '=output;']);

        end
               
    end

    eval(['tf_responses_dyn.cell' num2str(i) '=tf_responses_cell;']);
    eval(['fitted_curves.cell' num2str(i) '=fitted_curves_i;']);
    clear k j col_i TF_psths TFvalues_lin initial_values tf_responses tf_response_sem TFvalues_octv fittedmodel gof output tf_responses_cell fitted_curves_i mean_baseline upper_limits lower_limits lat tf_column_i
    
    %save(['/Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/TF_tuningDynamics_windowSize' num2str(window_size) '_800.mat']);
    
    disp(['elapsed time for cell ' num2str(i) ' = ' num2str(toc) ' seconds'])
end

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

save(['/Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/TF_tuningDynamics_windowSize' num2str(window_size) '_800.mat']);

% for i=[1:4,6:12,14:ncells]
%     tic
%     SFvalues_lin=eval(['TFvalues_linScale.complexCell' num2str(i)]);
%     if i==1
%         for k=68:analysis_period
%             disp(['computing LSFVs of cell #' num2str(i) ', time window ' num2str(k) '...'])
%             if ~isnan(pref_SFs_dyn(i,k))
%                 fittedmodel=eval(['fittedmodels.cell' num2str(i) '.window' num2str(k)]);
%                 try
%                     [LSFVs_dyn(i,k)]=calculate_lsfv(fittedmodel,TFvalues_lin(1),'off');
%                 catch
%                     LSFVs_dyn(i,k)=nan;
%                 end
%             end
%         end
%     else
%         for k=1:analysis_period
%             disp(['computing LSFVs of cell #' num2str(i) ', time window ' num2str(k) '...'])
%             if ~isnan(pref_SFs_dyn(i,k))
%                 fittedmodel=eval(['fittedmodels.cell' num2str(i) '.window' num2str(k)]);
%                 try
%                     [LSFVs_dyn(i,k)]=calculate_lsfv(fittedmodel,TFvalues_lin(1),'off');
%                 catch
%                     LSFVs_dyn(i,k)=nan;
%                 end
%             end
%         end
%     end
%     clear TFvalues_lin fittedmodel
%     save(['/Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/TF_tuningDynamics_windowSize' num2str(window_size) '.mat']);
%     disp(['elapsed time for cell ' num2str(i) ' = ' num2str(fixdec(toc/60,1)) ' minutes'])
% end
% 
% non_empty_cols_lsfv=nan(1,analysis_period);
% 
% for k=1:analysis_period
%     non_empty_cols_lsfv(i)=numel(find(~isnan(LSFVs_dyn(:,k))));
% end
% 
% mean_lsfv_dyn=nanmean(LSFVs_dyn);
% sem_lsfv_dyn=nanstd(LSFVs_dyn)./sqrt(non_empty_cols_lsfv);
% 
% save(['/Users/lucaspinto/Documents/Lab/ProjectBooks/CV&Dyn-Book/Analyses/SF_tuningDynamics_windowSize' num2str(window_size) '.mat']);